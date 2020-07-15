//! Constructor
// mpm::XMPMExplicit<Tdim>::XMPMExplicit(const std::shared_ptr<IO>& io)
//     : mpm::MPMBase<Tdim>(io) {
template <unsigned Tdim>
mpm::XMPMExplicit<Tdim>::XMPMExplicit(const std::shared_ptr<IO>& io)
    : mpm::MPMBase<Tdim>(io) {
  //! Logger
  console_ = spdlog::get("XMPMExplicit");
}

//! MPM Explicit pressure smoothing
template <unsigned Tdim>
void mpm::XMPMExplicit<Tdim>::pressure_smoothing(unsigned phase) {
  // Assign pressure to nodes
  mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::map_pressure_to_nodes, std::placeholders::_1));

#ifdef USE_MPI
  int mpi_size = 1;

  // Get number of MPI ranks
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  // Run if there is more than a single MPI task
  if (mpi_size > 1) {
    // MPI all reduce nodal pressure
    mesh_->template nodal_halo_exchange<double, 1>(
        std::bind(&mpm::NodeBase<Tdim>::pressure, std::placeholders::_1, phase),
        std::bind(&mpm::NodeBase<Tdim>::assign_pressure, std::placeholders::_1,
                  phase, std::placeholders::_2));
  }
#endif

  // Smooth pressure over particles
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::compute_pressure_smoothing,
                std::placeholders::_1));
}

//! MPM Explicit compute stress strain
template <unsigned Tdim>
void mpm::XMPMExplicit<Tdim>::compute_stress_strain(unsigned phase) {
  // Iterate over each particle to calculate strain
  mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::compute_strain, std::placeholders::_1, dt_));

  // Iterate over each particle to update particle volume
  mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::update_volume, std::placeholders::_1));

  // Pressure smoothing
  if (pressure_smoothing_) this->pressure_smoothing(phase);

  // Iterate over each particle to compute stress
  mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::compute_stress, std::placeholders::_1));
}

// Initialise discontinuities
template <unsigned Tdim>
bool mpm::XMPMExplicit<Tdim>::initialise_discontinuities() {
  bool status = true;
  try {
      // Get discontinuities data
    auto json_discontinuities = io_->json_object("discontinuity");

    for (const auto discontinuity_props : json_discontinuities) {
      // Get discontinuity type
      const std::string discontunity_type =
          discontinuity_props["type"].template get<std::string>();

      // Get discontinuity id
      auto discontinuity_id = discontinuity_props["id"].template get<unsigned>();

        // Get discontinuity  input type
      auto io_type = discontinuity_props["io_type"].template get<std::string>();

      // discontinuity file
      std::string discontinuity_file =
          io_->file_name(discontinuity_props["file"].template get<std::string>());

      auto discontinuity_frictional_coef = discontinuity_props["frictional_coefficient"].template get<double>();

          // Create a mesh reader
      auto discontunity_io = Factory<mpm::IOMesh<Tdim>>::instance()->create(io_type);
      
      // Create a new discontinuity surface from JSON object
      auto discontinuity =
          Factory<mpm::DiscontinuityBase<Tdim>>::instance()
              ->create(discontunity_type);

      bool status = 
        discontinuity->initialize(discontunity_io->read_mesh_nodes(discontinuity_file),discontunity_io->read_mesh_cells(discontinuity_file));

      discontinuity->set_frictional_coef(discontinuity_frictional_coef);
       // Create points from file
  
      // Add discontinuity to list
      auto result = discontinuities_.insert(std::make_pair(discontinuity_id, discontinuity));

      // If insert discontinuity failed
      if (!result.second) {
        status = false;
        throw std::runtime_error(
            "New discontinuity cannot be added, insertion failed");
      }
    }
  } catch (std::exception& exception) {
    console_->error("#{}: Reading discontinuities: {}", __LINE__, exception.what());
    status = false;
  }
  return status;
}
//! MPM Explicit solver
template <unsigned Tdim>
bool mpm::XMPMExplicit<Tdim>::solve() {
  bool status = true;

  console_->info("MPM analysis type {}", io_->analysis_type());

  // Initialise MPI rank and size
  int mpi_rank = 0;
  int mpi_size = 1;

#ifdef USE_MPI
  // Get MPI rank
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  // Get number of MPI ranks
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

  // Phase
  const unsigned phase = 0;

  // Test if checkpoint resume is needed
  bool resume = false;
  if (analysis_.find("resume") != analysis_.end())
    resume = analysis_["resume"]["resume"].template get<bool>();

  // Pressure smoothing
  if (analysis_.find("pressure_smoothing") != analysis_.end())
    pressure_smoothing_ =
        analysis_.at("pressure_smoothing").template get<bool>();

  // Interface
  if (analysis_.find("interface") != analysis_.end())
    interface_ = analysis_.at("interface").template get<bool>();

  // Discontinuity
  if (!io_->json_object("discontinuity").empty())
    discontinuity_ = true;
  // Initialise material
  bool mat_status = this->initialise_materials();
  if (!mat_status) {
    status = false;
    throw std::runtime_error("Initialisation of materials failed");
  }

  // Initialise mesh
  bool mesh_status = this->initialise_mesh();
  if (!mesh_status) {
    status = false;
    throw std::runtime_error("Initialisation of mesh failed");
  }

  // Initialise particles
  bool particle_status = this->initialise_particles();
  if (!particle_status) {
    status = false;
    throw std::runtime_error("Initialisation of particles failed");
  }

  // Initialise loading conditions
  bool loading_status = this->initialise_loads();
  if (!loading_status) {
    status = false;
    throw std::runtime_error("Initialisation of loading failed");
  }

  // Initialise discontinuity
  if (discontinuity_)
  {
    bool discontinuity_status = this->initialise_discontinuities();
    if (!discontinuity_status) {
      status = false;
      throw std::runtime_error("Initialisation of discontinuities failed");
    }

    bool initialise_lsm_status = this->initialise_lsm();
    if (!discontinuity_status) {
      status = false;
      throw std::runtime_error("Initialisation of level set method failed");
    }    
  }

  // Compute mass
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::compute_mass, std::placeholders::_1));

  // Domain decompose
  bool initial_step = true;
  this->mpi_domain_decompose(initial_step);

  // Check point resume
  if (resume) this->checkpoint_resume();

  auto solver_begin = std::chrono::steady_clock::now();
  // Main loop
  for (; step_ < nsteps_; ++step_) {

    if (mpi_rank == 0) console_->info("Step: {} of {}.\n", step_, nsteps_);

    // Inject particles
    mesh_->inject_particles(this->step_ * this->dt_);

    // Create a TBB task group
    tbb::task_group task_group;

    // Spawn a task for initialising nodes and cells
    task_group.run([&] {
      // Initialise nodes
      mesh_->iterate_over_nodes(
          std::bind(&mpm::NodeBase<Tdim>::initialise, std::placeholders::_1));

      mesh_->iterate_over_cells(
          std::bind(&mpm::Cell<Tdim>::activate_nodes, std::placeholders::_1));
    });

    // Spawn a task for particles
    task_group.run([&] {
      // Iterate over each particle to compute shapefn
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::compute_shapefn, std::placeholders::_1));
    });

    task_group.wait();

    // Assign material ids to node
    if (interface_)
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::append_material_id_to_nodes,
                    std::placeholders::_1));

    // Assign mass and momentum to nodes
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_mass_momentum_to_nodes,
                  std::placeholders::_1));

#ifdef USE_MPI
    // Run if there is more than a single MPI task
    if (mpi_size > 1) {
      // MPI all reduce nodal mass
      mesh_->template nodal_halo_exchange<double, 1>(
          std::bind(&mpm::NodeBase<Tdim>::mass, std::placeholders::_1, phase),
          std::bind(&mpm::NodeBase<Tdim>::update_mass, std::placeholders::_1,
                    false, phase, std::placeholders::_2));
      // MPI all reduce nodal momentum
      mesh_->template nodal_halo_exchange<Eigen::Matrix<double, Tdim, 1>, Tdim>(
          std::bind(&mpm::NodeBase<Tdim>::momentum, std::placeholders::_1,
                    phase),
          std::bind(&mpm::NodeBase<Tdim>::update_momentum,
                    std::placeholders::_1, false, phase,
                    std::placeholders::_2));
    }
#endif

    // Compute nodal velocity
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_velocity,
                  std::placeholders::_1),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Update stress first
    if (this->stress_update_ == mpm::StressUpdate::USF)
      this->compute_stress_strain(phase);

    // Spawn a task for external force
    task_group.run([&] {
      // Iterate over each particle to compute nodal body force
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_body_force,
                    std::placeholders::_1, this->gravity_));

      // Apply particle traction and map to nodes
      mesh_->apply_traction_on_particles(this->step_ * this->dt_);

      // Iterate over each node to add concentrated node force to external force
      if (set_node_concentrated_force_)
        mesh_->iterate_over_nodes(
            std::bind(&mpm::NodeBase<Tdim>::apply_concentrated_force,
                      std::placeholders::_1, phase, (this->step_ * this->dt_)));
    });

    // Spawn a task for internal force
    task_group.run([&] {
      // Iterate over each particle to compute nodal internal force
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::map_internal_force, std::placeholders::_1));
    });
    task_group.wait();

#ifdef USE_MPI
    // Run if there is more than a single MPI task
    if (mpi_size > 1) {
      // MPI all reduce external force
      mesh_->template nodal_halo_exchange<Eigen::Matrix<double, Tdim, 1>, Tdim>(
          std::bind(&mpm::NodeBase<Tdim>::external_force, std::placeholders::_1,
                    phase),
          std::bind(&mpm::NodeBase<Tdim>::update_external_force,
                    std::placeholders::_1, false, phase,
                    std::placeholders::_2));
      // MPI all reduce internal force
      mesh_->template nodal_halo_exchange<Eigen::Matrix<double, Tdim, 1>, Tdim>(
          std::bind(&mpm::NodeBase<Tdim>::internal_force, std::placeholders::_1,
                    phase),
          std::bind(&mpm::NodeBase<Tdim>::update_internal_force,
                    std::placeholders::_1, false, phase,
                    std::placeholders::_2));
    }
#endif

    //intergrate momentum Iterate over
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::intergrate_momentum,
                  std::placeholders::_1, phase, this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Iterate over each particle to compute updated position
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::compute_updated_position,
                  std::placeholders::_1, this->dt_, this->velocity_update_));

    // Apply particle velocity constraints
    mesh_->apply_particle_velocity_constraints();

    // Update Stress Last
    if (this->stress_update_ == mpm::StressUpdate::USL)
      this->compute_stress_strain(phase);

    // Locate particles
    auto unlocatable_particles = mesh_->locate_particles_mesh();

    if (!unlocatable_particles.empty() && this->locate_particles_)
      throw std::runtime_error("Particle outside the mesh domain");
    // If unable to locate particles remove particles
    if (!unlocatable_particles.empty() && !this->locate_particles_)
      for (const auto& remove_particle : unlocatable_particles)
        mesh_->remove_particle(remove_particle);

#ifdef USE_MPI
#ifdef USE_GRAPH_PARTITIONING
    // Run load balancer at a specified frequency
    if (step_ % nload_balance_steps_ == 0 && step_ != 0)
      this->mpi_domain_decompose(false);

    mesh_->transfer_nonrank_particles();
#endif
#endif

    if (step_ % output_steps_ == 0) {
      // HDF5 outputs
      this->write_hdf5(this->step_, this->nsteps_);
#ifdef USE_VTK
      // VTK outputs
      this->write_vtk(this->step_, this->nsteps_);
#endif
#ifdef USE_PARTIO
      // Partio outputs
      this->write_partio(this->step_, this->nsteps_);
#endif
    }
  }
  auto solver_end = std::chrono::steady_clock::now();
  console_->info(
      "Rank {}, Explicit {} solver duration: {} ms", mpi_rank,
      (this->stress_update_ == mpm::StressUpdate::USL ? "USL" : "USF"),
      std::chrono::duration_cast<std::chrono::milliseconds>(solver_end -
                                                            solver_begin)
          .count());

  return status;
}


// Initialise particles
template <unsigned Tdim>
bool mpm::XMPMExplicit<Tdim>::initialise_lsm() {

  bool status = true;

  try {

    for(mpm::Index i=0; i<ndiscontinuities();++i){

      std::vector<double> phi_list(mesh_->nparticles());

      discontinuities_[i]->compute_phi(mesh_->particle_coordinates(),  phi_list);

      mesh_->set_particle_phi(phi_list);      
    }

  } catch (std::exception& exception) {
    console_->error("#{}: XMPM initialise particles LSM: {}", __LINE__,
                    exception.what());
    status = false;
  }

  if (!status) throw std::runtime_error("Initialisation of particles LSM values failed");

  return status;
}