//! Constructor with id, coordinates and dof
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
mpm::NodeXMPM<Tdim, Tdof, Tnphases>::NodeXMPM(
    Index id, const Eigen::Matrix<double, Tdim, 1>& coord)
    : NodeBase<Tdim>(id, coord) {
  // Check if the dimension is between 1 & 3
  static_assert((Tdim >= 1 && Tdim <= 3), "Invalid global dimension");
  id_ = id;
  coordinates_ = coord;
  dof_ = Tdof;

  //! Logger
  std::string logger =
      "node" + std::to_string(Tdim) + "dxmpm::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);

  // Clear any velocity constraints
  velocity_constraints_.clear();
  concentrated_force_.setZero();
  this->initialise();
}

//! Initialise nodal properties
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::NodeXMPM<Tdim, Tdof, Tnphases>::initialise() noexcept {
  mass_.setZero();
  volume_.setZero();
  external_force_.setZero();
  internal_force_.setZero();
  pressure_.setZero();
  velocity_.setZero();
  momentum_.setZero();
  acceleration_.setZero();
  status_ = false;
  material_ids_.clear();
  mass_h_.setZero();
  momentum_h_.setZero();
  internal_force_h_.setZero();
  external_force_h_.setZero();
  enrich_h_ = true;
  direction_discontinuity_.setZero();
  direction_discontinuity_(0,0) = 0.4472136;
  direction_discontinuity_(1,0)= 0.0;
  direction_discontinuity_(2,0) = 0.8944272;
  frictional_coef = 0.75;
}

//! Initialise shared pointer to nodal properties pool
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::NodeXMPM<Tdim, Tdof, Tnphases>::initialise_property_handle(
    unsigned prop_id,
    std::shared_ptr<mpm::NodalProperties> property_handle) noexcept {
  // the property handle and the property id is set in the node
  this->property_handle_ = property_handle;
  this->prop_id_ = prop_id;
}

//! Update mass at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::NodeXMPM<Tdim, Tdof, Tnphases>::update_mass(bool update, unsigned phase,
                                                  double mass, double phi) noexcept {
  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign mass
  std::lock_guard<std::mutex> guard(node_mutex_);
  mass_(phase) = (mass_(phase) * factor) + mass;
  
  if(enrich_h_)
    mass_h_(phase) = (mass_h_(phase) * factor) + mass * sgn(phi);
}

//! Update volume at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::NodeXMPM<Tdim, Tdof, Tnphases>::update_volume(bool update, unsigned phase,
                                                    double volume) noexcept {
  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign volume
  std::lock_guard<std::mutex> guard(node_mutex_);
  volume_(phase) = volume_(phase) * factor + volume;
}

// Assign concentrated force to the node
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::NodeXMPM<Tdim, Tdof, Tnphases>::assign_concentrated_force(
    unsigned phase, unsigned direction, double concentrated_force,
    const std::shared_ptr<FunctionBase>& function) {
  bool status = false;
  try {
    if (phase >= Tnphases || direction >= Tdim) {
      throw std::runtime_error(
          "Cannot assign nodal concentrated forcey: Direction / phase is "
          "invalid");
    }
    // Assign concentrated force
    concentrated_force_(direction, phase) = concentrated_force;
    status = true;
    this->force_function_ = function;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Apply concentrated force to the node
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::NodeXMPM<Tdim, Tdof, Tnphases>::apply_concentrated_force(
    unsigned phase, double current_time) {
  const double scalar =
      (force_function_ != nullptr) ? force_function_->value(current_time) : 1.0;
  this->update_external_force(true, phase,
                              scalar * concentrated_force_.col(phase));
}

//! Update external force (body force / traction force)
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::NodeXMPM<Tdim, Tdof, Tnphases>::update_external_force(
    bool update, unsigned phase,
    const Eigen::Matrix<double, Tdim, 1>& force) noexcept {
  // Assert
  assert(phase < Tnphases);

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign external force
  std::lock_guard<std::mutex> guard(node_mutex_);
  external_force_.col(phase) = external_force_.col(phase) * factor + force;
}

//! Update external force (body force / traction force)
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::NodeXMPM<Tdim, Tdof, Tnphases>::update_external_force(
    bool update, unsigned phase,
    const Eigen::Matrix<double, Tdim, 1>& force, double phi) noexcept {
  // Assert
  assert(phase < Tnphases);

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign external force
  std::lock_guard<std::mutex> guard(node_mutex_);
  external_force_.col(phase) = external_force_.col(phase) * factor + force;

  if(enrich_h_)
    external_force_h_.col(phase) = external_force_h_.col(phase) * factor + sgn(phi)*force;

}

//! Update internal force (body force / traction force)
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::NodeXMPM<Tdim, Tdof, Tnphases>::update_internal_force(
    bool update, unsigned phase,
    const Eigen::Matrix<double, Tdim, 1>& force, double phi) noexcept {
  // Assert
  assert(phase < Tnphases);

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign internal force
  std::lock_guard<std::mutex> guard(node_mutex_);
  internal_force_.col(phase) = internal_force_.col(phase) * factor + force;
  if(enrich_h_)
    internal_force_h_.col(phase) = internal_force_h_.col(phase) * factor + force*sgn(phi);
}

//! Assign nodal momentum
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::NodeXMPM<Tdim, Tdof, Tnphases>::update_momentum(
    bool update, unsigned phase,
    const Eigen::Matrix<double, Tdim, 1>& momentum, double phi) noexcept {
  // Assert
  assert(phase < Tnphases);

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign momentum
  std::lock_guard<std::mutex> guard(node_mutex_);
  momentum_.col(phase) = momentum_.col(phase) * factor + momentum;

  if(enrich_h_)
    momentum_h_.col(phase) = momentum_h_.col(phase) * factor + momentum * sgn(phi);

}

//! Update pressure at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::NodeXMPM<Tdim, Tdof, Tnphases>::update_mass_pressure(
    unsigned phase, double mass_pressure) noexcept {
  // Assert
  assert(phase < Tnphases);

  const double tolerance = 1.E-16;
  // Compute pressure from mass*pressure
  if (mass_(phase) > tolerance) {
    std::lock_guard<std::mutex> guard(node_mutex_);
    pressure_(phase) += mass_pressure / mass_(phase);
  }
}

//! Assign pressure at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::NodeXMPM<Tdim, Tdof, Tnphases>::assign_pressure(unsigned phase,
                                                      double pressure) {
  const double tolerance = 1.E-16;

  // Compute pressure from mass*pressure
  std::lock_guard<std::mutex> guard(node_mutex_);
  pressure_(phase) = pressure;
}

//! Compute velocity from momentum
//! velocity = momentum / mass
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::NodeXMPM<Tdim, Tdof, Tnphases>::compute_velocity() {
  const double tolerance = 1.E-16;
  for (unsigned phase = 0; phase < Tnphases; ++phase) {
    if (mass_(phase) > tolerance) {
      velocity_.col(phase) = momentum_.col(phase) / mass_(phase);

      // Check to see if value is below threshold
      for (unsigned i = 0; i < velocity_.rows(); ++i)
        if (std::abs(velocity_.col(phase)(i)) < 1.E-15)
          velocity_.col(phase)(i) = 0.;
    }
  }

  // Apply velocity constraints, which also sets acceleration to 0,
  // when velocity is set.
  this->apply_velocity_constraints();
}

//! Update nodal acceleration
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::NodeXMPM<Tdim, Tdof, Tnphases>::update_acceleration(
    bool update, unsigned phase,
    const Eigen::Matrix<double, Tdim, 1>& acceleration) noexcept {
  assert(phase < Tnphases);

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  //! Update/assign acceleration
  std::lock_guard<std::mutex> guard(node_mutex_);
  acceleration_.col(phase) = acceleration_.col(phase) * factor + acceleration;
}

//! Compute acceleration and velocity
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::NodeXMPM<Tdim, Tdof, Tnphases>::compute_acceleration_velocity(
    unsigned phase, double dt) noexcept {
  bool status = false;
  const double tolerance = 1.0E-15;
  if (mass_(phase) > tolerance) {
    // acceleration = (unbalaced force / mass)
    this->acceleration_.col(phase) =
        (this->external_force_.col(phase) + this->internal_force_.col(phase)) /
        this->mass_(phase);

    // Apply friction constraints
    this->apply_friction_constraints(dt);

    // Velocity += acceleration * dt
    this->velocity_.col(phase) += this->acceleration_.col(phase) * dt;
    // Apply velocity constraints, which also sets acceleration to 0,
    // when velocity is set.
    this->apply_velocity_constraints();

    // Set a threshold
    for (unsigned i = 0; i < Tdim; ++i)
      if (std::abs(velocity_.col(phase)(i)) < tolerance)
        velocity_.col(phase)(i) = 0.;
    for (unsigned i = 0; i < Tdim; ++i)
      if (std::abs(acceleration_.col(phase)(i)) < tolerance)
        acceleration_.col(phase)(i) = 0.;
    status = true;
  }
  return status;
}

//! Compute acceleration and velocity with cundall damping factor
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::NodeXMPM<Tdim, Tdof, Tnphases>::compute_acceleration_velocity_cundall(
    unsigned phase, double dt, double damping_factor) noexcept {
  bool status = false;
  const double tolerance = 1.0E-15;
  if (mass_(phase) > tolerance) {
    // acceleration = (unbalaced force / mass)
    auto unbalanced_force =
        this->external_force_.col(phase) + this->internal_force_.col(phase);
    this->acceleration_.col(phase) =
        (unbalanced_force - damping_factor * unbalanced_force.norm() *
                                this->velocity_.col(phase).cwiseSign()) /
        this->mass_(phase);

    // Apply friction constraints
    this->apply_friction_constraints(dt);

    // Velocity += acceleration * dt
    this->velocity_.col(phase) += this->acceleration_.col(phase) * dt;
    // Apply velocity constraints, which also sets acceleration to 0,
    // when velocity is set.
    this->apply_velocity_constraints();
    // Set a threshold
    for (unsigned i = 0; i < Tdim; ++i)
      if (std::abs(velocity_.col(phase)(i)) < tolerance)
        velocity_.col(phase)(i) = 0.;
    for (unsigned i = 0; i < Tdim; ++i)
      if (std::abs(acceleration_.col(phase)(i)) < tolerance)
        acceleration_.col(phase)(i) = 0.;
    status = true;
  }
  return status;
}

//! Compute acceleration and velocity with cundall damping factor
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::NodeXMPM<Tdim, Tdof, Tnphases>::intergrate_momentum(
    unsigned phase, double dt) noexcept {
  

  momentum_.col(phase) = momentum_.col(phase) 
    + (internal_force_.col(phase)  + external_force_.col(phase)) * dt;
  if(enrich_h_)
    momentum_h_.col(phase) = momentum_h_.col(phase) 
      + (internal_force_h_.col(phase) + external_force_h_.col(phase)) * dt;
  // Apply velocity constraints, which also sets acceleration to 0,
  // when velocity is set.
  this->apply_velocity_constraints();

  this->self_contact_discontinuity(dt);

  this->apply_velocity_constraints();

 
  return true;
}

//! Assign velocity constraint
//! Constrain directions can take values between 0 and Dim * Nphases
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::NodeXMPM<Tdim, Tdof, Tnphases>::assign_velocity_constraint(
    unsigned dir, double velocity) {
  bool status = true;
  try {
    //! Constrain directions can take values between 0 and Dim * Nphases
    if (dir < (Tdim * Tnphases))
      this->velocity_constraints_.insert(std::make_pair<unsigned, double>(
          static_cast<unsigned>(dir), static_cast<double>(velocity)));
    else
      throw std::runtime_error("Constraint direction is out of bounds");

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Apply velocity constraints
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::NodeXMPM<Tdim, Tdof, Tnphases>::apply_velocity_constraints() {
  // Set velocity constraint
  for (const auto& constraint : this->velocity_constraints_) {
    // Direction value in the constraint (0, Dim * Nphases)
    const unsigned dir = constraint.first;
    // Direction: dir % Tdim (modulus)
    const auto direction = static_cast<unsigned>(dir % Tdim);
    // Phase: Integer value of division (dir / Tdim)
    const auto phase = static_cast<unsigned>(dir / Tdim);

    if (!generic_boundary_constraints_) {
      // Velocity constraints are applied on Cartesian boundaries
      this->velocity_(direction, phase) = constraint.second;
      this->momentum_(direction, phase) = this->mass(phase) * constraint.second;
      this->momentum_h_(direction, phase) = this->mass_h_(phase) * constraint.second;
      // Set acceleration to 0 in direction of velocity constraint
      this->acceleration_(direction, phase) = 0.;
    } else {
      // Velocity constraints on general boundaries
      // Compute inverse rotation matrix
      const Eigen::Matrix<double, Tdim, Tdim> inverse_rotation_matrix =
          rotation_matrix_.inverse();
      // Transform to local coordinate
      Eigen::Matrix<double, Tdim, Tnphases> local_velocity =
          inverse_rotation_matrix * this->velocity_;
      Eigen::Matrix<double, Tdim, Tnphases> local_acceleration =
          inverse_rotation_matrix * this->acceleration_;
      // Apply boundary condition in local coordinate
      local_velocity(direction, phase) = constraint.second;
      local_acceleration(direction, phase) = 0.;
      // Transform back to global coordinate
      this->velocity_ = rotation_matrix_ * local_velocity;
      this->acceleration_ = rotation_matrix_ * local_acceleration;
    }
  }
}

//! Apply velocity constraints
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::NodeXMPM<Tdim, Tdof, Tnphases>::self_contact_discontinuity(double dt) {

  if(!enrich_h_)
    return;
  
  unsigned  phase = 0; 
  const double tolerance = 1.0E-15;

  auto mass_positive = mass_.col(phase) + mass_h_.col(phase);
  auto mass_negative = mass_.col(phase) - mass_h_.col(phase);

  if(mass_positive(phase) < tolerance || mass_negative(phase) < tolerance)
    return;

  auto velocity_positive = (momentum_.col(phase) + momentum_h_.col(phase)) / mass_positive(phase);
  auto velocity_negative = (momentum_.col(phase) - momentum_h_.col(phase)) / mass_negative(phase);

  if((velocity_positive - velocity_negative).col(phase).dot(direction_discontinuity_.col(phase)) >= 0)
    return;

  auto momentum_contact = (mass_h_(phase)*momentum_.col(phase) - mass_(phase)*momentum_h_.col(phase)) / mass_(phase);
  auto force_contact = momentum_contact/dt;

  //! frictional_coef < 0: move together without slide
  if(frictional_coef < 0)
  {
    momentum_h_.col(phase) = momentum_h_.col(phase) + momentum_contact.col(phase);
    internal_force_h_.col(phase) = internal_force_h_.col(phase) + force_contact.col(phase);
  }
  else
  {
    double momentum_contact_norm = momentum_contact.col(phase).dot(direction_discontinuity_.col(phase));
    double force_contact_norm = momentum_contact_norm/dt;

    double max_frictional_force = frictional_coef * abs(force_contact_norm);

    auto momentum_tangential = momentum_contact.col(phase) - momentum_contact_norm*direction_discontinuity_.col(phase);
    auto force_tangential = momentum_tangential/dt;

    double force_tangential_value = force_tangential.norm();

    double frictional_force = force_tangential_value < max_frictional_force? force_tangential_value : max_frictional_force;

    //!adjust the momentum and force
    momentum_h_.col(phase) = momentum_h_.col(phase) + momentum_contact_norm*direction_discontinuity_.col(phase);
    internal_force_h_.col(phase) = internal_force_h_.col(phase) + force_contact_norm*direction_discontinuity_.col(phase);

    momentum_h_.col(phase) = momentum_h_.col(phase) + frictional_force*force_tangential.col(phase).normalized()*dt;
    internal_force_h_.col(phase) = internal_force_h_.col(phase) + frictional_force*force_tangential.col(phase).normalized();

  }
}

//! Assign friction constraint
//! Constrain directions can take values between 0 and Dim * Nphases
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::NodeXMPM<Tdim, Tdof, Tnphases>::assign_friction_constraint(
    unsigned dir, int sign_n, double friction) {
  bool status = true;
  try {
    //! Constrain directions can take values between 0 and Dim * Nphases
    if (dir < Tdim) {
      this->friction_constraint_ =
          std::make_tuple(static_cast<unsigned>(dir), static_cast<int>(sign_n),
                          static_cast<double>(friction));
      this->friction_ = true;
    } else
      throw std::runtime_error("Constraint direction is out of bounds");

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Apply friction constraints
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::NodeXMPM<Tdim, Tdof, Tnphases>::apply_friction_constraints(double dt) {
  if (friction_) {
    auto sign = [](double value) { return (value > 0.) ? 1. : -1.; };

    // Set friction constraint
    // Direction value in the constraint (0, Dim)
    const unsigned dir_n = std::get<0>(this->friction_constraint_);

    // Normal direction of friction
    const double sign_dir_n = sign(std::get<1>(this->friction_constraint_));

    // Friction co-efficient
    const double mu = std::get<2>(this->friction_constraint_);

    const unsigned phase = 0;

    // Acceleration and velocity
    double acc_n, acc_t, vel_t;

    if (Tdim == 2) {
      // tangential direction to boundary
      const unsigned dir_t = (Tdim - 1) - dir_n;

      if (!generic_boundary_constraints_) {
        // Cartesian case
        // Normal and tangential acceleration
        acc_n = this->acceleration_(dir_n, phase);
        acc_t = this->acceleration_(dir_t, phase);
        // Velocity tangential
        vel_t = this->velocity_(dir_t, phase);
      } else {
        // General case, transform to local coordinate
        // Compute inverse rotation matrix
        const Eigen::Matrix<double, Tdim, Tdim> inverse_rotation_matrix =
            rotation_matrix_.inverse();
        // Transform to local coordinate
        Eigen::Matrix<double, Tdim, Tnphases> local_acceleration =
            inverse_rotation_matrix * this->acceleration_;
        Eigen::Matrix<double, Tdim, Tnphases> local_velocity =
            inverse_rotation_matrix * this->velocity_;
        // Normal and tangential acceleration
        acc_n = local_acceleration(dir_n, phase);
        acc_t = local_acceleration(dir_t, phase);
        // Velocity tangential
        vel_t = local_velocity(dir_t, phase);
      }

      if ((acc_n * sign_dir_n) > 0.0) {
        if (vel_t != 0.0) {  // kinetic friction
          const double vel_net = dt * acc_t + vel_t;
          const double vel_frictional = dt * mu * std::abs(acc_n);
          if (std::abs(vel_net) <= vel_frictional)
            acc_t = -vel_t / dt;
          else
            acc_t -= sign(vel_net) * mu * std::abs(acc_n);
        } else {  // static friction
          if (std::abs(acc_t) <= mu * std::abs(acc_n))
            acc_t = 0.0;
          else
            acc_t -= sign(acc_t) * mu * std::abs(acc_n);
        }

        if (!generic_boundary_constraints_) {
          // Cartesian case
          this->acceleration_(dir_t, phase) = acc_t;
        } else {
          // Local acceleration in terms of tangential and normal
          Eigen::Matrix<double, Tdim, Tnphases> acc;
          acc(dir_t, phase) = acc_t;
          acc(dir_n, phase) = acc_n;

          // General case, transform to global coordinate
          this->acceleration_.col(phase) = rotation_matrix_ * acc.col(phase);
        }
      }
    } else if (Tdim == 3) {
      Eigen::Matrix<int, 3, 2> dir;
      dir(0, 0) = 1;
      dir(0, 1) = 2;  // tangential directions for dir_n = 0
      dir(1, 0) = 0;
      dir(1, 1) = 2;  // tangential directions for dir_n = 1
      dir(2, 0) = 0;
      dir(2, 1) = 1;  // tangential directions for dir_n = 2

      const unsigned dir_t0 = dir(dir_n, 0);
      const unsigned dir_t1 = dir(dir_n, 1);

      Eigen::Matrix<double, Tdim, 1> acc, vel;
      if (!generic_boundary_constraints_) {
        // Cartesian case
        acc = this->acceleration_.col(phase);
        vel = this->velocity_.col(phase);
      } else {
        // General case, transform to local coordinate
        // Compute inverse rotation matrix
        const Eigen::Matrix<double, Tdim, Tdim> inverse_rotation_matrix =
            rotation_matrix_.inverse();
        // Transform to local coordinate
        acc = inverse_rotation_matrix * this->acceleration_.col(phase);
        vel = inverse_rotation_matrix * this->velocity_.col(phase);
      }

      const auto acc_n = acc(dir_n);
      auto acc_t =
          std::sqrt(acc(dir_t0) * acc(dir_t0) + acc(dir_t1) * acc(dir_t1));
      const auto vel_t =
          std::sqrt(vel(dir_t0) * vel(dir_t0) + vel(dir_t1) * vel(dir_t1));

      if (acc_n * sign_dir_n > 0.0) {
        // kinetic friction
        if (vel_t != 0.0) {
          Eigen::Matrix<double, 2, 1> vel_net;
          // friction is applied opposite to the vel_net
          vel_net(0) = vel(dir_t0) + acc(dir_t0) * dt;
          vel_net(1) = vel(dir_t1) + acc(dir_t1) * dt;
          const double vel_net_t =
              sqrt(vel_net(0) * vel_net(0) + vel_net(1) * vel_net(1));
          const double vel_fricion = mu * std::abs(acc_n) * dt;

          if (vel_net_t <= vel_fricion) {
            acc(dir_t0) = -vel(dir_t0);  // To set particle velocity to zero
            acc(dir_t1) = -vel(dir_t1);
          } else {

            acc(dir_t0) -= mu * std::abs(acc_n) * (vel_net(0) / vel_net_t);
            acc(dir_t1) -= mu * std::abs(acc_n) * (vel_net(1) / vel_net_t);
          }
        } else {                                // static friction
          if (acc_t <= mu * std::abs(acc_n)) {  // since acc_t is positive
            acc(dir_t0) = 0;
            acc(dir_t1) = 0;
          } else {
            acc_t -= mu * std::abs(acc_n);
            acc(dir_t0) -= mu * std::abs(acc_n) * (acc(dir_t0) / acc_t);
            acc(dir_t1) -= mu * std::abs(acc_n) * (acc(dir_t1) / acc_t);
          }
        }

        if (!generic_boundary_constraints_) {
          // Cartesian case
          this->acceleration_.col(phase) = acc;
        } else {
          // General case, transform to global coordinate
          this->acceleration_.col(phase) = rotation_matrix_ * acc;
        }
      }
    }
  }
}

//! Add material id from material points to material_ids_
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::NodeXMPM<Tdim, Tdof, Tnphases>::append_material_id(unsigned id) {
  std::lock_guard<std::mutex> guard(node_mutex_);
  material_ids_.emplace(id);
}

// Assign MPI rank to node
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::NodeXMPM<Tdim, Tdof, Tnphases>::mpi_rank(unsigned rank) {
  std::lock_guard<std::mutex> guard(node_mutex_);
  auto status = this->mpi_ranks_.insert(rank);
  return status.second;
}

//! Update nodal property at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::NodeXMPM<Tdim, Tdof, Tnphases>::update_property(
    bool update, const std::string& property,
    const Eigen::MatrixXd& property_value, unsigned mat_id,
    unsigned nprops) noexcept {
  // Update/assign property
  std::lock_guard<std::mutex> guard(node_mutex_);
  property_handle_->update_property(property, prop_id_, mat_id, property_value,
                                    nprops);
}