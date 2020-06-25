#include "discontinuity_base.h"
#include "factory.h"
#include "discontinuity_3d.h"

// Triangle 3-noded element
static Register<mpm::DiscontinuityBase<3>, mpm::Discontinuity_3D<3>> tri3d("tri3d");
