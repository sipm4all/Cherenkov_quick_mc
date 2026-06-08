// SPDX-License-Identifier: MIT
//
// cherenkov_mc/cherenkov_mc.h — top-level umbrella header.
//
// Include this to pull in the whole library; prefer the targeted includes
// below in real code.
//
//   #include <cherenkov_mc/cherenkov_mc.h>        // everything
//   #include <cherenkov_mc/physics_vector.h>      // four-vector only
//   #include <cherenkov_mc/fast_mc.h>             // radiators + propagator
//
// Cascade:
//   cherenkov_mc.h
//   ├── constants.h
//   ├── physics_vector.h
//   ├── particle.h
//   ├── dynamics.h
//   ├── cherenkov.h            (cherenkov_mc::cherenkov)
//   ├── probability.h          (cherenkov_mc::probability)
//   ├── photon_interactions.h  (cherenkov_mc::photon_interactions)
//   ├── geometry.h             (box, cylinder)
//   └── fast_mc.h              (ch_radiator, fast_mc_event, propagator)
//
#pragma once

#include "cherenkov_mc/constants.h"
#include "cherenkov_mc/physics_vector.h"
#include "cherenkov_mc/particle.h"
#include "cherenkov_mc/dynamics.h"
#include "cherenkov_mc/cherenkov.h"
#include "cherenkov_mc/probability.h"
#include "cherenkov_mc/photon_interactions.h"
#include "cherenkov_mc/geometry.h"
#include "cherenkov_mc/fast_mc.h"
