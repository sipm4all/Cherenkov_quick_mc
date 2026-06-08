# cherenkov_mc

A small **Cherenkov fast Monte Carlo**: it shoots charged particles through
optical radiators (e.g. aerogel), emits Cherenkov photons via a simplified
Frank–Tamm yield, and tracks them through absorption / scattering until they
are detected or lost.

It is a ROOT-backed C++20 static library built on top of
[`mist`](../mist) (logger / utilities) and [`mist-hep`](../mist-hep)
(ROOT helpers, process-wide RNG). All textual output goes through
`mist::logger`; all randomness is drawn from `mist::hep::rng()`.

## Layout

```
include/cherenkov_mc/      public headers  (namespace cherenkov_mc)
  cherenkov_mc.h           umbrella header
  constants.h              physical constants / particle masses
  physics_vector.h         (t,x,y,z) four-vector with Minkowski metric
  particle.h               relativistic particle + straight-line stepping
  dynamics.h               beta / gamma helpers
  cherenkov.h              emission angle, thresholds, Frank–Tamm yield
  probability.h            Poisson samplers (Knuth / Atkinson)
  photon_interactions.h    absorption / scattering / reflection models
  geometry.h               box / cylinder volumes
  fast_mc.h                radiators, event, propagator
src/cherenkov_mc/          implementation translation units (.cxx)
example/                   run, fit_optical_curves (compiled)
test/                      ctest testers
```

The directory under `include/cherenkov_mc/` maps to the `cherenkov_mc`
namespace; helper namespaces are `cherenkov_mc::cherenkov`,
`::dynamics`, `::probability`, `::photon_interactions`.

## Dependencies

| Dependency | Provides                              | How found            |
|------------|---------------------------------------|----------------------|
| ROOT ≥ 6.24 | `TF1`, `TGraph`, histograms, fitting  | via `mist-hep`       |
| `mist`     | `mist::logger`, utilities             | `find_package(mist)` |
| `mist-hep` | `mist::hep::rng()`, ROOT linkage      | `find_package(mist-hep)` |

`mist-hep` re-resolves `mist` and ROOT itself, so consuming `mist-hep` is
enough to pull in the whole stack.

**C++ standard.** The library is built as C++20. If your ROOT was configured
for C++17 (the Homebrew default), every translation unit that includes a ROOT
header emits a one-line standard-mismatch warning. It is benign — the code
links and the test suite passes — but to silence it, build ROOT against C++20
or set `CMAKE_CXX_STANDARD=17` for this project.

## Build

```sh
cmake -B build -S . \
  -DCMAKE_PREFIX_PATH="$HOME/.local;$HOME/Analysis/mist-hep/build;$HOME/Analysis/mist/build" \
  -DCHERENKOV_MC_BUILD_TESTS=ON \
  -DCHERENKOV_MC_BUILD_EXAMPLES=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

(`mist` is consumed from its install prefix `~/.local`; `mist-hep` from its
build tree. Adjust `CMAKE_PREFIX_PATH` to wherever they live.)

## Run the examples

```sh
# Full simulation: 11.5 GeV e- through two aerogel radiators -> spectra
./build/bin/run 100 data/PDE.root run.root

# Calibrate the optical models against a measured "clarity" dataset
./build/bin/fit_optical_curves data/clarity_jap2022_AG22J003.root
```

## Use as a library

```cmake
find_package(cherenkov_mc REQUIRED)
target_link_libraries(my_target PRIVATE cherenkov_mc::cherenkov_mc)
```

```cpp
#include <cherenkov_mc/cherenkov_mc.h>

cherenkov_mc::fast_mc_event event;
event.set_geometry(geometry);
cherenkov_mc::propagator(event).run_event();
```

## Roadmap

Deferred past v1.0.0 (tracked as `TODO` markers in the headers):

- **Geometry rotation / translation** (`include/cherenkov_mc/geometry.h`) —
  the box and cylinder volumes are axis-aligned; full rigid-body placement
  with matching `is_inside` / step-size handling is not yet wired up. The
  position and axes infrastructure is present but unused.
- **PDG-driven particle kinematics** (`include/cherenkov_mc/fast_mc.h`) —
  `ch_radiator` does not yet derive particle mass / charge from a supplied
  PDG id, and there is no per-event / per-particle time model. The current
  propagator works for the single hard-coded particle hypothesis.
