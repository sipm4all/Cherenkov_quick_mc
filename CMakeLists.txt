cmake_minimum_required(VERSION 3.19)

project(cherenkov_mc
    VERSION     1.0.0
    DESCRIPTION "Cherenkov fast Monte Carlo — ROOT-backed, built on mist / mist-hep"
    LANGUAGES   CXX
)

# ------------------------------------------------------------------
# Dependencies
#
# cherenkov_mc consumes:
#   - mist        (logger / utilities)        — via mist::mist
#   - mist-hep    (ROOT-backed helpers, RNG)  — via mist::hep
#
# mist::hep links mist::mist and ROOT publicly, so pulling in mist-hep is
# enough to get the logger, the process-wide RNG, and ROOT itself.
#
# Two consumption modes:
#   1. Standalone: find_package(mist-hep) (which re-finds mist + ROOT).
#   2. As a sibling subdir under a parent build that already added mist-hep
#      (mist::hep target already exists — skip the find_package).
# ------------------------------------------------------------------
if(NOT TARGET mist::hep AND NOT TARGET mist::mist-hep)
    find_package(mist-hep 1 REQUIRED)   # SameMajorVersion: any mist-hep 1.x; transitively brings mist + ROOT
endif()

# Normalise the target name: find_package() exports mist::mist-hep, while
# add_subdirectory() exposes the mist::hep alias. Link whichever exists.
if(TARGET mist::mist-hep)
    set(CHERENKOV_MC_MISTHEP mist::mist-hep)
else()
    set(CHERENKOV_MC_MISTHEP mist::hep)
endif()

# ------------------------------------------------------------------
# Library target
# ------------------------------------------------------------------
add_library(cherenkov_mc STATIC
    src/cherenkov_mc/physics_vector.cxx
    src/cherenkov_mc/particle.cxx
    src/cherenkov_mc/cherenkov.cxx
    src/cherenkov_mc/probability.cxx
    src/cherenkov_mc/photon_interactions.cxx
    src/cherenkov_mc/geometry.cxx
    src/cherenkov_mc/fast_mc.cxx
)

target_include_directories(cherenkov_mc
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_compile_features(cherenkov_mc PUBLIC cxx_std_20)

target_link_libraries(cherenkov_mc
    PUBLIC
        ${CHERENKOV_MC_MISTHEP}   # mist::mist-hep (find_package) or mist::hep (subdir) -> logger + ROOT
)

# Namespaced alias — consumers write `cherenkov_mc::cherenkov_mc` regardless
# of whether they pulled us in via find_package() or add_subdirectory().
add_library(cherenkov_mc::cherenkov_mc ALIAS cherenkov_mc)

# ------------------------------------------------------------------
# Install rules — mirrors the mist / mist-hep shape
# ------------------------------------------------------------------
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(TARGETS cherenkov_mc
    EXPORT      cherenkov_mcTargets
    ARCHIVE     DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY     DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME     DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT cherenkov_mcTargets
    FILE        cherenkov_mcTargets.cmake
    NAMESPACE   cherenkov_mc::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cherenkov_mc
)

export(EXPORT cherenkov_mcTargets
    NAMESPACE cherenkov_mc::
    FILE      "${CMAKE_CURRENT_BINARY_DIR}/cherenkov_mcTargets.cmake"
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/cherenkov_mcConfigVersion.cmake"
    VERSION       ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cherenkov_mcConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cherenkov_mcConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cherenkov_mc
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/cherenkov_mcConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/cherenkov_mcConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cherenkov_mc
)

# ------------------------------------------------------------------
# Example(s)
# ------------------------------------------------------------------
option(CHERENKOV_MC_BUILD_EXAMPLES "Build cherenkov_mc example executables" OFF)

if(CHERENKOV_MC_BUILD_EXAMPLES)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
    function(cherenkov_mc_add_example name)
        add_executable(${name} example/${name}.cxx)
        target_link_libraries(${name} PRIVATE cherenkov_mc::cherenkov_mc)
        target_compile_features(${name} PRIVATE cxx_std_20)
    endfunction()

    cherenkov_mc_add_example(run)
    cherenkov_mc_add_example(fit_optical_curves)
    message(STATUS "[cherenkov_mc] Examples enabled — binaries in ${CMAKE_BINARY_DIR}/bin")
endif()

# ------------------------------------------------------------------
# Tests
# ------------------------------------------------------------------
option(CHERENKOV_MC_BUILD_TESTS "Build cherenkov_mc test executables" OFF)
enable_testing()

if(CHERENKOV_MC_BUILD_TESTS)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

    function(cherenkov_mc_add_test name)
        add_executable(test_${name} test/tester_${name}.cxx)
        target_link_libraries(test_${name} PRIVATE cherenkov_mc::cherenkov_mc)
        target_compile_features(test_${name} PRIVATE cxx_std_20)
        add_test(NAME test_${name} COMMAND test_${name})
    endfunction()

    cherenkov_mc_add_test(physics_vector)
    cherenkov_mc_add_test(geometry)
    cherenkov_mc_add_test(propagation)

    message(STATUS "[cherenkov_mc] Tests enabled — binaries in ${CMAKE_BINARY_DIR}/bin")
endif()
