cmake_minimum_required(VERSION 3.13.2)

project(1dim_sic_simulator)
enable_language(Fortran)

set(EXECUTABLE main.exe)

add_executable(
    ${EXECUTABLE}
    main.f90
    parameters.f90
    commons.f90
    efield.f90
    particles.f90
    boundary.f90
    progressbar.f90
    status.f90
    utils.f90
    argparse.f90
    fft.f90
    sic.f90
)

find_package(OpenMP REQUIRED)

target_link_libraries(
    ${EXECUTABLE}
    PUBLIC
    OpenMP::OpenMP_Fortran
)