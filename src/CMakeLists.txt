cmake_minimum_required(VERSION 3.0)
#set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

find_package(RapidJSON)

project(Apothesis)

set(header_files
    apothesis.h
    pointers.h
    register.h
    processes/io.h
    processes/abstract_process.h
    lattice/lattice.h
    processes/process.h
    lattice/site.h
    processes/adsorption.h 
    processes/diffusion.h
    processes/factory_process.h
    processes/SurfaceReaction.h
    error/errorhandler.h
    processes/parameters.h
    IO/read.h
    species/species.h
)
set(essential_src_files
    apothesis.cpp
    register.cpp
)
set(IO_files
    IO/read.cpp
)
set(process_files
    processes/adsorption.cpp
    processes/diffusion.cpp
    processes/factory_process.cpp
    processes/process.cpp
    processes/SurfaceReaction.cpp
    processes/parameters.cpp
    processes/io.cpp
)
set(error_files
    error/errorhandler.cpp 
)

set(lattice_files
    lattice/site.cpp
    lattice/lattice.cpp
)

set(species_files
    species/species.cpp
)
add_executable(${PROJECT_NAME} "main.cpp"
    ${header_files}
    ${process_files}
    ${error_files}
    ${IO_files}
    ${lattice_files}
    ${species_files}
    ${essential_src_files}
)

target_include_directories(${PROJECT_NAME} PUBLIC
    .
    error
    processes
    IO
    lattice
    species
)
