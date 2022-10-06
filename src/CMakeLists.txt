cmake_minimum_required(VERSION 3.0)
#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

find_package(RapidJSON)

project(Apothesis)

set(header_files
    apothesis.h
    pointers.h
    IO/io.h
    processes/abstract_process.h
    lattice/lattice.h
    processes/process.h
    lattice/site.h
    lattice/FCC.h
    lattice/BCC.h
    processes/adsorption.h
    processes/adsorption_simple_cubic.h
    processes/diffusion_simple_cubic.h
    processes/diffusion.h
    processes/factory_process.h
    processes/desorption.h
    processes/desorption_simple_cubic.h
    processes/abstract_process.h
    processes/reaction_new.h
    error/errorhandler.h
    processes/parameters.h
    IO/read.h
    species/species.h
    species/species_new.h
    properties.h
)
set(essential_src_files
    properties.cpp
    apothesis.cpp
)
set(IO_files
    IO/read.cpp
)
set(aux_files
    aux/random_generator.h
    aux/randomc.h
    aux/random_generator.cpp
    aux/mersenne.cpp
)
set(process_files
    processes/adsorption.cpp
    processes/diffusion.cpp
    processes/adsorption_simple_cubic.cpp
    processes/diffusion_simple_cubic.cpp
    processes/factory_process.cpp
    processes/process.cpp
    processes/parameters.cpp
    IO/io.cpp
    processes/desorption.cpp
    processes/desorption_simple_cubic.cpp
    processes/process.cpp
    processes/abstract_process.cpp
    processes/reaction_new.h
    processes/reaction_new.cpp
)
set(error_files
    error/errorhandler.cpp 
)

set(lattice_files
    lattice/site.cpp
    lattice/lattice.cpp
    lattice/FCC.cpp
    lattice/BCC.cpp
)

set(species_files
    species/species.cpp
    species/species_new.cpp
)
add_executable(${PROJECT_NAME} "main.cpp"
    ${header_files}
    ${process_files}
    ${error_files}
    ${IO_files}
    ${lattice_files}
    ${species_files}
    ${aux_files}
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
