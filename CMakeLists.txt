cmake_minimum_required(VERSION 3.0)
#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

project(apothesis)

set(header_files
    ./src/apothesis.h
    ./src/pointers.h
    ./src/IO/io.h
    ./src/processes/abstract_process.h
    ./src/lattice/lattice.h
    ./src/processes/process.h
    ./src/lattice/site.h
    ./src/lattice/FCC.h
    ./src/lattice/HCP.h
    ./src/lattice/SimpleCubic.h
    ./src/processes/adsorption.h
    ./src/processes/adsorption.h
    ./src/processes/diffusion.h
    ./src/processes/diffusion.h
    ./src/processes/factory_process.h
    ./src/processes/desorption.h
    ./src/processes/desorption.h
    ./src/processes/abstract_process.h
    ./src/processes/reaction.h
    ./src/error/errorhandler.h
    ./src/processes/parameters.h
    ./src/IO/xyz_reader.h
    ./src/IO/cml_reader.h
    ./src/IO/reader.h
    ./src/IO/io.h
    ./src/properties.h
    ./src/extLibs/random_generator.h
    ./src/extLibs/randomc.h
)
set(essential_src_files
    ./src/main.cpp
    ./src/properties.cpp
    ./src/apothesis.cpp
)
set(IO_files
    ./src/IO/xyz_reader.cpp
    ./src/IO/cml_reader.cpp
    ./src/IO/reader.cpp
    ./src/IO/io.cpp
 )
set(extLibs_files
    ./src/extLibs/random_generator.cpp
    ./src/extLibs/mersenne.cpp
)
set(process_files
    ./src/processes/adsorption.cpp
    ./src/processes/diffusion.cpp
    ./src/processes/adsorption.cpp
    ./src/processes/diffusion.cpp
    ./src/processes/factory_process.cpp
    ./src/processes/process.cpp
    ./src/processes/parameters.cpp
    ./src/IO/io.cpp
    ./src/IO/xyz_reader.cpp
    ./src/IO/cml_reader.cpp
    ./src/IO/reader.cpp
    ./src/processes/desorption.cpp
    ./src/processes/process.cpp
    ./src/processes/abstract_process.cpp
    ./src/processes/reaction.h
    ./src/processes/reaction.cpp
)
set(error_files
    ./src/error/errorhandler.cpp 
)

set(lattice_files
    ./src/lattice/site.cpp
    ./src/lattice/lattice.cpp
    ./src/lattice/FCC.cpp
    ./src/lattice/HCP.cpp
    ./src/lattice/SimpleCubic.cpp
)

add_executable(${PROJECT_NAME} ./src/main.cpp
    ${header_files}
    ${process_files}
    ${error_files}
    ${IO_files}
    ${lattice_files}
    ${species_files}
    ${extLibs_files}
    ${essential_src_files}
)

target_include_directories(${PROJECT_NAME} PUBLIC
    .
    ./src/
    ./src/error
    ./src/processes
    ./src/IO
    ./src/lattice
    ./src/species
)
