cmake_minimum_required(VERSION 3.10)
#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if (WIN32 AND NOT CMAKE_GENERATOR MATCHES "MinGW")
    message(WARNING "On Windows, consider using MinGW or MSYS2 for building.")
endif()

project(apothesis)

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()


set(header_files
    ./src/apothesis.h
    ./src/pointers.h
    ./src/IO/io.h
    ./src/processes/abstract_process.h
    ./src/lattice/lattice.h
    ./src/processes/process.h
    ./src/lattice/site.h
    ./src/lattice/diamond.h
    ./src/lattice/FCC.h
    ./src/lattice/HCP.h
    ./src/lattice/SimpleCubic.h
    ./src/lattice/diamond.h
    ./src/processes/adsorption.h
    ./src/processes/diffusion.h
    ./src/processes/factory_process.h
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
    ./src/processes/adsorption_perform.h 
    ./src/processes/adsorption_rules.h 
    ./src/processes/adsorption_types.h
    ./src/processes/adsorption_perform.h 
    ./src/processes/adsorption_rules.h 
    ./src/processes/adsorption_types.h 
    ./src/processes/desorption_perform.h 
    ./src/processes/desorption_rules.h 
    ./src/processes/desorption_types.h 
    ./src/processes/diffusion_perform.h 
    ./src/processes/diffusion_rules.h 
    ./src/processes/diffusion_types.h

)
set(essential_src_files
    ./src/main.cpp
    ./src/properties.cpp
    ./src/apothesis.cpp
)
set(IO_files
    ./src/config/configuration.cpp
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
    ./src/processes/atomistic_kmc.cpp
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
    ./src/processes/adsorption_perform.cpp 
    ./src/processes/adsorption_rules.cpp 
    ./src/processes/adsorption_types.cpp 
    ./src/processes/desorption_perform.cpp 
    ./src/processes/desorption_rules.cpp 
    ./src/processes/desorption_types.cpp 
    ./src/processes/diffusion_perform.cpp 
    ./src/processes/diffusion_rules.cpp 
    ./src/processes/diffusion_types.cpp
)
set(error_files
    ./src/error/errorhandler.cpp 
)

set(lattice_files
    ./src/lattice/atomistic_surface.cpp
    ./src/lattice/site.cpp
    ./src/lattice/lattice.cpp
    ./src/lattice/diamond.cpp
    ./src/lattice/FCC.cpp
    ./src/lattice/HCP.cpp
    ./src/lattice/SimpleCubic.cpp
    ./src/lattice/diamond.cpp
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

include(CTest)
if(BUILD_TESTING)
    add_executable(input_tests tests/input_tests.cpp src/config/configuration.cpp src/lattice/atomistic_surface.cpp src/processes/atomistic_kmc.cpp)
    add_executable(atomistic_tests tests/atomistic_tests.cpp src/config/configuration.cpp src/lattice/atomistic_surface.cpp src/processes/atomistic_kmc.cpp)
    add_executable(atomistic_kmc_tests tests/atomistic_kmc_tests.cpp src/config/configuration.cpp src/lattice/atomistic_surface.cpp src/processes/atomistic_kmc.cpp)
    add_test(NAME atomistic_kmc COMMAND atomistic_kmc_tests ${CMAKE_CURRENT_SOURCE_DIR}/examples/atomistic-kmc.json)
    add_test(NAME atomistic_surfaces COMMAND atomistic_tests ${CMAKE_CURRENT_SOURCE_DIR}/examples/atomistic-cu.json)
    add_test(NAME json_input COMMAND input_tests ${CMAKE_CURRENT_SOURCE_DIR}/examples/height-field.json)
    find_package(PythonInterp 3 QUIET)
    if(PYTHONINTERP_FOUND)
        add_test(NAME input_integration
            COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/input_integration.py
                    $<TARGET_FILE:apothesis> ${CMAKE_CURRENT_SOURCE_DIR})
    endif()
endif()
