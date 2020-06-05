cmake_minimum_required(VERSION 3.0)
#set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

project(Apothesis)

# Recurse into the "test" subdirectory. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.


set(header_files
    apothesis.h
    pointers.h
    register.h
    IO/io.h
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
)
set(essential_src_files
    apothesis.cpp
    register.cpp
)
set(IO_files
    IO/io.cpp
)
set(process_files
    processes/adsorption.cpp
    processes/diffusion.cpp
    processes/factory_process.cpp
    processes/process.cpp
    processes/SurfaceReaction.cpp
)
set(error_files
    error/errorhandler.cpp 
)

set(lattice_files
    lattice/site.cpp
    lattice/lattice.cpp
)


#add_library(Apo 
#    ${header_files}
#    ${error_files}
#    ${process_files}
#    ${IO_files}
#    ${lattice_files}
#    ${essential_src_files}
#)
#
#target_include_directories(Apo PRIVATE
#    .
#    error
#    lattice
#    processes
#    IO
#)

add_executable(${PROJECT_NAME} "main.cpp"
    ${header_files}
    ${process_files}
    ${error_files}
    ${IO_files}
    ${lattice_files}
    ${essential_src_files}
)

target_include_directories(${PROJECT_NAME} PUBLIC
    .
    error
    processes
    IO
    lattice
)

#target_link_libraries(${PROJECT_NAME} Apo)