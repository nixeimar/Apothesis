cmake_minimum_required(VERSION 3.16)
project(apothesis VERSION 1.0 LANGUAGES C CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(QT NAMES Qt5 Qt6 REQUIRED COMPONENTS Core)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Gui)

qt_standard_project_setup()

qt_add_executable(apothesis WIN32 MACOSX_BUNDLE
    src/IO/cml_reader.cpp src/IO/cml_reader.h
    src/IO/io.cpp src/IO/io.h
    src/IO/xyz_reader.cpp src/IO/xyz_reader.h
    src/apothesis.cpp src/apothesis.h
    src/error/errorhandler.cpp src/error/errorhandler.h
    src/extLibs/mersenne.cpp
    src/extLibs/random_generator.cpp src/extLibs/random_generator.h
    src/extLibs/randomc.h
    src/lattice/FCC.cpp src/lattice/FCC.h
    src/lattice/HCP.cpp src/lattice/HCP.h
    src/lattice/SimpleCubic.cpp src/lattice/SimpleCubic.h
    src/lattice/diamond.cpp src/lattice/diamond.h
    src/lattice/lattice.cpp src/lattice/lattice.h
    src/lattice/site.cpp src/lattice/site.h
    src/macroprocesses/ald.cpp src/macroprocesses/ald.h
    src/macroprocesses/ale.cpp src/macroprocesses/ale.h
    src/macroprocesses/catalysis.cpp src/macroprocesses/catalysis.h
    src/macroprocesses/cvd.cpp src/macroprocesses/cvd.h
    src/macroprocesses/etching.cpp src/macroprocesses/etching.h
    src/macroprocesses/macroprocess.cpp src/macroprocesses/macroprocess.h
    src/macroprocesses/pvd.cpp src/macroprocesses/pvd.h
    src/main.cpp
    src/microprocesses/abstract_process.cpp src/microprocesses/abstract_process.h
    src/microprocesses/adsorption.cpp src/microprocesses/adsorption.h
    src/microprocesses/adsorption_perform.cpp src/microprocesses/adsorption_perform.h
    src/microprocesses/adsorption_rules.cpp src/microprocesses/adsorption_rules.h
    src/microprocesses/adsorption_types.cpp src/microprocesses/adsorption_types.h
    src/microprocesses/desorption.cpp src/microprocesses/desorption.h
    src/microprocesses/desorption_perform.cpp src/microprocesses/desorption_perform.h
    src/microprocesses/desorption_rules.cpp src/microprocesses/desorption_rules.h
    src/microprocesses/desorption_types.cpp src/microprocesses/desorption_types.h
    src/microprocesses/diffusion.cpp src/microprocesses/diffusion.h
    src/microprocesses/diffusion_perform.cpp src/microprocesses/diffusion_perform.h
    src/microprocesses/diffusion_rules.cpp src/microprocesses/diffusion_rules.h
    src/microprocesses/diffusion_types.cpp src/microprocesses/diffusion_types.h
    src/microprocesses/factory_process.cpp src/microprocesses/factory_process.h
    src/microprocesses/io.h
    src/microprocesses/parameters.cpp src/microprocesses/parameters.h
    src/microprocesses/process.cpp src/microprocesses/process.h
    src/microprocesses/reaction.cpp src/microprocesses/reaction.h
    src/pointers.h
    src/properties.cpp src/properties.h
    src/register.h
)
target_include_directories(apothesis PRIVATE
    .
    src
    src/IO
    src/error
    src/lattice
    src/macroprocesses
    src/microprocesses
)

target_link_libraries(apothesis PRIVATE
    Qt::Core
    Qt::Gui
)

target_compile_options(apothesis
    -std=c++17
)

install(TARGETS apothesis
    BUNDLE DESTINATION .
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

qt_generate_deploy_app_script(
    TARGET apothesis
    FILENAME_VARIABLE deploy_script
    NO_UNSUPPORTED_PLATFORM_ERROR
)
install(SCRIPT ${deploy_script})
