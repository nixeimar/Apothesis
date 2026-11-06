//============================================================================
//    Apothesis: A kinetic Monte Calro (KMC) code for deposition processes.
//    Copyright (C) 2019  Nikolaos (Nikos) Cheimarios
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.

//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.
//============================================================================

/*! \mainpage My Personal Index Page
 *
 * \section intro_sec Introduction
 * Apothesis is a generalized kinteic Monte Carlo Code for tackling realistic deposition porcesses
 * and more specificaly Chemical Vapor Deposition (CVD) and Atomic Layer Deposition (ALD) processes.
 * It is developed in C++ and it is disrtibuted under GNU license.
 *
 *
 * \section install_sec Installation
 * Currently no special installations instructions are required.
 * A simple "make" should create the executable "Apothesis".
 * The Makefile was generated with qmake but the qt libradies are excluded.
 *
 */
#include <iostream>
#include <list>
#include <memory>
#include <filesystem>
#include "config/configuration.h"
#include "site.h"
#include "lattice.h"
#include "process.h"
#include "apothesis.h"

using namespace std;
using namespace MicroProcesses;

int main(int argc, char* argv[])
{
    try {
        if (argc == 2 && string(argv[1]) == "--help") {
            cout << "Usage: apothesis [--input|-i] <input.json>\n"
                 << "       apothesis --validate-input <input.json>\n"
                 << "       apothesis --build-surface <input.json> [--output <surface.xyz>]\n"
                 << "With no arguments, reads input.json in the working directory.\n";
            return 0;
        }
        if ((argc == 3 || argc == 5) && string(argv[1]) == "--build-surface") {
            if (argc == 5 && string(argv[3]) != "--output")
                throw std::runtime_error("Usage: apothesis --build-surface <input.json> [--output <surface.xyz>]");
            const auto config = Config::readJson(argv[2]);
            if (!config.atomisticSurface && !config.diamondCVD && !(config.surfaceKMC && config.surfaceKMC->customGeometry))
                throw std::runtime_error("--build-surface requires lattice.representation: atomistic");
            const auto surface = config.diamondCVD ? config.diamondCVD->lattice : config.surfaceKMC && config.surfaceKMC->customGeometry ? *config.surfaceKMC->customGeometry : Atomistic::build(*config.atomisticSurface);
            const string output = argc == 5 ? argv[4] : "surface.xyz";
            surface.writeXYZ(output);
            cout << "Created " << surface.atomCount() << " atoms on " << surface.sites.size()
                 << " crystal sites in " << output << " (angstrom).\n";
            return 0;
        }
        if (argc == 3 && string(argv[1]) == "--validate-input") {
            Config::readJson(argv[2]);
            cout << "Input is valid.\n";
            return 0;
        }
        string input="input.json";
        if (argc==2 && argv[1][0]!='-') input=argv[1];
        else if (argc==3 && (string(argv[1])=="--input" || string(argv[1])=="-i")) input=argv[2];
        else if (argc!=1) throw std::runtime_error("Usage: apothesis [--input|-i] <input.json>");
        const auto configuration=Config::readJson(input);
        if (configuration.diamondCVD) {
            cout<<"Running diamond (111) covalent C/H growth ..."<<endl;
            DiamondCVD::runToFiles(*configuration.diamondCVD);
            cout<<"Finished. See diamond_summary.json, diamond.log and diamond_*.xyz in "<<std::filesystem::current_path()<<endl;
            return 0;
        }
        if (configuration.surfaceKMC) {
            cout<<"Running explicit surface processes on 3D site occupancy ..."<<endl;
            SurfaceKMC::runToFiles(*configuration.surfaceKMC);
            cout<<"Finished. Output directory: "<<std::filesystem::current_path()<<endl;
            cout<<"See surface_summary.json and surface.log";
            if (configuration.surfaceKMC->writeXYZ) cout<<", surface_initial.xyz, surface_final.xyz and surface_trajectory.xyz";
            cout<<'.'<<endl;return 0;
        }
        if (configuration.zgb3d) {
            cout<<"Running ZGB chemistry on explicit 3D site occupancy ..."<<endl;
            ZGB3D::runToFiles(*configuration.zgb3d);
            cout<<"Finished. See zgb3d_summary.json and zgb3d.log."<<endl;
            return 0;
        }
        if (configuration.zgb) {
            cout<<"Running Cheimarios 2023 height-field ZGB model ..."<<endl;
            ZGB::runToFiles(*configuration.zgb);
            cout<<"Finished. See zgb_summary.json and zgb.log."<<endl;
            return 0;
        }
        if (configuration.atomisticSurface) {
            if (!configuration.atomisticSimulation)
                throw std::runtime_error("Atomistic geometry-only input: add system, simulation, processes and output to run KMC, or use --build-surface.");
            cout<<"Running atomistic KMC ..."<<endl;
            Atomistic::runToFiles(Atomistic::build(*configuration.atomisticSurface),*configuration.atomisticSimulation);
            cout<<"Atomistic KMC finished. See atomistic.log and atomistic_trajectory.xyz."<<endl;
            return 0;
        }
        unique_ptr<Apothesis> apothesis(new Apothesis(argc, argv));
        cout << "Initiating Apothesis" << endl;
        apothesis->init();
        cout << "Apothesis running ..." << endl;
        apothesis->exec();
        cout << "Apothesis finished successfully." << endl;
        return 0;
    } catch (const std::exception& error) {
        cerr << "Input/simulation error: " << error.what() << endl;
        return 1;
    }
}
