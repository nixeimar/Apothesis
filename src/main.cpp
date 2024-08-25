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
#include "site.h"
#include "lattice.h"
#include "process.h"
#include "apothesis.h"
#include "src/IO/io.h"

using namespace std;
using namespace MicroProcesses;

int main( int argc, char* argv[] )
{

    //Checking commit
    Apothesis* apothesis = new Apothesis( argc, argv );

    cout << "Initiating Apothesis" << endl;
    apothesis->init();

    cout << "Apothesis runnning ..." << endl;
    apothesis->exec();
    cout << "Apothesis finished succesfully." << endl;

    int run = apothesis->pParameters->getCurrRunNum();
    int total_runs = apothesis->pParameters->getRuns();
    while(run!=total_runs){
      //reinitialize heights and species
      //current height file
      // string latestHeightFileLocation = apothesis->pParameters->getLatestHeightFileLocation(); 
      // vector<vector<int>>currHeights = apothesis->pIO->readHeightFile(latestHeightFileLocation);
      // apothesis->pLattice->setInitialHeightAllSites(latestHeightFileLocation);

      // string latestSpeciesFileLocation = apothesis->pParameters->getLatestSpeciesFileLocation();
      // vector<vector<string>>

      // //apothesis->init()->initnewSurface();
      apothesis->exec();
      apothesis->pParameters->setCurrRunNum(run++);
    }

    if ( apothesis )
      delete apothesis;
}


