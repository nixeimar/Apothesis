//============================================================================
//    Apothesis: A kinetic Monte Calro (KMC) code for deposotion processes.
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
//#include "site.h"
//#include "lattice.h"
//#include "process.h"
#include "apothesis.h"

/////////////////////////
//#include "SurfaceReaction.h"
/////////////////////////
using namespace std;
using namespace MicroProcesses;

int main( int argc, char* argv[] )
{
    //SurfaceReaction sr;
    //string  writeOut;
    //
    //sr.setMessage("Hello I am MySurfaceReaction class");
    //writeOut = sr.getMessage();

    //cout << "Class Output: " << endl;
    //cout  << writeOut << endl;

    Apothesis* apothesis = new Apothesis( argc, argv );
    apothesis->init();
    apothesis->exec();

    if ( apothesis )
      delete apothesis;
}


