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
#include<time.h>

using namespace std;
using namespace MicroProcesses;

#define CLK CLOCK_MONOTONIC
/* Function to compute the difference between two points in time */
struct timespec diff(struct timespec start, struct timespec end){
	struct timespec temp;
	if((end.tv_nsec-start.tv_nsec)<0){
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	}
	else{
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return temp;
}

int main( int argc, char* argv[] )
{
    struct timespec start_e2e, end_e2e, e2e;
	  clock_gettime(CLK, &start_e2e);

    //Checking commit
    Apothesis* apothesis = new Apothesis( argc, argv );

    cout << "Initiating Apothesis" << endl;
    apothesis->init();

    cout << "Apothesis runnning ..." << endl;
    apothesis->exec();
    cout << "Apothesis finished succesfully." << endl;

    clock_gettime(CLK, &end_e2e);
	  e2e = diff(start_e2e, end_e2e);
	  cout << "Time taken for the simulation is : ";
    cout << e2e.tv_sec<<" seconds "<< endl;

    if ( apothesis )
      delete apothesis;
}


