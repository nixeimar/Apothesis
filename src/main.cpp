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

#include <iostream>
#include <list>
#include "site.h"
#include "lattice.h"
#include "process.h"
#include "apothesis.h"
#include "io.h"

#include "catalysis.h"
#include "pvd.h"
#include "cvd.h"
#include "etching.h"
#include "ald.h"
#include "ale.h"


using namespace std;
using namespace MicroProcesses;

int main( int argc, char* argv[] )
{
    //The primary communication channel is reading from a file
    IO* io =  new IO();

    //Pass any arguments given by the user
    io->init(argc, argv);

    //Read the input file and proceed accordingly based on the process read
    io->readInputFile();

    if ( io->getParameters()->sProcess().compare("catalysis") == 0 ) {
        Catalysis* catalysis = new Catalysis();
        catalysis->setIO( io );
        catalysis->init();
        catalysis->perform();
    }
    else if (io->getParameters()->sProcess().compare("PVD") == 0 || io->getParameters()->sProcess().compare("pvd") == 0)
    {
        PVD* pvd = new PVD();
        pvd->setIO( io );
        pvd->init();
        pvd->perform();

    }
    else if (io->getParameters()->sProcess().compare("CVD") == 0 || io->getParameters()->sProcess().compare("cvd") == 0)
    {
        CVD* cvd = new CVD();
        cvd->setIO( io );
        cvd ->init();
        cvd ->perform();

    }
    else if (io->getParameters()->sProcess().compare("etching") == 0)
    {
        Etching* etching = new Etching();
        etching->setIO( io );
        etching->init();
        etching->perform();

    }
    else if (io->getParameters()->sProcess().compare("ALD") == 0 || io->getParameters()->sProcess().compare("ald") == 0)
    {
        ALD* ald = new ALD();
        ald->setIO( io );
        ald->init();
        ald->perform();

    }
    else if (io->getParameters()->sProcess().compare("ALE") == 0 || io->getParameters()->sProcess().compare("ale") == 0)
    {
        ALE* ale = new ALE();
        ale->setIO( io );
        ale->init();
        ale->perform();
    }
    else {
        cout << "Error: Unknown or nor defined process. Please make sure that the input file contains the \"process\" keyword. "
                "Available processes are: catalysis, PVD, CVD, etching, ALD and ALE" << endl;
    }

}


