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

#include "pointers.h"
#include "apothesis.h"
#include "lattice/lattice.h"
#include "FCC.h"
#include "io.h"
#include "read.h"
#include "errorhandler.h"
#include "parameters.h"
#include "properties.h"
#include "process.h"
#include "species.h"
#include "string.h"
#include "reaction_new.h"
#include "aux/random_generator.h"
#include <bits/stdc++.h>

#include "factory_process.h"

#include <numeric>

using namespace MicroProcesses;

typedef rapidjson::Document Document;
typedef rapidjson::Value Value;
typedef rapidjson::SizeType SizeType;

//using namespace Utils;

Apothesis::Apothesis(int argc, char *argv[])
    : pLattice(0),
      pRead(0),
      m_dProcTime(0.0),
      m_dRTot(0.0),
      m_dProcRate(0.0),
      m_debugMode(false)
{
    m_iArgc = argc;
    m_vcArgv = argv;

    pParameters = new Utils::Parameters(this);
    pProperties = new Utils::Properties(this);
    pRandomGen = new RandomGen::RandomGenerator( this );

    /* This must be constructed before the input */

    // Create input instance
    pIO = new IO(this);
    pRead = new Read(this);

    vector<string> pName = pRead->getSpeciesNames();

    // Build the lattice. This should always follow the read input
    std::cout << "Building the lattice" << std::endl;
    pLattice->setOrientation("110");
    pLattice->build();

    //For building with steps surface. Works only for simple cubic
    //pLattice->buildSteps( 20, 1, 0);

    std::cout << "Finished building the lattice" << std::endl;

    // initialize number of species
    m_nSpecies = 0;
}

Apothesis::~Apothesis()
{
    delete pIO;
    delete pRead;
    delete pLattice;
}

void Apothesis::init()
{
    cout << "Opening output file" << endl;
    /// The output file name will come from the user and will have the extenstion .log
    /// This would come as a parameter from the user from the args (also the input).
    /// Now both are hard copied.

    if (!pIO->outputOpen())
        pIO->openOutputFile("Output");

    // Initialize Random generator.
    pRandomGen->init( 213212 );

    //To Deifilia: This must be created for each process in order to pass
    //the parameters from the input file to the porcess
    map<string, any> params;
    params.insert( {"T", 500.} );
    params.insert( {"P", 101325.} );
    params.insert( {"f", 2.0e-3} );
    params.insert( {"C_tot", 1.0e+19} );
    params.insert( {"s0", 0.1} );
    params.insert( {"k", 1.3806503e-23} );
    params.insert( {"Na", 6.0221417930e+23} );

    double Ed = 7.14e+4/6.0221417930e+23;
    params.insert( {"E_d", Ed } );
    params.insert( {"Na", 6.0221417930e+23} );

    double Em = 4.28e+4/6.0221417930e+23;
    params.insert( {"E_m", Em } );

    set< Site* > emptySet;
    //    set< Site* > tempSet;
    //    for ( Site* s:pLattice->getSites() )
    //       tempSet.insert( s );

    auto pos = m_processMap.insert( { FactoryProcess::createProcess("AdsorptionFCC1102SSimple"), emptySet } );
    pos.first->first->setName("AdsorptionFCC1102SSimple");
    pos.first->first->init( params );
    pos.first->first->setLattice( pLattice );
    pos.first->first->setRandomGen( pRandomGen );
    pos.first->first->init( params );


/*    for (int i = 4; i < 9; i++){
        auto des = m_processMap.insert( { FactoryProcess::createProcess("DesorptionFCC110Simple"), emptySet } );
        string name = "Desortpion " + std::to_string( i );
        des.first->first->setName( name );
        params.insert( {"neighs", i } );
        des.first->first->init( params );
        params.erase( "neighs" );
    }*/

    for ( auto &p:m_processMap){
        for ( Site* s:pLattice->getSites() ){
            if ( p.first->rules( s ) )
                p.second.insert( s );
        }
    }

//    pLattice->print();
    cout << " Lets see! " << endl;
}

void Apothesis::exec()
{
    Site* tempSite = 0;
    //--------------- Open files for writting ---------------------->
    pIO->openRoughnessFile( "testRough" );
    pIO->writeLogOutput("Running " + to_string( m_dEndTime ) + " sec");

    //Calculate first time the total probability (R) --------------------------//
    m_dRTot = 0.0;
    for (pair<Process*, set< Site* > > p:m_processMap)
        m_dRTot += p.first->getProbability()*(double)p.second.size();

    m_dEndTime = 0.01;

    pIO->writeInOutput( "\n" );
    pIO->writeInOutput( "********************************************************************" );
    string output = "Time (s)"s + '\t' + "Micro roughness (-)" + '\t' + "RMS (-)" + '\t' ;

    for ( auto &p:m_processMap)
        output += p.first->getName() + '\t';

    for ( auto &p:m_processMap)
        output += "Coverage " + p.first->getName() + '\t';
    pIO->writeInOutput( output );

    int iTimeStep = 0;
    pIO->writeLatticeHeights( m_dProcTime, iTimeStep );

    double writeLatHeigsEvery = 1e-3; //in s
    double timeToWrite = 0.0;

    output = std::to_string(m_dProcTime) + '\t' + std::to_string( pProperties->getMicroroughness() ) + '\t' + std::to_string( pProperties->getRMS() )  + '\t' ;
    for ( auto &p:m_processMap)
        output += std::to_string( p.first->getNumEventHappened() ) + '\t';

    for ( auto &p:m_processMap)
        output += std::to_string( (double)p.second.size()/(double)pLattice->getSize() ) + '\t';
    pIO->writeInOutput( output );


    pLattice->writeXYZ( "lattice.xyz" );

    while ( m_dProcTime <= m_dEndTime ){
        //1. Get a random numbers
        m_dSum = 0.0;
        m_dRandom = pRandomGen->getDoubleRandom();

        for ( auto &p:m_processMap){
            m_dProcRate = p.first->getProbability()*p.second.size();
            m_dSum += m_dProcRate/m_dRTot;

            //2. Pick a process according to the rates
            if ( m_dRandom <= m_dSum ){
                //Get a random number which is the ID of the site where this process can performed
                m_iSiteNum = pRandomGen->getIntRandom(0, p.second.size() - 1 );

                //3. From this process pick the random site with id and perform it:
                Site* s = *next( p.second.begin(), m_iSiteNum );

                p.first->perform( s );
                tempSite = s;
                //Count the event for this class
                p.first->eventHappened();

                // Check if an affected site must enter to a class or not
                for (Site* affectedSite:p.first->getAffectedSites() ){
                    //Erase the affected site from the processes
                    for (auto &p2:m_processMap){
                        if ( !p2.first->isUncoAccepted() ) {
                            //Added if it obeys the rules of this process
                            if ( p2.first->rules( affectedSite ) ) {
                                if (p2.second.find( affectedSite ) == p2.second.end() )
                                    p2.second.insert( affectedSite );
                            }
                            else
                                p2.second.erase( affectedSite );


                            if ( p2.second.empty() && p2.first->getName() == "AdsorptionFCC1102SSimple") {
                                cout << "No place to adsorb! " << endl;
                                exit(0);
                            }

                        }
                    }
                }


                //4. Re-compute the processes rates and re-compute Rtot (see ppt).
                m_dRTot = 0.0;
                for (pair<Process*, set< Site* > > p3:m_processMap)
                    m_dRTot += p3.first->getProbability()*p3.second.size();

                //5. Compute dt = -ln(ksi)/Rtot
                m_dt = -log( pRandomGen->getDoubleRandom()  )/m_dRTot;
                break;
            }
        }

        //6. advance time: time += dt;
        m_dProcTime += m_dt;
        timeToWrite += m_dt;

        cout << m_dProcTime << endl;

        pLattice->writeXYZ( "test.xzy" );

        //Write the lattice heights
        iTimeStep++;

//        cout << "******************" << endl;
//        cout <<" Site " << tempSite->getID() << endl;
//        cout << "******************" << endl;
//        pLattice->print();
//        cout << "******************" << endl;

        if ( timeToWrite >= writeLatHeigsEvery ) {
            pIO->writeLatticeHeights( m_dProcTime, iTimeStep );

            output = std::to_string(m_dProcTime) + '\t' + std::to_string( pProperties->getMicroroughness() ) + '\t' + std::to_string( pProperties->getRMS() )  + '\t' ;
            for ( auto &p:m_processMap)
                output += std::to_string( p.first->getNumEventHappened() ) + '\t';

            for ( auto &p:m_processMap)
                output += std::to_string( (double)p.second.size() ) + '\t';
            pIO->writeInOutput( output );

            timeToWrite = 0.0;
        }
    }
}

void Apothesis::logSuccessfulRead(bool read, string parameter)
{
    if (!pIO->outputOpen())
    {
        pIO->openOutputFile("Output");
    }

    read ? pIO->writeLogOutput("Reading " + parameter)
         : pErrorHandler->error_simple_msg("No " + parameter + " found in input file");
}

map<string, Species *> Apothesis::getAllSpecies()
{
    return m_species;
}

Species *Apothesis::getSpecies(string species)
{
    return m_species[species];
}


