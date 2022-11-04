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

#include "pointers.h"
#include "apothesis.h"
#include "lattice/lattice.h"
#include "FCC.h"
#include "io.h"
#include "errorhandler.h"
#include "parameters.h"
#include "properties.h"
#include "process.h"
#include "species.h"
#include "string.h"
#include "reaction.h"
#include "aux/random_generator.h"
#include <bits/stdc++.h>
#include "reader.h"

#include "factory_process.h"

#include <numeric>

using namespace MicroProcesses;

//using namespace Utils;

Apothesis::Apothesis(int argc, char *argv[])
    : pLattice(0),
      pReader(0),
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
    pLattice = new SimpleCubic(this);

    // Create input instance
    pIO = new IO(this);
    pIO->init(m_iArgc, m_vcArgv);

    // initialize number of species
    m_nSpecies = 0;
}

Apothesis::~Apothesis()
{
    delete pIO;
    delete pReader;
    delete pLattice;
}

void Apothesis::init()
{
    //Read the input file
    pIO->readInputFile();

    //Open the output file
    if ( !pIO->outputOpen() )
        pIO->openOutputFile("Output");

    // Initialize Random generator
    pRandomGen->init( pParameters->getRandGenInit() );

    //Create the processes or get the processes

    //To Deifilia: This must be created for each process in order to pass
    //the parameters from the input file to the porcess
    map<string, any> params;
    params.insert( {"T", 533.} );
    params.insert( {"f", 0.0000340636296296296} ); //Arrhenius
    params.insert( {"s0", 0.137295210293006} ); //Arrhenius

    params.insert( {"P", 1333.22} );
    params.insert( {"C_tot", 2.0e+19} ); //For Cu
    params.insert( {"k", 1.3806503e-23} );
    params.insert( {"Na", 6.0221417930e+23} );

    double Ed = 7.14e+4/6.0221417930e+23;
    params.insert( {"E_d", Ed } );
    params.insert( {"Na", 6.0221417930e+23} );

    double Em = 4.28e+4/6.0221417930e+23;
    params.insert( {"E_m", Em } );

    set< Site* > emptySet;

    //FOr IKY ---->
/*    auto pos = m_processMap.insert( { FactoryProcess::createProcess("AdsroptionSimpleCubic4sMulti"), emptySet } );
    pos.first->first->setName("AdsortpionSimpleCubinc4SMulti");
    pos.first->first->init( params );
    pos.first->first->setLattice( pLattice );
    pos.first->first->setRandomGen( pRandomGen );
    auto des = m_processMap.insert( { FactoryProcess::createProcess("DesorptionSimpleCubic4sMulti"), emptySet } );
    des.first->first->setName("DesorptionSimpleCubic4sMulti");

//    des.first->first->setConstant( true );
    des.first->first->init( params );
    des.first->first->setLattice( pLattice );
    des.first->first->setRandomGen( pRandomGen );
    for ( Site* s:pLattice->getSites() )
        s->setLabel("Cu"); // in all cases we start with Cu
    for ( auto &p:m_processMap){
        cout << p.first->getName() << endl;
        for ( Site* s:pLattice->getSites() ){
            if ( p.first->rules( s ) )
                p.second.insert( s );
        }
    } <---------------- */

    for ( auto &p:m_processMap){
        cout << p.first->getName() << endl;
        for ( Site* s:pLattice->getSites() ){
            if ( p.first->rules( s ) )
                p.second.insert( s );
        }
    }
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

    m_dEndTime = 10;

    pIO->writeInOutput( "\n" );
    pIO->writeInOutput( "********************************************************************" );
    string output = "Time (s)"s + '\t' + "Growth rate (Ang/min)" + '\t' + "RMS (-)" + '\t' ;

    for ( auto &p:m_processMap)
        output += p.first->getName() + '\t';

    for ( auto &p:m_processMap)
        output += "Coverage " + p.first->getName() + '\t';
    pIO->writeInOutput( output );

    int iTimeStep = 0;
    pIO->writeLatticeHeights( m_dProcTime, iTimeStep );

    double writeLatHeigsEvery = 0.001; //in s
    double timeToWrite = 0;

    pLattice->writeXYZ( "initial.xzy" );

    // The average height for the first time
    double aveDH1 = pProperties->getMeanDH();

    output = std::to_string(m_dProcTime) + '\t' + std::to_string( fabs(pProperties->getMeanDH() - aveDH1)*2.55*60 ) + '\t' + std::to_string( pProperties->getRMS() )  + '\t' ;
    for ( auto &p:m_processMap)
        output += std::to_string( p.first->getNumEventHappened() ) + '\t';

    for ( auto &p:m_processMap)
        output += std::to_string( (double)p.second.size()/(double)pLattice->getSize() ) + '\t';
    pIO->writeInOutput( output );

    double oldTime = 0.0;
    double oldAve = 0.0;

    while ( m_dProcTime <= m_dEndTime ){
        //1. Get a random numbers
        m_dSum = 0.0;
        m_dRandom = pRandomGen->getDoubleRandom();
        oldTime = m_dProcTime;
        oldAve = pProperties->getMeanDH();
        for ( auto &p:m_processMap){
            if ( !p.first->isConstant() )
                m_dProcRate = p.first->getProbability()*(double)p.second.size();
            else {
                if ( !p.second.empty() )
                    m_dProcRate = p.first->getProbability();
            }

            m_dSum += m_dProcRate/m_dRTot;

            //2. Pick a process according to the rates
            if ( m_dRandom <= m_dSum ){

                // Calculate the average Height before the
                aveDH1 = pProperties->getMeanDH();

                //Get a random number which is the ID of the site where this process can performed
                m_iSiteNum = pRandomGen->getIntRandom(0, p.second.size() - 1 );

                //3. From this process pick the random site with id and perform it:
                Site* s = *next( p.second.begin(), m_iSiteNum );

                cout << "Performing: " << p.first->getName() << endl;

                p.first->perform( s );
                tempSite = s;
                //Count the event for this class
                p.first->eventHappened();

                // Check if an affected site must enter tob a class or not
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
                        }
                    }
                }

                //4. Re-compute the processes rates and re-compute Rtot (see ppt).
                m_dRTot = 0.0;
                for (pair<Process*, set< Site* > > p3:m_processMap){
                    if ( p3.first->isConstant() ){
                        if ( p3.second.size() != 0 )
                            m_dRTot += p3.first->getProbability();
                    }
                    else
                        m_dRTot += p3.first->getProbability()*(double)p3.second.size();
                }

               pLattice->print();

               for (auto &p2:m_processMap)
                   cout << p2.first->getName() << " " << p2.second.size() << endl;


                cout << "Timestep: " << iTimeStep << endl;

                //5. Compute dt = -ln(ksi)/Rtot
                m_dt = -log( pRandomGen->getDoubleRandom()  )/m_dRTot;
                break;
            }
        }

        //6. advance time: time += dt;
        m_dProcTime += m_dt;
        timeToWrite += m_dt;

        cout << m_dProcTime << endl;

 //       for (auto &p2:m_processMap)
 //           cout << p2.first->getName() << " " << p2.second.size() << endl;


        //Write the lattice heights
        iTimeStep++;

        if ( timeToWrite >= writeLatHeigsEvery ) {
            string latName = "lattice_" + to_string( m_dProcTime )+".xyz";
 //           pLattice->writeXYZ( latName );
 //          pLattice->writeLatticeHeights( m_dProcTime, iTimeStep );

            output = std::to_string(m_dProcTime) + '\t' + std::to_string( 2.55*fabs(pProperties->getMeanDH() ) ) + '\t' + std::to_string( pProperties->getRMS() )  + '\t' ;
            for ( auto &p:m_processMap)
                output += std::to_string( p.first->getNumEventHappened() ) + '\t';

            for ( auto &p:m_processMap)
                output += p.first->getName() + '\t' + std::to_string( (double)p.first->getProbability()*p.second.size() ) + '\t';
            pIO->writeInOutput( output );

            timeToWrite = 0.0;
        }
    }

    string latName = "lattice_" + to_string( m_dProcTime )+".xyz";
    pLattice->writeXYZ( latName );

    pIO->writeLatticeHeights( m_dProcTime, iTimeStep );

    //    output = std::to_string(m_dProcTime) + '\t' + std::to_string( 2.55*60.*fabs(pProperties->getMeanDH() - aveDH1)/m_dt ) + '\t' + std::to_string( pProperties->getRMS() )  + '\t' ;
    output = std::to_string(m_dProcTime) + '\t' + std::to_string( 2.55*fabs(pProperties->getMeanDH() ) ) + '\t' + std::to_string( pProperties->getRMS() )  + '\t' ;
    for ( auto &p:m_processMap)
        output += std::to_string( p.first->getNumEventHappened() ) + '\t';

    for ( auto &p:m_processMap)
        output += p.first->getName() + '\t' + std::to_string( (double)p.first->getProbability() ) + '\t';
    pIO->writeInOutput( output );

    //    pLattice->print();
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
