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
    delete pParameters;
    delete pErrorHandler;
    delete pRandomGen;
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

    //Create the lattice
    pLattice->build();
    if ( pLattice->hasSteps() )
        pLattice->buildSteps();
    pLattice->printInfo();
    pLattice->print();

    //Prin parameters if you want
    pParameters->printInfo();

    //An empty set used for the initialization of the processMap
    set< Site* > emptySet;

    //Create the processes
    for ( auto proc:pParameters->getProcessesInfo() ){
        vector<string> params;
        if ( proc.first.compare( "Desorption" ) == 0 || proc.first.compare( "Diffusion" ) == 0 ){
            for ( int neighs = 0; neighs < pLattice->getNumFirstNeihgs(); neighs++) {
                auto pos = m_processMap.insert( { FactoryProcess::createProcess(proc.first), emptySet } );
                params = proc.second;
                params.push_back( to_string(neighs + 1) ); //Since neighbors start from 1 and not 0
                pos.first->first->setName( proc.first + " N:" +  to_string( neighs + 1 ) ); //The name of the actuall class used
                pos.first->first->setLattice( pLattice );
                pos.first->first->setRandomGen( pRandomGen );
                pos.first->first->setSysParams( pParameters );
                pos.first->first->setErrorHandler( pErrorHandler );

                pos.first->first->init( params );
            }
        }
        else {
            auto pos = m_processMap.insert( { FactoryProcess::createProcess(proc.first), emptySet } );
            pos.first->first->setName( proc.first ); //The name of the actual class used
            pos.first->first->setLattice( pLattice );
            pos.first->first->setRandomGen( pRandomGen );
            pos.first->first->setErrorHandler( pErrorHandler );
            pos.first->first->setSysParams( pParameters );

            //Check how we will impose that
            pos.first->first->setUncoAccepted(true);

            pos.first->first->init( proc.second );
        }
    }

    for ( auto &p:m_processMap ){
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
    //pIO->openRoughnessFile( "testRough" );
    pIO->writeLogOutput("Running " + to_string( m_dEndTime ) + " sec");

    //Calculate first time the total probability (R) --------------------------//
    m_dRTot = 0.0;
    for (pair<Process*, set< Site* > > p:m_processMap)
        m_dRTot += p.first->getProbability()*(double)p.second.size();

    m_dEndTime = pParameters->getEndTime();

    pIO->writeInOutput( "\n" );
    pIO->writeInOutput( "********************************************************************" );
    string output = "Time (s)"s + '\t' +     "Growth rate (ML/s)" + '\t' + "RMS (-)" + '\t' ;

    for ( auto &p:m_processMap)
        output += p.first->getName() + '\t';

    for ( auto &p:m_processMap)
        output +=  p.first->getName() + " (class size)" + '\t';

    pIO->writeInOutput( output );
    pIO->writeLatticeHeights( m_dProcTime );

    double timeToWriteLog = 0;
    double timeToWriteLattice = 0;

//    pLattice->writeXYZ( "initial.xzy" );

    // The average height for the first time
    double aveDH1 = pProperties->getMeanDH();

    output = std::to_string(m_dProcTime) + '\t' + std::to_string( fabs(pProperties->getMeanDH() - aveDH1) ) + '\t' + std::to_string( pProperties->getRMS() )  + '\t' ;
    for ( auto &p:m_processMap)
        output += std::to_string( p.first->getNumEventHappened() ) + '\t';

    for ( auto &p:m_processMap)
        output += std::to_string( p.second.size() ) + '\t';

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

//                cout << "Performing: " << p.first->getName() << endl;

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

                //To print additional info
        //        for (auto &p2:m_processMap)
        //            cout << p2.first->getName() << " " << p2.second.size() << endl;

                //5. Compute dt = -ln(ksi)/Rtot
                m_dt = -log( pRandomGen->getDoubleRandom()  )/m_dRTot;
                break;
            }
        }

        //6. advance time: time += dt;
        m_dProcTime += m_dt;

        //Here compute the time for writing
        timeToWriteLog += m_dt;
        timeToWriteLattice += m_dt;

        cout << "Time: " << m_dProcTime << endl;

        if ( timeToWriteLog >= pParameters->getWriteLogTimeStep() ){ // writeLatHeigsEvery ) {
            output = std::to_string(m_dProcTime) + '\t' + std::to_string( fabs(pProperties->getMeanDH() ) ) + '\t' + std::to_string( pProperties->getRMS() )  + '\t' ;
            for ( auto &p:m_processMap)
                output += std::to_string( p.first->getNumEventHappened() ) + '\t';

            for ( auto &p:m_processMap)
                output += std::to_string( p.second.size() ) + '\t';

            pIO->writeInOutput( output );
            timeToWriteLog = 0.0;
        }

        if ( timeToWriteLattice >= pParameters->getWriteLatticeTimeStep() ) {
            string latName = "lattice_" + to_string( m_dProcTime )+".xyz";
            pIO->writeLatticeHeights( m_dProcTime  );
            timeToWriteLattice = 0.0;
        }
    }

    string latName = "lattice_" + to_string( m_dProcTime )+".xyz";
    //pLattice->writeXYZ( latName );
    pIO->writeLatticeHeights( m_dProcTime  );

    output = std::to_string(m_dProcTime) + '\t' + std::to_string( fabs(pProperties->getMeanDH() ) ) + '\t' + std::to_string( pProperties->getRMS() )  + '\t' ;
    for ( auto &p:m_processMap)
        output += std::to_string( p.first->getNumEventHappened() ) + '\t';

    for ( auto &p:m_processMap)
        output += std::to_string( p.second.size() ) + '\t';

//    for ( auto &p:m_processMap)
//        output += std::to_string( (double)p.first->getProbability() ) + '\t';
    pIO->writeInOutput( output );

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
