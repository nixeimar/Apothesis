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
#include "reaction.h"

#include "factory_process.h"

#include <numeric>
#include <algorithm>

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
    pLattice->setLabels( pParameters->getLatticeLabels() );
    pLattice->build();

    // TODO: Here we must take into account the case of two or more species participating in the film growth
    // and the user should give the per cent of each species in t=0s e.g. 0.8Ga 0.2As
    for ( Site* s:pLattice->getSites() ){
        s->setLabel(  pParameters->getLatticeLabels() );
        s->setBelowLabel(  pParameters->getLatticeLabels() );
    }

    if ( pLattice->hasSteps() )
        pLattice->buildSteps();

    //Print lattice info: To be move in debug version
    pLattice->printInfo();

    pLattice->print();

    //Print parameters to check: To be move in debug version
    pParameters->printInfo();

    //An empty set used for the initialization of the processMap
    set< Site* > emptySet;

    //Create the processes
    for ( auto proc:pParameters->getProcessesInfo() ){

        string process = m_fAnalyzeProc( proc.first );

        map<string, double> reactants;
        for (string react: pIO->getReactants( proc.first ) )
            reactants.insert( pIO->analyzeCompound( react ) );

        map<string, double> products;
        for (string react: pIO->getProducts( proc.first ) )
            products.insert( pIO->analyzeCompound( react ) );


        if ( process.compare("Adsorption") == 0 ){

            Adsorption* a = new Adsorption();
            //The name of the actual class used

            for ( pair<string, int> s: reactants) {
                if ( s.first.compare("*") != 0 )
                    a->setAdrorbed( s.first );
                else
                    a->setNumSites( s.second );
            }

            a->setName( proc.first );
            a->setLattice( pLattice );
            a->setRandomGen( pRandomGen );
            a->setErrorHandler( pErrorHandler );
            a->setSysParams( pParameters ); //These are the systems and constants parameters
            a->init( proc.second ); //These are the process per se parameters

            m_processMap.insert( {a, emptySet} );
        }
        else if ( process.compare("Reaction") == 0 ){

            Reaction* r = new Reaction();
            //The name of the actual class used
            r->setName( proc.first );
            r->setLattice( pLattice );
            r->setRandomGen( pRandomGen );
            r->setErrorHandler( pErrorHandler );
            r->setSysParams( pParameters ); //These are the systems and constants parameters
            r->init( proc.second ); //These are the process per se parameters

            m_processMap.insert( {r, emptySet} );
        }
        else if ( process.compare("Desorption") == 0 ){

            if (proc.second.at( proc.second.size() - 1 ).compare("all") != 0){

                proc.second.push_back( to_string(1) );

                Desorption* des = new Desorption();
                des->setName( proc.first );
                des->setLattice( pLattice );
                des->setRandomGen( pRandomGen );
                des->setErrorHandler( pErrorHandler );
                des->setSysParams( pParameters ); //These are the systems and constants parameters
                des->init( proc.second ); //These are the process per se parameters

                m_processMap.insert( {des, emptySet} );

            } else {
                for ( int neighs = 0; neighs < pLattice->getNumFirstNeihgs(); neighs++) {

                    proc.second.pop_back();
                    proc.second.push_back( to_string(neighs + 1) );

                    Desorption* des = new Desorption();
                    des->setName( proc.first );
                    des->setLattice( pLattice );
                    des->setRandomGen( pRandomGen );
                    des->setErrorHandler( pErrorHandler );
                    des->setSysParams( pParameters ); //These are the systems and constants parameters
                    des->init( proc.second ); //These are the process per se parameters

                    m_processMap.insert( {des, emptySet} );
                }
            }
        }
        else if ( process.compare("Diffusion") == 0 ){

            if (proc.second.at( proc.second.size() - 1 ).compare("all") != 0){

                proc.second.push_back( to_string(1) );

                Diffusion* dif = new Diffusion();
                dif->setName( proc.first );
                dif->setLattice( pLattice );
                dif->setRandomGen( pRandomGen );
                dif->setErrorHandler( pErrorHandler );
                dif->setSysParams( pParameters ); //These are the systems and constants parameters
                dif->init( proc.second ); //These are the process per se parameters

                m_processMap.insert( {dif, emptySet} );

            } else {

                for ( int neighs = 0; neighs < pLattice->getNumFirstNeihgs(); neighs++) {

                    proc.second.pop_back();
                    proc.second.push_back( to_string(neighs + 1) );

                    Diffusion* dif = new Diffusion();
                    dif->setName( proc.first + "( " + to_string(neighs + 1) + " )"  );
                    dif->setLattice( pLattice );
                    dif->setRandomGen( pRandomGen );
                    dif->setErrorHandler( pErrorHandler );
                    dif->setSysParams( pParameters ); //These are the systems and constants parameters
                    dif->init( proc.second ); //These are the process per se parameters

                    m_processMap.insert( {dif, emptySet} );
                }
            }
        }

        /*        vector<string> params;
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
        }*/
    }

    //Partition the lattice sites depending on the rules of each process
    for ( auto &p:m_processMap ){
        for ( Site* s:pLattice->getSites() ){
            if ( p.first->rules( s ) )
                p.second.insert( s );
        }
    }

    //The end time of the simulation
    m_dEndTime = pParameters->getEndTime();

    //Calculate first time the total probability (R) for apothesis to start --------------------------//
    m_dRTot = 0.0;
    for (pair<Process*, set< Site* > > p:m_processMap)
        m_dRTot += p.first->getProbability()*(double)p.second.size();

    //Start writing in the output log
    //Write initialization info to log
    pIO->writeLogOutput("Apothesis build on " __TIMESTAMP__);
    pIO->writeLogOutput("-------------------------------------------------");
    pIO->writeLogOutput("");
    pIO->writeLogOutput("End time " + to_string( m_dEndTime ) + " sec");
    pIO->writeLogOutput("Temperature " + to_string( pParameters->getTemperature() ) + " K");
    pIO->writeLogOutput("Pressure " + to_string( pParameters->getPressure() ) + " P");
    pIO->writeLogOutput("Random init num " + to_string( pParameters->getRandGenInit() ) );

    string toWrite = "\n";
    toWrite = "Lattice " +  pLattice->getTypeAsString() + " ";
    toWrite += to_string( pLattice->getX() ) + " ";
    toWrite += to_string( pLattice->getY() ) + " ";

    if ( pLattice->hasSteps() ) {
        toWrite += "stepped  ";
        toWrite += to_string( pLattice->getNumSteps() ) + " ";
        toWrite += to_string( pLattice->getStepHeight() ) + " ";
    }
    pIO->writeInOutput( toWrite );


    pIO->writeInOutput(" ");
    pIO->writeLogOutput("Processes");
    for (auto proc:pParameters->getProcessesInfo() ) {
        toWrite = "";
        toWrite = proc.first + " ";
        for ( string str:proc.second ) {
            toWrite += str + " ";
        }
        pIO->writeLogOutput( toWrite );
    }

    pIO->writeInOutput( "\n" );
    pIO->writeInOutput( "********************************************************************" );

    string output = "Time (s)"s + '\t' + "Growth rate (ML/s)" + '\t' + "RMS (-)" + '\t' + "Micro-roughness (-)" + '\t';

    for ( auto &p:m_processMap)
        output += p.first->getName() + '\t';

    for ( auto &p:m_processMap)
        output +=  p.first->getName() + " (class size)" + '\t';

    pIO->writeInOutput( output );
    pIO->writeLatticeHeights( m_dProcTime );
}

void Apothesis::exec()
{
    double timeToWriteLog = 0;
    double timeToWriteLattice = 0;

    string output ="";

    //    pLattice->writeXYZ( "initial.xzy" );

    // The average height for the first time
    double timeGrowth = 0;
    double meanDHPrevStep = pProperties->getMeanDH();
    output = std::to_string(m_dProcTime) + '\t'
            + std::to_string( 0.0  ) + '\t'
            + std::to_string( pProperties->getRMS() )  + '\t'
            + std::to_string( pProperties->getMicroroughness() )  + '\t';

    for ( auto &p:m_processMap)
        output += std::to_string( p.first->getNumEventHappened() ) + '\t';

    for ( auto &p:m_processMap)
        output += std::to_string( p.second.size() ) + '\t';

    pIO->writeInOutput( output );

    while ( m_dProcTime <= m_dEndTime ){
        //1. Get a random numbers
        m_dSum = 0.0;
        m_dRandom = pRandomGen->getDoubleRandom();

        for ( auto &p:m_processMap){
            m_dProcRate = p.first->getProbability()*(double)p.second.size();
            m_dSum += m_dProcRate/m_dRTot;

            //2. Pick a process according to the rates
            if ( m_dRandom <= m_dSum ){

                // Calculate the average Height before
                //                aveDH1 = pProperties->getMeanDH();

                //Get a random number which is the ID of the site where this process can performed
                m_iSiteNum = pRandomGen->getIntRandom(0, p.second.size() - 1 );

                //3. From this process pick the random site with id and perform it:
                Site* s = *next( p.second.begin(), m_iSiteNum );

                //Compute the average height before performing the process to measure the growth rate
                meanDHPrevStep = pProperties->getMeanDH();
                timeGrowth = m_dProcTime;

                p.first->perform( s );

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


        double test = (pProperties->getMeanDH() - meanDHPrevStep) / ((pLattice->getSize()*(m_dProcTime - timeGrowth) ) );

        if ( timeToWriteLog >= pParameters->getWriteLogTimeStep() ){
            output = std::to_string(m_dProcTime) + '\t'
                    + std::to_string( (pProperties->getMeanDH() - meanDHPrevStep) / ( (pLattice->getSize()*(m_dProcTime - timeGrowth) ) ) )+ '\t'
                    + std::to_string( pProperties->getRMS() )  + '\t'
                    + std::to_string( pProperties->getMicroroughness() )  + '\t';

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

    pIO->writeLatticeHeights( m_dProcTime  );
    output = std::to_string(m_dProcTime) + '\t'
            + std::to_string( (pProperties->getMeanDH() - meanDHPrevStep)/ (pLattice->getSize())*(m_dProcTime - timeGrowth)  ) + '\t'
            + std::to_string( pProperties->getRMS() )  + '\t'
            + std::to_string( pProperties->getMicroroughness() )  + '\t';

    for ( auto &p:m_processMap)
        output += std::to_string( p.first->getNumEventHappened() ) + '\t';

    for ( auto &p:m_processMap)
        output += std::to_string( p.second.size() ) + '\t';

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

string Apothesis::m_fAnalyzeProc(string process){

    vector<string> parts = pIO->split( process, "->");

    //It is reaction or adsoprtion
    if ( pIO->contains( parts[0], "+" ) ){

        vector<string> reactants = pIO->split( parts[ 0 ], "+" );
        for (string s:reactants){
            pIO->trim(s);
            if (  s.compare("*") == 0 )
                return "Adsorption";
        }

        return "Reaction";
    }
    else {

        vector<string> products = pIO->split( parts[ 1 ], "+" );
        for (string s:products){
            pIO->trim(s);
            if (  s.compare("*") == 0 )
                return "Desorption";
        }

        return "Diffusion";
    }
}
