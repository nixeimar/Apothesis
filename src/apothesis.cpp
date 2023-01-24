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
        s->setOccupied( false ); //Start from clear surface
    }

    if ( pLattice->hasSteps() )
        pLattice->buildSteps();

    //Print lattice info: To be move in debug version
    pLattice->printInfo();

    //pLattice->print();

    //Print parameters to check: To be move in debug version
    pParameters->printInfo();

    //An empty set used for the initialization of the processMap
    set< Site* > emptySet;

    //Create the processes
    for ( auto proc:pParameters->getProcessesInfo() ){

        string process = mf_analyzeProc( proc.first );

        if ( process.compare("Adsorption") == 0 ){

            unordered_map<string, int> reactants;
            for (string react: pIO->getReactants( proc.first ) )
                reactants.insert( pIO->analyzeCompound( react ) );

            unordered_map<string, int> products;
            for (string prod: pIO->getProducts( proc.first ) )
                products.insert( pIO->analyzeCompound( prod ) );


            if (proc.second.at( proc.second.size() - 1 ).compare("all") != 0){

                Adsorption* a = new Adsorption();
                for ( pair<string, int> s: products) {
                    a->setAdrorbed( s.first );
                    a->setNumSites( s.second );
                }

                a->setName( proc.first );
                a->setLattice( pLattice );
                a->setRandomGen( pRandomGen );
                a->setErrorHandler( pErrorHandler );
                a->setSysParams( pParameters ); //These are the systems and constants parameters
                a->init( proc.second ); //These are the process per se parameters

                m_processMap.insert( {a, emptySet} );

            } else {

                for ( int neighs = 0; neighs < pLattice->getNumFirstNeihgs(); neighs++) {

                    proc.second.pop_back();
                    proc.second.push_back( to_string(neighs) );

                    Adsorption* a = new Adsorption();

                    for ( pair<string, int> s: products) {
                        a->setAdrorbed( s.first );
                        a->setNumSites( s.second );
                    }

                    a->setAllNeighs(true);
                    a->setName( proc.first + " (" + to_string(neighs) + " N)" );
                    a->setLattice( pLattice );
                    a->setRandomGen( pRandomGen );
                    a->setErrorHandler( pErrorHandler );
                    a->setSysParams( pParameters ); //These are the systems and constants parameters
                    a->init( proc.second ); //These are the process per se parameters

                    m_processMap.insert( {a, emptySet} );
                }
            }
        }
        else if ( process.compare("Reaction") == 0 ){

            vector<string> reactants;
            vector<int> coefReactants;

            for (string react: pIO->getReactants( proc.first ) ){
                reactants.push_back(  pIO->analyzeCompound( react ).first  );
                coefReactants.push_back(  pIO->analyzeCompound( react ).second  );
            }

            vector<string> products;
            vector<int> coefProducts;

            for (string react: pIO->getProducts( proc.first ) ){
                products.push_back(  pIO->analyzeCompound( react ).first  );
                coefProducts.push_back(  pIO->analyzeCompound( react ).second  );
            }

            unordered_map<string, int> reactantsmap;
            for (string react: pIO->getReactants( proc.first ) )
                reactantsmap.insert( pIO->analyzeCompound( react ) );

            unordered_map<string, int> productsmap;
            for (string prod: pIO->getProducts( proc.first ) )
                productsmap.insert( pIO->analyzeCompound( prod ) );

            Reaction* r = new Reaction();
            //The name of the actual class used
            r->setName( proc.first );
            r->setLattice( pLattice );
            r->setRandomGen( pRandomGen );
            r->setErrorHandler( pErrorHandler );
            r->setSysParams( pParameters );
            r->setReactants( reactantsmap );
            r->setReactants( reactants );
            r->setProducts( productsmap);
            r->setProducts( products);
            r->setCoefReactants( coefReactants );
            r->setCoefProducts( coefProducts );

            r->init( proc.second ); //These are the process per se parameters

            m_processMap.insert( {r, emptySet} );
        }
        else if ( process.compare("Desorption") == 0 ){

            unordered_map<string, int> reactants;
            for (string react: pIO->getReactants( proc.first ) )
                reactants.insert( pIO->analyzeCompound( react ) );

            unordered_map<string, int> products;
            for (string prod: pIO->getProducts( proc.first ) )
                products.insert( pIO->analyzeCompound( prod ) );


            if (proc.second.at( proc.second.size() - 1 ).compare("all") != 0){

                proc.second.push_back( to_string(1) );

                Desorption* des = new Desorption();

                for ( pair<string, int> s: products) {
                    if ( s.first.compare("*") != 0 )
                        des->setDesorbed( s.first );
                }

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
                    proc.second.push_back( to_string(neighs) );

                    Desorption* des = new Desorption();

                    for ( pair<string, int> s: products) {
                        if ( s.first.compare("*") != 0 )
                            des->setDesorbed( s.first );
                    }

                    des->setAllNeighs(true);
                    des->setName( proc.first + " (" + to_string(neighs + 1) + " N)" );
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

            unordered_map<string, int> reactants;
            for (string react: pIO->getReactants( proc.first ) )
                reactants.insert( pIO->analyzeCompound( react ) );

            unordered_map<string, int> products;
            for (string prod: pIO->getProducts( proc.first ) )
                products.insert( pIO->analyzeCompound( prod ) );


            if (proc.second.at( proc.second.size() - 1 ).compare("all") != 0){

                proc.second.push_back( to_string(1) );

                Diffusion* dif = new Diffusion();

                for ( pair<string, int> s: products) {
                    if ( s.first.compare("*") != 0 ) {

                        std::string::size_type i = s.first.find("*");
                        if (i != std::string::npos)
                            dif->setDiffused( s.first.erase(i, s.first.length() ) );
                    }
                }

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
                    proc.second.push_back( to_string(neighs) );

                    Diffusion* dif = new Diffusion();

                    for ( pair<string, int> s: products) {
                        if ( s.first.compare("*") != 0 ) {

                            std::string::size_type i = s.first.find("*");
                            if (i != std::string::npos)
                                dif->setDiffused( s.first.erase(i, s.first.length() ) );

                        }
                    }

                    dif->setAllNeighs(true);
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
    }


    /*    pLattice->getSite( 1)->setOccupied(true);
    pLattice->getSite( 1)->setLabel("CO*");
    pLattice->getSite( 2)->setOccupied(true);
    pLattice->getSite( 2)->setLabel("O*");

    pLattice->getSite( 7)->setOccupied(true);
    pLattice->getSite( 7)->setLabel("O*");
    pLattice->getSite( 12)->setOccupied(true);
    pLattice->getSite( 12)->setLabel("CO*");

    pLattice->getSite( 23)->setOccupied(true);
    pLattice->getSite( 23)->setLabel("O*");
    pLattice->getSite( 24)->setOccupied(true);
    pLattice->getSite( 24)->setLabel("O*");


    pLattice->getSite( 18)->setOccupied(true);
    pLattice->getSite( 18)->setLabel("O*");
    pLattice->getSite( 19)->setOccupied(true);
    pLattice->getSite( 19)->setLabel("CO*");*/

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

    m_bHasGrowth = pParameters->getGrowthSpecies().size() > 0 ? true : false;
    m_bReportCoverages = pParameters->getCoverageSpecies().size() > 0 ? true : false;

    // If the user wants the coverages to be reported
    if ( m_bReportCoverages ){
        unordered_map<string, double> covs = pLattice->computeCoverages( pParameters->getCoverageSpecies() );
        for ( auto &p:covs)
            output +=  p.first + " (coverage)" + '\t';
    }

    pIO->writeInOutput( output );

    if ( m_bHasGrowth )
        pIO->writeLatticeHeights( m_dProcTime );

    if ( m_bReportCoverages )
        pIO->writeLatticeSpecies( m_dProcTime  );
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

    ostringstream streamObj;
    streamObj.precision(15);
    streamObj << m_dProcTime;
    //output = std::to_string( m_dProcTime ) + '\t'
    output = streamObj.str() + '\t'
            + std::to_string( 0.0  ) + '\t'
            + std::to_string( pProperties->getRMS() )  + '\t'
            + std::to_string( pProperties->getMicroroughness() )  + '\t';

    for ( auto &p:m_processMap)
        output += std::to_string( p.first->getNumEventHappened() ) + '\t';

    for ( auto &p:m_processMap)
        output += std::to_string( p.second.size() ) + '\t';

    if ( m_bReportCoverages ) {
        unordered_map<string, double> covs = pLattice->computeCoverages( pParameters->getCoverageSpecies() );

        for ( auto &p:covs)
            output += std::to_string( p.second ) + '\t';
    }

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

                //4. Re-compute the processes rates and re-compute Rtot (see ppt)
                m_dRTot = 0.0;
                for (pair<Process*, set< Site* > > p3:m_processMap)
                    m_dRTot += p3.first->getProbability()*(double)p3.second.size();

                //5. Compute dt = -ln(ksi)/Rtot
                m_dt = -log( pRandomGen->getDoubleRandom()  )/m_dRTot;
                //                cout << m_dt << endl;
                break;
            }
        }

        //6. advance time: time += dt;
        m_dProcTime += m_dt;

        //Here compute the time for writing
        timeToWriteLog += m_dt;
        timeToWriteLattice += m_dt;

        if ( timeToWriteLog >= pParameters->getWriteLogTimeStep() ){

            ostringstream streamObj;
            streamObj.precision(15);
            streamObj << m_dProcTime;
            //            output = std::to_string( m_dProcTime ) + '\t'
            output = streamObj.str() + '\t'
                    + std::to_string( (pProperties->getMeanDH() - meanDHPrevStep) / ( (pLattice->getSize()*(m_dProcTime - timeGrowth) ) ) )+ '\t'
                    + std::to_string( pProperties->getRMS() )  + '\t'
                    + std::to_string( pProperties->getMicroroughness() )  + '\t';

            for ( auto &p:m_processMap)
                output += std::to_string( p.first->getNumEventHappened() ) + '\t';

            for ( auto &p:m_processMap)
                output += std::to_string( p.second.size() ) + '\t';

            if ( m_bReportCoverages ) {
                unordered_map<string, double> covs = pLattice->computeCoverages( pParameters->getCoverageSpecies() );

                for ( auto &p:covs)
                    output += std::to_string( p.second ) + '\t';
            }

            pIO->writeInOutput( output );
            timeToWriteLog = 0.0;
        }

        if ( timeToWriteLattice >= pParameters->getWriteLatticeTimeStep() ) {

            if ( m_bHasGrowth )
                pIO->writeLatticeHeights( m_dProcTime );

            if ( m_bReportCoverages )
                pIO->writeLatticeSpecies( m_dProcTime  );

            timeToWriteLattice = 0.0;
        }
    }

    ostringstream streamObjEnd;
    streamObjEnd.precision(15);
    streamObjEnd << m_dProcTime;
    //            output = std::to_string( m_dProcTime ) + '\t'
    output = streamObjEnd.str() + '\t'
            + std::to_string( (pProperties->getMeanDH() - meanDHPrevStep)/ (pLattice->getSize())*(m_dProcTime - timeGrowth)  ) + '\t'
            + std::to_string( pProperties->getRMS() )  + '\t'
            + std::to_string( pProperties->getMicroroughness() )  + '\t';

    for ( auto &p:m_processMap)
        output += std::to_string( p.first->getNumEventHappened() ) + '\t';

    for ( auto &p:m_processMap)
        output += std::to_string( p.second.size() ) + '\t';

    if ( m_bReportCoverages ) {
        unordered_map<string, double> covs = pLattice->computeCoverages( pParameters->getCoverageSpecies() );

        for ( auto &p:covs)
            output += std::to_string( p.second ) + '\t';
    }

    pIO->writeInOutput( output );

    if ( m_bHasGrowth )
        pIO->writeLatticeHeights( m_dProcTime );

    if ( m_bReportCoverages )
        pIO->writeLatticeSpecies( m_dProcTime  );
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

string Apothesis::mf_analyzeProc(string process){

    vector<string> parts = pIO->split( process, "->");

    //It is reaction or adsoprtion
    if ( pIO->contains( parts[0], "+" ) ){

        vector<string> reactants = pIO->split( parts[ 0 ], "+" );
        for (string s:reactants){
            pIO->trim(s);
            if (  pIO->analyzeCompound( s ).first.compare("*") == 0 )
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
