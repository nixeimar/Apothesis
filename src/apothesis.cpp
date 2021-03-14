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
#include "processpool.h"

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


    // ----- Validating with the case of Lam & Vlachos ------
    set< Site* > emptySet;
    set< Site* > tempSet;
    for ( Site* s:pLattice->getSites() )
        tempSet.insert( s );

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

    cout << "HERE: " << any_cast<double>(params["T"]) << endl;
    auto pos = m_processMap.insert( { FactoryProcess::createProcess("AdsorptionSimpleCubic"), tempSet } );
    pos.first->first->setName("Adsorption");
    pos.first->first->setParams( params );

    cout << pos.first->first->getName() << " " << pos.first->first->getProbability()  << endl; // << any_cast<int>(pos.first->first->getParameter("neighs") );

    params.insert( {"v0", 1.0e+13} );
    params.insert( {"E", 1.0e+13} );

    params.insert( {"neighs", 1} );
    pos = m_processMap.insert( { FactoryProcess::createProcess("DesorptionSimpleCubic"), emptySet } );
    pos.first->first->setName("Desorption1N");
    pos.first->first->setParams( params );

    cout << pos.first->first->getName() << " " << pos.first->first->getProbability()  << " " << any_cast<int>(pos.first->first->getParameter("neighs") ) << endl;

    pos = m_processMap.insert( { FactoryProcess::createProcess("DiffusionSimpleCubic"), emptySet } );
    pos.first->first->setName("Diffusion1N");
    pos.first->first->setParams( params );

    cout << pos.first->first->getName() << " " << pos.first->first->getProbability()  << " " << any_cast<int>(pos.first->first->getParameter("neighs") ) << endl;

    params.erase("neighs");
    params.insert( {"neighs", 2} );
    pos = m_processMap.insert( { FactoryProcess::createProcess("DesorptionSimpleCubic"), emptySet } );
    pos.first->first->setName("Desorption2N");
    pos.first->first->setParams( params );


    cout << pos.first->first->getName() << " " << pos.first->first->getProbability()  << " " << any_cast<int>(pos.first->first->getParameter("neighs") ) << endl;

    pos = m_processMap.insert( { FactoryProcess::createProcess("DiffusionSimpleCubic"), emptySet } );
    pos.first->first->setName("Diffusion2N");
    pos.first->first->setParams( params );

    cout << pos.first->first->getName() << " " << pos.first->first->getProbability()  << " " << any_cast<int>(pos.first->first->getParameter("neighs") ) << endl;
    params.erase("neighs");
    params.insert( {"neighs", 3} );
    pos = m_processMap.insert( { FactoryProcess::createProcess("DesorptionSimpleCubic"), emptySet } );
    pos.first->first->setName("Desorption3N");
    pos.first->first->setParams( params );


    cout << pos.first->first->getName() << " " << pos.first->first->getProbability()  << " " << any_cast<int>(pos.first->first->getParameter("neighs") ) << endl;

    pos = m_processMap.insert( { FactoryProcess::createProcess("DiffusionSimpleCubic"), emptySet } );
    pos.first->first->setName("Diffusion3N");
    pos.first->first->setParams( params );

    cout << pos.first->first->getName() << " " << pos.first->first->getProbability()  << " " << any_cast<int>(pos.first->first->getParameter("neighs") ) << endl;

    params.erase("neighs");
    params.insert( {"neighs", 4} );
    pos = m_processMap.insert( { FactoryProcess::createProcess("DesorptionSimpleCubic"), emptySet } );
    pos.first->first->setName("Desorption4N");
    pos.first->first->setParams( params );

    cout << pos.first->first->getName() << " " << pos.first->first->getProbability()  << " " << any_cast<int>(pos.first->first->getParameter("neighs") ) << endl;

    pos = m_processMap.insert( { FactoryProcess::createProcess("DiffusionSimpleCubic"), emptySet } );
    pos.first->first->setName("Diffusion4N");
    pos.first->first->setParams( params );

    cout << pos.first->first->getName() << " " << pos.first->first->getProbability()  << " " << any_cast<int>(pos.first->first->getParameter("neighs") ) << endl;

    params.erase("neighs");
    params.insert( {"neighs", 5} );
    pos = m_processMap.insert( { FactoryProcess::createProcess("DesorptionSimpleCubic"), tempSet } );
    pos.first->first->setName("Desorption5N");
    pos.first->first->setParams( params );

    cout << pos.first->first->getName() << " " << pos.first->first->getProbability()  << " " << any_cast<int>(pos.first->first->getParameter("neighs") ) << endl;

    pos = m_processMap.insert( { FactoryProcess::createProcess("DiffusionSimpleCubic"), tempSet } );
    pos.first->first->setName("Diffusion5N");
    pos.first->first->setParams( params );

    cout << pos.first->first->getName() << " " << pos.first->first->getProbability()  << " " << any_cast<int>(pos.first->first->getParameter("neighs") ) << endl;

//    for (auto& [key, value]: m_processMap)
 //       cout << key->getName() << " " << key->getProbability() << any_cast<int>(key->getParameter("neighs") );

    //    cout << m_processMap.size() << endl;

    // Initialize Random generator.
    pRandomGen->init( 0 );
    MicroProcesses::ProrcessPool* procPool = new MicroProcesses::ProrcessPool();

    //Give initial height for lattice. This is for FCC(110).
    for ( Site* s:pLattice->getSites() )
        if ( s->getID()%2 != 0)
            s->setHeight( 9 );
        else
            s->setHeight( 10 );

    pLattice->print();
}

void Apothesis::exec()
{
    //-------------------- This is an example for FCC(110) lattice  ------------------//

    Process* adsosption = FactoryProcess::createProcess("Adsorption");
    //    Adsorption* adsosption = new Adsorption_new();
    adsosption->setName("Simple_FCC_110");
    adsosption->setID( 0 );

    pair<string, set<int> > p;
    p.first = adsosption->getName();
    set< int > ids;

    m_procMap.insert(  p );
    //    procPool->addProcess( adsosption->getName(), adsosption );
    //   procPool->addProcess( adsosption->getID(),  adsosption);

    //We always start from the even numbers in FCC 110
    for (Site* s:pLattice->getSites() )
        if ( s->getID()%2 == 0)
            m_procMap[ adsosption->getName() ].insert( s->getID() );

    // Here we set the process map to the lattice in order for the lattice to be able to modified according to the structural properties of the lattice.
    //(must be transferred to init)
    pLattice->setProcMap( &m_procMap );

    adsosption->perform( 15 );

    EXIT;

    //-------------------- This is an example for FCC(110) lattice  ------------------//


    //---------------------- Creation of the process map & initialization (must be transferred to init) - These are to reproduce Lam & Vlachos (2000) results  ------------------------------>//
    /*    Adsorption_new* adsosption = new Adsorption_new();
    adsosption->setActivationEnergy( 12.0 );
    adsosption->setName("Adsoprtion");
    adsosption->setID( 0 );

    pair<string, set<int> > p;
    p.first = adsosption->getName();
    set< int > ids;

    m_procMap.insert(    p );
    procPool->addProcess( adsosption->getName(), adsosption );
    procPool->addProcess( adsosption->getID(),  adsosption);

    for (Site* s:pLattice->getSites() )
        m_procMap[ adsosption->getName() ].insert( s->getID() );

    Desorption_new* desorption_1N = new Desorption_new();
    desorption_1N->setName("Desorption 1N");
    desorption_1N->setNumNeigh( 1 );
    desorption_1N->setID( 1 );

    p.first = desorption_1N->getName();
    m_procMap.insert( p );
    procPool->addProcess( desorption_1N->getName(), desorption_1N );
    procPool->addProcess( desorption_1N->getID(), desorption_1N );

    Desorption_new* desorption_2N = new Desorption_new();
    desorption_2N->setName("Desorption 2N");
    desorption_2N->setNumNeigh( 2 );
    desorption_2N->setID( 2 );

    p.first = desorption_2N->getName();
    m_procMap.insert( p );
    procPool->addProcess( desorption_2N->getName(), desorption_2N );
    procPool->addProcess( desorption_2N->getID(), desorption_2N );

    Desorption_new* desorption_3N = new Desorption_new();
    desorption_3N->setName("Desorption 3N");
    desorption_3N->setNumNeigh( 3 );
    desorption_3N->setID( 3 );

    p.first = desorption_3N->getName();
    m_procMap.insert( p );
    procPool->addProcess( desorption_3N->getName(), desorption_3N );
    procPool->addProcess( desorption_3N->getID(), desorption_3N );

    Desorption_new* desorption_4N = new Desorption_new();
    desorption_4N->setName("Desorption 4N");
    desorption_4N->setNumNeigh( 4 );
    desorption_4N->setID( 4 );

    p.first = desorption_4N->getName();
    m_procMap.insert( p );
    procPool->addProcess( desorption_4N->getName(), desorption_4N );
    procPool->addProcess( desorption_4N->getID(), desorption_4N );

    Desorption_new* desorption_5N = new Desorption_new();
    desorption_5N->setName("Desorption 5N");
    desorption_5N->setNumNeigh( 5 );
    desorption_5N->setID( 5 );

    p.first = desorption_5N->getName();
    m_procMap.insert( p );
    procPool->addProcess( desorption_5N->getName(), desorption_5N );
    procPool->addProcess( desorption_5N->getID(), desorption_5N );

    for (Site* s:pLattice->getSites() )
        m_procMap[ desorption_5N->getName() ].insert( s->getID() );

    Diffusion_new* diffusion_1N = new Diffusion_new();
    diffusion_1N->setName("Diffusion 1N");
    diffusion_1N->setNeigh(1);
    diffusion_1N->setID(6);

    p.first = diffusion_1N->getName();
    m_procMap.insert( p );
    procPool->addProcess( diffusion_1N->getName(), diffusion_1N );
    procPool->addProcess( diffusion_1N->getID(), diffusion_1N );

    Diffusion_new* diffusion_2N = new Diffusion_new();
    diffusion_2N->setName("Diffusion 2N");
    diffusion_2N->setID(7);

    p.first = diffusion_2N->getName();
    m_procMap.insert( p );
    procPool->addProcess( diffusion_2N->getName(), diffusion_2N );
    procPool->addProcess( diffusion_2N->getID(), diffusion_2N );

    Diffusion_new* diffusion_3N = new Diffusion_new();
    diffusion_3N->setName("Diffusion 3N");
    diffusion_3N->setNeigh(3);
    diffusion_3N->setID(8);

    p.first = diffusion_3N->getName();
    m_procMap.insert( p );
    procPool->addProcess( diffusion_3N->getName(), diffusion_3N );
    procPool->addProcess( diffusion_3N->getID(), diffusion_3N );

    Diffusion_new* diffusion_4N = new Diffusion_new();
    diffusion_4N->setName("Diffusion 4N");
    diffusion_4N->setNeigh(4);
    diffusion_4N->setID(9);

    p.first = diffusion_4N->getName();
    m_procMap.insert( p );
    procPool->addProcess( diffusion_4N->getName(), diffusion_4N );
    procPool->addProcess( diffusion_4N->getID(), diffusion_4N );

    Diffusion_new* diffusion_5N = new Diffusion_new();
    diffusion_5N->setName("Diffusion 5N");
    diffusion_5N->setNeigh(5);
    diffusion_5N->setID(10);

    p.first = diffusion_5N->getName();
    m_procMap.insert( p );
    procPool->addProcess( diffusion_5N->getName(), diffusion_5N );
    procPool->addProcess( diffusion_5N->getID(), diffusion_5N );

    for (Site* s:pLattice->getSites() )
        m_procMap[ diffusion_5N->getName() ].insert( s->getID() );

    //Set reference lattice to each process
    for (pair<string, set<int> > p:m_procMap)
        procPool->getProcessByName( p.first )->setLattice( pLattice ); */
    //<---------------------- End creation of the process map & initialization for Lam & Vlachos (2000) resutls  ------------------------------//


    ProrcessPool* procPool = new ProrcessPool();
    for (pair<string, set<int> > p:m_procMap)
        procPool->getProcessByName( p.first )->setLattice( pLattice );


    //--------------- Open files for writting ---------------------->
    pIO->openRoughnessFile( "testRough" );

    // Here we set the process map to the lattice in order for the lattice to be able to modified according to the structural properties of the lattice.
    //(must be transferred to init)
    pLattice->setProcMap( &m_procMap );

    //Perform the number of KMC steps read from the input.
    m_dEndTime = pParameters->getEndTime();

    if (m_dEndTime == 0.0){
        pErrorHandler->error_simple_msg("Zero iterations found.");
        exit(0);
    }
    else
        pIO->writeLogOutput("Running " + to_string( m_dEndTime ) + " sec");

    //Calculate the total probability (R) --------------------------//
    m_dRTot = 0.0;
    for (pair<string, set<int> > p:m_procMap)
        m_dRTot += (double)procPool->getProcessByName( p.first )->getProbability()*p.second.size();

    m_dEndTime = 0.01;

    pProperties->calculateRoughness();
    pIO->writeRoughness( m_dProcTime, pProperties->getRoughness() );

    while ( m_dProcTime <= m_dEndTime ){
        //1. Get a random numbers
        m_dRandom = pRandomGen->getDoubleRandom();

        m_dSum = 0.0;
        for (pair<string, set<int> > p:m_procMap) {
            m_dProcRate = (double)procPool->getProcessByName( p.first )->getProbability()*p.second.size();
            m_dSum += m_dProcRate/m_dRTot;

            //2. Pick a process according to the rates
            if ( m_dRandom <= m_dSum ){

                //Get a random number which is the ID of the site where this process can performed
                m_iSiteNum = pRandomGen->getIntRandom(0, m_procMap[ p.first ].size() - 1 );

                //3. From this process pick the random site with id and perform it:
                procPool->getProcessByName( p.first )->perform( *next(m_procMap[ p.first ].begin(), m_iSiteNum) );
                procPool->getProcessByName( p.first )->eventHappened();

                //4. Re-compute the processes rates and re-compute Rtot (see ppt).
                m_dRTot = 0.0;
                for (pair<string, set<int> > p:m_procMap)
                    m_dRTot += (double)procPool->getProcessByName( p.first )->getProbability()*p.second.size();

                //5. Compute dt = -ln(ksi)/Rtot
                m_dt = -log( pRandomGen->getDoubleRandom()  )/m_dRTot;
                break;
            }
        }
        //6. advance time: time += dt;
        m_dProcTime += m_dt;

        //Calulate the roughness
        pProperties->calculateRoughness();

        //Then write it
        pIO->writeRoughness( m_dProcTime, pProperties->getRoughness() );
    }

    delete procPool;
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


