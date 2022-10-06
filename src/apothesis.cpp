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
#include "errorhandler.h"
#include "parameters.h"
#include "properties.h"
#include "process.h"
#include "species.h"
#include "string.h"
#include "reaction_new.h"
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

    // Create input instance
    pIO = new IO(this);
    pReader = new Reader(this);
    pReader->setInputPath("./input.txt");
    pReader->parseFile();

    //FOR the FCC case
    // Build the lattice. This should always follow the read input
    //    std::cout << "Building the lattice" << std::endl;
    //   pLattice->setOrientation("110");
    //  pLattice->build();

    //For building with steps surface. Works only for simple cubic
    //pLattice->buildSteps( 20, 1, 0);
    pLattice->build();

    std::cout << "Finished building the lattice" << std::endl;

    std::cout << "Building the lattice" << std::endl;
    pLattice->build();
//    pIO->writeLatticeHeights();
    std::cout << "Finished building the lattice" << std::endl;

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
    cout << "Opening output file" << endl;
    /// The output file name will come from the user and will have the extenstion .log
    /// This would come as a parameter from the user from the args (also the input).
    /// Now both are hard copied.

    if (!pIO->outputOpen())
        pIO->openOutputFile("Output");

    // Initialize Random generator.
    pRandomGen->init( 21321200 );

    //To Deifilia: This must be created for each process in order to pass
    //the parameters from the input file to the porcess
    map<string, any> params;

//        params.insert( {"T", 473.} );
//        params.insert( {"f", 0.0000593894333333333} ); //
//       params.insert( {"s0",  0.0262442285521079} ); //

    params.insert( {"T", 533.} );
    params.insert( {"f", 0.0000340636296296296} ); //Arrhenius
    params.insert( {"s0", 0.137295210293006} ); //Arrhenius

//    params.insert( {"f", 0.0000397675925925926} ); //LH
//    params.insert( {"s0", 0.0976139043395484} ); //LH

    //    params.insert( {"T", 573.} );
    //   params.insert( {"f", 0.0000399932074074074} ); //LH
    //    params.insert( {"s0", 0.0946803966934485} ); //LH

  //  params.insert( {"T", 623.} );
  //  params.insert( {"f", 0.000012289962962963} );//Arrhenius
//   params.insert( {"s0", 0.670471054574958} ); //Arrhenius

//      params.insert( {"f", 0.0000401917740740741} ); //LH
//      params.insert( {"s0", 0.0913294779849493} ); //LH

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
=======
  cout << "Opening output file" << endl;
  /// The output file name will come from the user and will have the extenstion .log
  /// This would come as a parameter from the user from the args (also the input).
  /// Now both are hard copied.

  if (!pIO->outputOpen())
    pIO->openOutputFile("Output-700K");

  // Processes in this case
  vector<string> pProc = m_processes;

  cout << pProc[0] << endl;

  Document &doc = pRead->getDoc();

  pIO->writeLogOutput("Initializing instances of species");
  if (std::find(pProc.begin(), pProc.end(), "Reaction") == pProc.end())
  {
    vector<double> mws = pRead->getMWs();
    vector<string> names = pRead->getSpeciesNames();

    for (int i = 0; i < mws.size(); ++i)
    {
      Species *s = new Species(names[i], mws[i], m_nSpecies);
      m_nSpecies++;
      m_species[names[i]] = s;
    }
  }

  // Initializing interactions between species
  pIO->writeLogOutput("Reading interactions between species");
  Value &speciesName = doc["Species"];
  for (Value::ConstMemberIterator itr = speciesName.MemberBegin(); itr != speciesName.MemberEnd(); ++itr)
  {
    const char *name = itr->name.GetString();
    Value &singleSpecies = speciesName[name];
    if (singleSpecies.HasMember("Interactions"))
    {
      Value &interactions = singleSpecies["Interactions"];

      pIO->writeLogOutput("Initializing interactions between species...");

      int numInters = interactions.Size();
      for (int i = 0; i < numInters; i++)
      {
        m_interactions.push_back(make_tuple(name, interactions[i].GetString()));
        pIO->writeLogOutput("(" + get<0>(m_interactions.back()) + ", " + get<1>(m_interactions.back()) + ")");
      }
    }
  }

  pIO->writeLogOutput("Initializing processes");

  if (std::find(pProc.begin(), pProc.end(), "Adsorption") != pProc.end())
  {
    pIO->writeLogOutput("Initializing Adsorption");

    // Read parameters for Adsorption
    Value &specie = doc["Process"]["Adsorption"]["Species"];
    Value &stick = doc["Process"]["Adsorption"]["Sticking"];
    Value &mFrac = doc["Process"]["Adsorption"]["MassFraction"];
    Value &ads = doc["Process"]["Adsorption"];

    // Verify presence of each parameter in input file
    logSuccessfulRead(specie.IsArray(), "Adsorption species");
    logSuccessfulRead(stick.IsArray(), "Adsorption sticking coefficients");
    logSuccessfulRead(mFrac.IsArray(), "Adsorption mass fraction");

    // Initialize vectors
    vector<string> species;
    vector<double> sticking;
    vector<double> massFraction;

    // Sum of mass fraction. Later used to normalize.
    double sum = 0;

    for (SizeType i = 0; i < specie.Size(); i++)
    {
      // Output possible errors
      if (!specie[i].IsString())
        pErrorHandler->error_simple_msg("Species format is not a string");
      if (!stick[i].IsNumber())
        pErrorHandler->error_simple_msg("Sticking coefficient format is not a double");
      if (!mFrac[i].IsNumber())
        pErrorHandler->error_simple_msg("Mass fraction format is not a double");

      // Push values to corresponding vectors
      species.push_back(specie[i].GetString());
      sticking.push_back(stick[i].GetDouble());
      massFraction.push_back(mFrac[i].GetDouble());
      sum += massFraction[i];
    }

    // Normalize the values of the mass fraction
    for (vector<double>::iterator itr = massFraction.begin(); itr != massFraction.end(); ++itr)
    {
      *itr = *itr / sum;
    }

    for (int i = 0; i < species.size(); ++i)
    {
      bool direct = false;
      if (ads.HasMember("Direct"))
      {
        string dir = ads["Direct"].GetString();
        if (dir.compare("Yes") == 0)
        {
          direct = true;
>>>>>>> master
        }
    } <---------------- */


        auto pos = m_processMap.insert( { FactoryProcess::createProcess("AdsroptionSimpleCubic4sMulti"), emptySet } );
        pos.first->first->setName("AdsortpionSimpleCubinc4SMulti");
        pos.first->first->init( params );
        pos.first->first->setLattice( pLattice );
        pos.first->first->setRandomGen( pRandomGen );

        auto des = m_processMap.insert( { FactoryProcess::createProcess("DesorptionSimpleCubic4sMulti"), emptySet } );
        des.first->first->setName("0Neighs");
    //    des.first->first->setConstant( true );
        des.first->first->init( params );
        des.first->first->setLattice( pLattice );
        des.first->first->setRandomGen( pRandomGen );
        des.first->first->setNeighs( 0 );

        auto des1 = m_processMap.insert( { FactoryProcess::createProcess("DesorptionSimpleCubic4sMulti"), emptySet } );
        des1.first->first->setName("1Neighs");
    //    des.first->first->setConstant( true );
        des1.first->first->init( params );
        des1.first->first->setLattice( pLattice );
        des1.first->first->setRandomGen( pRandomGen );
        des1.first->first->setNeighs( 1 );

        auto des2 = m_processMap.insert( { FactoryProcess::createProcess("DesorptionSimpleCubic4sMulti"), emptySet } );
        des2.first->first->setName("2Neighs");
    //    des.first->first->setConstant( true );
        des2.first->first->init( params );
        des2.first->first->setLattice( pLattice );
        des2.first->first->setRandomGen( pRandomGen );
        des2.first->first->setNeighs( 2 );

        auto des3 = m_processMap.insert( { FactoryProcess::createProcess("DesorptionSimpleCubic4sMulti"), emptySet } );
        des3.first->first->setName("3Neighs");
    //    des.first->first->setConstant( true );
        des3.first->first->init( params );
        des3.first->first->setLattice( pLattice );
        des3.first->first->setRandomGen( pRandomGen );
        des3.first->first->setNeighs( 3 );

        auto des4 = m_processMap.insert( { FactoryProcess::createProcess("DesorptionSimpleCubic4sMulti"), emptySet } );
        des4.first->first->setName("4Neighs");
    //    des.first->first->setConstant( true );
        des4.first->first->init( params );
        des4.first->first->setLattice( pLattice );
        des4.first->first->setRandomGen( pRandomGen );
        des4.first->first->setNeighs( 4 );

        auto des5 = m_processMap.insert( { FactoryProcess::createProcess("DesorptionSimpleCubic4sMulti"), emptySet } );
        des5.first->first->setName("5Neighs");
    //    des.first->first->setConstant( true );
        des5.first->first->init( params );
        des5.first->first->setLattice( pLattice );
        des5.first->first->setRandomGen( pRandomGen );
        des5.first->first->setNeighs( 5 );

        auto des6 = m_processMap.insert( { FactoryProcess::createProcess("DesorptionSimpleCubic4sMulti"), emptySet } );
        des6.first->first->setName("6Neighs");
    //    des.first->first->setConstant( true );
        des6.first->first->init( params );
        des6.first->first->setLattice( pLattice );
        des6.first->first->setRandomGen( pRandomGen );
        des6.first->first->setNeighs( 6 );

        for ( Site* s:pLattice->getSites() )
            s->setLabel("Cu"); // in all cases we start with Cu

        for ( auto &p:m_processMap){
            cout << p.first->getName() << endl;
            for ( Site* s:pLattice->getSites() ){
                if ( p.first->rules( s ) )
                    p.second.insert( s );
            }
        }
      }
      Adsorption *a = new Adsorption(this, species[i], m_species[species[i]], sticking[i], massFraction[i], direct);

      // Keep two separate vectors: one for all processes, one for adsorption processes only
      m_vAdsorption.push_back(a);
      m_vProcesses.push_back(a);
    }

    // Resolving Conflicts

    pIO->writeLogOutput("...Done initializing Adsorption process.");
  }
  if (std::find(pProc.begin(), pProc.end(), "Desorption") != pProc.end())
  {
    pIO->writeLogOutput("Initializing Desorption");
>>>>>>> master

    /*    auto pos = m_processMap.insert( { FactoryProcess::createProcess("AdsortpionFCC1102SMulti"), emptySet } );
    pos.first->first->setName("AdsortpionFCC1102SMulti");
    pos.first->first->init( params );
    pos.first->first->setLattice( pLattice );
    pos.first->first->setRandomGen( pRandomGen );

    auto des = m_processMap.insert( { FactoryProcess::createProcess("DesorptionFCC110Multi"), emptySet } );
    string name = "Desorption HAMD"; //+ std::to_string( i );
    des.first->first->setName( name );
    des.first->first->init( params );

    for ( Site* s:pLattice->getSites() )
        s->setLabel("Cu"); // in all cases we start with Cu

    for ( auto &p:m_processMap){
        cout << p.first->getName() << endl;
        for ( Site* s:pLattice->getSites() ){
            if ( p.first->rules( s ) )
                p.second.insert( s );
        }
    } -- > Enable for FCC */

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
