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
#include "lattice.h"
#include "io.h"
#include "errorhandler.h"
#include "parameters.h"

#include "process.h"

//using namespace MicroProcesses;

Apothesis::Apothesis( int argc, char* argv[] ):pLattice( 0 ),pIO( 0 )
  {
  m_iArgc = argc;
  m_vcArgv = argv;

  pParameters = new Parameters( this);

  /* This must be constructed before the input */
  pLattice = new Lattice( this );

  // Create input instance
  pIO = new IO(this);
  pIO->init( m_iArgc, m_vcArgv );

  // Read the input
  pIO->readInputFile();

  // Build the lattice. This should always follow the read input
  pLattice->build();
  }

Apothesis::~Apothesis()
  {
  delete pIO;
  delete pLattice;

  // Delete the processes created by the factory method
  for (vector<Process*>::iterator it = m_vProcesses.begin();
         it != m_vProcesses.end(); it++) {
       delete *it;
    }
  }

void Apothesis::init()
  {

  /// Get the processes read and create them
  map< string, vector<double> > tempMap = pParameters->getProcesses();
  // Contruct the process
  int procsCounter = 0;
  map< string, vector<double> >::iterator mapIt = tempMap.begin();
  for (; mapIt != tempMap.end(); mapIt++ ){
    Process* proc = FactoryProcess::createProcess( mapIt->first );
    if ( proc )
      m_vProcesses.push_back( proc );
    else {
      pErrorHandler->error_simple_msg("Unknown process->" + mapIt->first );
      EXIT;
      }
    }

  if ( m_vProcesses.empty() ){
    pErrorHandler->error_simple_msg("No processes found.");
    EXIT;
    }

  /// The output file name will come from the user and will have the extenstion .log
  /// This would come as a parameter from the user from the args (also the input).
  /// Now both are hard copied.
  pIO->openOutputFile( "Output");

  /// First the processes that participate in the simulation
  /// that were read from the file input and the I/O functionality
    m_vProcesses[ 0]->setInstance( this );
    m_vProcesses[ 0]->activeSites( pLattice );
    m_vProcesses[ 0]->setProcessMap( &m_processMap );
}

void Apothesis::exec()
  {
  ///Perform the number of KMC steps read from the input.
  int iterations = pParameters->getIterations();
  if ( iterations == 0){
    pErrorHandler->error_simple_msg("Zero iterations found.");
    EXIT;
    }else{
      cout << "adsorption process is being started..." << endl;
    }
    
  for ( int i = 0; i< iterations; i++)
    {
    
    ///Here a process should be selected in random. We have only one for now...

    /// Print to output
    pIO->writeLogOutput( "Time step: " + to_string( i ) );


    /// Get the active sites of the process
    list<Site*> lAdsList = m_vProcesses[ 0]->getActiveList();

    /// Check if there are available sites that it can be performed
    if (lAdsList.size() == 0){
      cout << "No more adsorption site is available. Exiting..." << endl;
      pErrorHandler->error_simple_msg( "No adsorption site is available. Last time step: " + to_string( i ) );
      EXIT;
      }

    /// Select randomly a site
    m_vProcesses[ 0 ]->selectSite();
    /// Perform the process
    m_vProcesses[ 0]->perform();

    /// The frequency that the various information are written in the file
    /// must befined by the user. Fix it ...
    /// The user should also check if the messages are written on the terminal or not.
    pIO->writeLogOutput( m_vProcesses[ 0]->getName() + " " );
    pIO->writeLatticeHeights();

    if(i==0){
      cout << "Adsorption is being performed..." << endl;
    }
    }

  }


