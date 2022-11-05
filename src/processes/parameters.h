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

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "pointers.h"
#include "apothesis.h"
#include <iostream>
#include <any>

using namespace std;

namespace Utils {

/** A class which hold all the parameters needed by KMC.
 * Other parameters needed by the individual processes can be defined there. */


class Parameters: public Pointers
  {
  public:
    /// Constructor.
    Parameters( Apothesis* apothesis );

    /// Destructor.
    ~Parameters(){;}

    /// Set the temperature value.
    inline void setTemperature( double T) { m_dT = T; }

    /// Get the temperature value.
    inline double getTemperature() { return m_dT; }

    /// Set the pressure value.
    inline void setPressure( double P) { m_dP = P; };

    /// Get the pressure value.
    inline double getPressure() { return m_dP; }

    /// Set the total number of KMC iterations
   // inline void setIterations( int iter ) { m_iIter = iter; }

    /// Set the total time for KMC
    inline void setEndTime( double time ) { m_dTime = time; }

    /// Get the total time
    inline double getEndTime() { return m_dTime; }

    /// The Avogadro number.
    const double dAvogadroNum = 6.022141793e+23;

    /// The boltzmann constant.
    const double dkBoltz = 1.3806503e-23;

    /// Pi.
    const double dPi = 3.14159265;

    /// R value (J/mol K)
    const double dUniversalGasConst = 8.3145;

    /// Store the processes to be created by the factory method.
    void setProcess(string, vector< any > );

    /// Store the initial value for the random generator
    inline void setRandGenInit( double val ) { m_dRand = val; }

    /// Store the initial value for the random generator
    inline double getRandGenInit(){return m_dRand; }

    /// Get the processes to be created.
    map< string,  vector< any > > getProcessesInfo() { return m_mProcs; }

    /// Set when to write the log
    inline void setWriteLogTimeStep( double val ) { m_dWriteLogEvery = val; }

    /// Get the time step to write lattice file
    inline double getWriteLogTimeStep() { return m_dWriteLogEvery; }

    /// Set when to write the lattice
    inline void setWriteLatticeTimeStep( double val ) { m_dWriteLatticeEvery = val; }

    /// Set when to write lattice file
    inline double getWriteLatticeTimeStep() { return m_dWriteLatticeEvery; }

    /// Print parameters info
    void printInfo();

  protected:
    /// The temperature.
    double m_dT;

    /// The pressure.
    double m_dP;

    /// The time to run kmc.
    double m_dTime;

    /// The random generator initializer
    double m_dRand;

    /// Stores the processes as read from the input file allong with their parameters.
    map< string,  vector< any > > m_mProcs;

    /// The time step to write to log
    double m_dWriteLogEvery;

    /// The time step to write the lattice
    double m_dWriteLatticeEvery;
  };

}

#endif // PARAMETERS_H
