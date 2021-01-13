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

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "pointers.h"
#include "apothesis.h"
#include <iostream>

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
    inline void setTime( double time ) 
    {
      m_time = time; 
    }

    /// Set the total number of KMC iterations
    inline double getTime() 
    {
      return m_time; 
    }

    /// The Avogadro number.
    const double dAvogadroNum = 6.022141793e+23;

    /// The boltzmann constant.
    const double dkBoltz = 1.3806503e-23;

    /// Pi.
    const double dPi = 3.14159265;

    /// R value (J/molK)
    const double dR = 8.3145;

    /// Store the processes to be created by the factory method.
    void setProcess(string, vector<double> );

    /// Get the processes to be created.
    map< string,  vector< double> > getProcesses() { return m_mProcs; }

  protected:
    /// The temperature.
    double m_dT;

    /// The pressure.
    double m_dP;

    /// The number of iterations to be performed.
    double m_time;

    /// Stores the processes as read from the input file allong with their parameters.
    map< string,  vector< double> > m_mProcs;
  };

}

#endif // PARAMETERS_H
