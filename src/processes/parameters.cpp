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

#include "parameters.h"

namespace Utils  
{

Parameters::Parameters(Apothesis* apothesis ):Pointers(apothesis), m_iRand(0), m_bReadHeightsFromFile(false),
    m_bReadSpeciesFromFile(false), m_dStartTime(0.0){}
  
  void Parameters::setProcess( string processName, vector< string > processParams )
  {
    m_mProcs[ processName ] = processParams;
  }

  void Parameters::printInfo()
  {
      cout << endl;
      cout << "--- start info simulation parameters -- " << endl;
      cout << "---------------------------------------- " << endl;
      cout << "Time "<< m_dTime << endl;
      cout << "Temperature "<< m_dT << endl;
      cout << "Pressure "<< m_dP << endl;
      cout << "Random gen init " << m_iRand << endl;
      cout << "Write in log every " << m_dWriteLogEvery << endl;
      cout << "Write lattice every " << m_dWriteLatticeEvery << endl;
      cout << "---------------------------------------- " << endl;
      cout << "--- end simulation parameters info ----- " << endl;
      cout << endl;
  }

}
