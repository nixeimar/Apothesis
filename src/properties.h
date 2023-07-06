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

#ifndef PROPERTIES_H
#define PROPERTIES_H

#include "pointers.h"
#include "apothesis.h"

#include "lattice/lattice.h"

using namespace std;


namespace Utils {

class Properties: public Pointers
{
public:
    Properties(Apothesis* apothesis );

    double getMicroroughness();
    double getRMS();
    double eventCountingGrowthRate( int, int, double );

    inline double getRoughness(){ return m_dRoughness; }
    double classCoverage();
    double getMeanDH();

private:
    //The roughness of the surface
    double m_dRoughness;

    //The RMS roughness
    double m_dRMS;

    //The event counting growth rate
    double m_dEvGrRate;
};

}

#endif // PROPERTIES_H
