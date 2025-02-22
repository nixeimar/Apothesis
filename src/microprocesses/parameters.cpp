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


Parameters& Parameters::operator=(const Parameters& other) {

    // Copy primitive data members
    m_dT = other.m_dT;
    m_dP = other.m_dP;
    m_dTime = other.m_dTime;
    m_iRand = other.m_iRand;
    m_dWriteLogEvery = other.m_dWriteLogEvery;
    m_dWriteLatticeEvery = other.m_dWriteLatticeEvery;
    m_iX = other.m_iX;
    m_iY = other.m_iY;
    m_iH = other.m_iH;
    m_dStartTime = other.m_dStartTime;
    m_bReadHeightsFromFile = other.m_bReadHeightsFromFile;
    m_bReadSpeciesFromFile = other.m_bReadSpeciesFromFile;
    m_bHasSteps = other.m_bHasSteps;
    m_iSteps = other.m_iSteps;
    m_iHeightStep = other.m_iHeightStep;
    m_iNumCycles = other.m_iNumCycles;

    // Copy string members
    m_sProcess = other.m_sProcess;
    m_sLatticeType = other.m_sLatticeType;
    m_sLatticeLabel = other.m_sLatticeLabel;

    // Copy pair
    m_pStopCov = other.m_pStopCov;

    // Deep copy vectors
    m_vsGrowthSpecies = other.m_vsGrowthSpecies;
    m_vsEtchedSpecies = other.m_vsEtchedSpecies;
    m_vCovSpecies = other.m_vCovSpecies;

    // Deep copy maps
    m_mProcs = other.m_mProcs;
    m_mReactants = other.m_mReactants;

    return *this;
}

void Parameters::setMircoProcess( string processName, vector< string > processParams )
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

string Parameters::sProcess() const
{
    return m_sProcess;
}

void Parameters::setSProcess(const string &newSProcess)
{
    m_sProcess = newSProcess;
}

double Parameters::getDurationTime() const
{
    return m_dDurationTime;
}

void Parameters::setDurationTime(double newDDurationTime)
{
    m_dDurationTime = newDDurationTime;
}

}
