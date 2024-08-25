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
#include "site.h"
#include <iostream>
#include <any>

using namespace std;
using namespace SurfaceTiles;

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

    /// Set the total time for KMC
    inline void setEndTime( double time ) { m_dTime = time; }

    /// Get the total time
    inline double getEndTime() { return m_dTime; }

    /// The Avogadro number [1/mol]
    const double dAvogadroNum = 6.022141793e+23;

    /// The boltzmann constant in [J/K]
    const double dkBoltz = 1.3806503e-23;

    /// Pi
    const double dPi = 3.14159265;

    /// R value (J/mol K)
    const double dUniversalGasConst = 8.3145;

    /// Store the processes to be created by the factory method.
    void setProcess(string, vector< string > );

    /// Store the initial value for the random generator
    inline void setRandGenInit( int val ) { m_iRand = val; }

    /// Store the initial value for the random generator
    inline int getRandGenInit(){return m_iRand; }

    /// Get the processes to be created.
    map< string,  vector< string > > getProcessesInfo() { return m_mProcs; }

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

    /// The label for the species in the lattice
    inline void setLatticeLabels( string label ){ m_sLatticeLabel = label; }

    /// The label for the species in the lattice
    inline string getLatticeLabels(){ return m_sLatticeLabel; }

    /// Instert a species that participates in the growth of the film
    inline void insertInGrowthSpecies( string s ){ m_vsGrowthSpecies.push_back( s ); }

    /// Returns the species that participates in the growth of the film
    inline vector<string> getGrowthSpecies(){ return m_vsGrowthSpecies; }

    /// The species to compute coverage for
    inline void setCoverageSpecies( vector<string> species ){ m_vCovSpecies = species; }

    /// Returns the species for which to compute coverage for
    inline vector<string> getCoverageSpecies(){ return m_vCovSpecies; }

    /// Sets the type of the lattice
    inline void setLatticeType( string type ){ m_sLatticeType = type; }

    /// Returns the type of the lattice
    inline string getLatticeType(){ return m_sLatticeType; }

    /// set the total number of runs
    inline void setRuns(int runs){m_sRuns = runs;}

    /// returns the number of runs remaining
    inline int getRuns(){return m_sRuns;}

    /// Sets the dimensions of the lattice
    inline void setLatticeXDim( int x ){ m_iX = x; }
    inline void setLatticeYDim( int y ){ m_iY = y; }
    inline void setLatticeHeight( int h ){ m_iH = h; }
    inline void setHeightFileExists(bool var){m_heightFileExist = var;}

    inline void setHeightData(vector<vector<int>> height){

        m_height.resize( height.size() );
        for ( int i = 0; i < height.size(); i++)
            m_height[ i ].resize( height[ i ].size() );

        for ( int i = 0; i < height.size(); i++)
            for ( int j = 0; j < height[ i ].size(); j++)
                m_height[ i ][ j ] = height[ i ][ j ];
    }

    inline int getLatticeXDim(){ return m_iX; }
    inline int getLatticeYDim(){ return m_iY; }
    inline int getLatticeHeight(){ return m_iH; }

    inline bool getHeightFileExist(){ return m_heightFileExist; }
    inline vector<vector<int>> getHeightData(){return m_height;}
   
  protected:

    /// Parameters of the lattice
    int m_iX, m_iY, m_iH;

    /// The type of the lattice
    string m_sLatticeType;

    /// The temperature [K].
    double m_dT;

    /// The pressure [Pascal].
    double m_dP;

    /// The time to run kmc [s].
    double m_dTime;

    /// The random generator initializer
    double m_iRand;

    /// Stores the processes as read from the input file allong with their parameters.
    map< string,  vector< string > > m_mProcs;

    /// The time step to write to log
    double m_dWriteLogEvery;

    /// The time step to write the lattice
    double m_dWriteLatticeEvery;

    /// The label of the lattice species
    string m_sLatticeLabel;

    /// The species participating in the growth of the surface
    vector<string> m_vsGrowthSpecies;

    /// Map of the reactions which holds the reactants enumerated
    map<string, int> m_mReactants;

    /// The species to compute coverage for
    vector<string> m_vCovSpecies;

    /// Height Data of all sites 
    vector<vector<int>> m_height;

    /// Bool to notify if height file exists or not
    bool m_heightFileExist;

    ///total number of runs
    int m_sRuns = 0;
  };

}

#endif // PARAMETERS_H
