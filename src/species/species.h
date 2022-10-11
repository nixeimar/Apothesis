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

#ifndef SPECIES_H
#define SPECIES_H

#include <string>

using namespace std;

/** The basic class of the kinetic monte carlo code. */

class Species
{
public:
    Species
    (
        string name,
        double mw,
        double stoicCoeff,
        int id
    );

    Species
    (
        string name,
        double mw,
        int id
    );
    
    virtual ~Species();

    string getName();
    
    double getMW();

    double getStoicCoeff();

    const int getId();

protected:

    // Name of species
    string m_name;

    // Molecular weight
    double m_mw;

    // Stoichiometric coefficient
    double m_stoicCoeff;

    // Species id
    int m_id;
};

#endif // SPECIES_H
