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

#ifndef REACTION_H
#define REACTION_H

#include <iostream>
#include <string>
#include <vector>

#include "process.h"

class species;

using namespace std;

class Reaction: public Process
{
public:
    Reaction();
    ~Reaction();

    void perform(Site *) override;
    bool rules(Site *) override;
    double getProbability() override;

    inline void setActivationEnergy( double Ea) { m_dEa = Ea; }
    inline void setPreExpFactor( double k0 ){ m_dK0 = k0; }

    inline void setReactants(list<string> reactants) {m_lReactants = reactants;}
    list<string> getReactants() { return m_lReactants; }

    void setStoichiometry(string species, double stoichCoeff) { m_mStoichiometry[species] = stoichCoeff; }
    map<string, double> getStoichiometry(string species, double stoichCoeff) { m_mStoichiometry; }

    double getStoichCoeff( string species) { return m_mStoichiometry[species]; }

    void print();

private:
    /// The activation energy of this reaction
    double m_dEa;

    /// The pre-exponential factors for this reaction
    double m_dK0;

    /// The reactants participating in this reaction
    list<string> m_lReactants;

    /// The stoichiometry of the reaction e.g. map<"Cu", 2>
    map<string, double > m_mStoichiometry;

};

#endif // REACTION_NEW_H
