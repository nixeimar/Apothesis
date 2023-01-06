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
    void init(vector<string> params) override;

    inline void setReactants( vector<pair<string, double> > reactants ) {m_vReactants = reactants;}
    inline void setProducts( vector<pair<string, double> > products ) {m_vProducts = products;}

private:
    /// The activation energy of this reaction
    double m_dEa;

    /// The pre-exponential factors for this reaction
    double m_dK0;

    /// The reactants participating in this reaction
    vector< pair<string, double> > m_vReactants;

    /// The products of this reaction
    vector< pair<string, double> > m_vProducts;
};

#endif // REACTION_NEW_H
