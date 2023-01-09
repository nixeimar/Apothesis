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
#include <deque>

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

    inline void setReactants( map<string, int> reactants ) {m_mReactants = reactants;}
    inline void setProducts(  map<string, int> products ) {m_mProducts = products;}

private:
    /// Pointers to functions in order to switch between different functions
    void (Reaction::*m_fType)();
    bool (Reaction::*m_fRules)(Site*);
    void (Reaction::*m_fPerform)(Site*);

    /// Arrhenius type rate
    void arrheniusType( double, double, double);

    /// Constant rate
    void constantType();

    /// Reactions without growth taken into account
    void catalysis(Site* s);

    /// The reactants participating in this reaction
    map<string, int> m_mReactants;

    /// The products of this reaction
    map<string, int> m_mProducts;

    /// Checks if s is a reactant
    bool isReactant(Site* s);

    double m_dReactionRate;

    vector< set<int> > m_idReacting;

};

#endif // REACTION_NEW_H
