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
    double getRateConstant() override;
    void init(vector<string> params) override;

    inline void setReactants( unordered_map<string, int> reactants ) {m_mReactants = reactants;}
    inline void setProducts(  unordered_map<string, int> products ) {m_mProducts = products;}

    inline void setReactants( vector<string> reactants ) { m_vReactants = reactants;}
    inline void setProducts( vector<string> products ) { m_vProducts = products;}

    inline void setCoefReactants( vector<int> coefReactants ) { m_vCoefReactants = coefReactants;}
    inline void setCoefProducts( vector<int> coefProducts ) { m_vCoefProducts = coefProducts;}

private:
    /// Pointers to functions in order to switch between different functions
    void (Reaction::*m_fType)();
    bool (Reaction::*m_fRules)(Site*);
    void (Reaction::*m_fPerform)(Site*);

    vector<string> m_vReactants;
    vector<int> m_vCoefReactants;

    vector<string> m_vProducts;
    vector<int> m_vCoefProducts;

    /// Arrhenius type rate
    void arrheniusType( double, double, double);

    /// Constant rate
    void constantType();

    /// Reactions without growth taken into account
    void catalysis(Site* s);

    /// 1-1 Reaction, e.g. A* + B* -> AB*
    void oneOneReaction(Site* s);
    bool oneOneRule(Site* s);

    bool simpleRule(Site* s);

    /// The reactants participating in this reaction
    unordered_map<string, int> m_mReactants;

    /// The products of this reaction
    unordered_map<string, int> m_mProducts;

    /// Checks if the site is a reactant
    bool isReactant(Site* s);

    double m_dReactionRate;

    /// If true it leads to growth.
    bool m_bLeadsToGrowth;

    /// Checks if all reactants stoichiometric coefficients are one.
    bool allReactCoeffOne();

    /// Holds the species what to be tranformed according to the reaction e.g. A + B -> C + D, then A will be replaced by C and b by D and so on and so forth.
    unordered_map<string, string> m_mTransformationMatrix;

    /// Constructs the transformation matrix.
    void buildTransformationMatrix();

    bool leadsToGrowth(Site* s);
};

#endif // REACTION_NEW_H
