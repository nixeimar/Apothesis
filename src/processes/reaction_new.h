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

#ifndef REACTION_NEW_H
#define REACTION_NEW_H

#include <iostream>
#include <string>
#include <vector>

#include "process.h"
#include "species_new.h"

using namespace std;

class reaction_new: public Process
{
public:
    reaction_new();
    ~reaction_new();

    inline vector< pair < double, species_new* > > getReactants(){ return m_vpReactants; }
    inline vector< pair < double, species_new* > > getproducts(){ return m_vpProducts; }

    void addReactants( const double coeff, species_new* species );
    void addProducts( const double coeff, species_new* species );

    inline void setActivationEnergy( double Ea) { m_dEa = Ea; }
    inline void setPreExpFactor( double k0 ){ m_dK0 = k0; }

    void print();

private:
    /// The reactants participating in this reaction
    vector< pair< double, species_new* > > m_vpReactants;

    /// The products formed by this reaction
    vector< pair< double, species_new* > > m_vpProducts;

    /// The activation energy of this reaction
    double m_dEa;

    /// The pre-exponential factors for this reaction
    double m_dK0;
};

#endif // REACTION_NEW_H
