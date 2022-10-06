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

#include "reaction_new.h"
//#include "species_new.cpp"

reaction_new::reaction_new(){}
reaction_new::~reaction_new(){}


void reaction_new::addReactants( const double coeff, species_new* species ){
    pair<double, species_new*> data;
    data.first = coeff;
    data.second = species;
    m_vpReactants.push_back( data );
}

void reaction_new::addProducts( const double coeff, species_new* species ){
    pair<double, species_new*> data;
    data.first = coeff;
    data.second = species;
    m_vpProducts.push_back( data );
}


void reaction_new::print()
{
    int iCount = 0;
    for ( pair<double, species_new*> &p:m_vpReactants){
        if ( iCount != m_vpReactants.size() - 1)
            if ( p.first != 1 )
                cout<< p.first << " " << p.second->getChemFormula() << " + " ;
            else
                cout<< p.second->getChemFormula() << " + ";
        else
            if ( p.first != 1 )
                cout<< p.first << " " << p.second->getChemFormula();
            else
                cout<< p.second->getChemFormula();

        iCount++;
    }

    cout<< " = ";

    iCount = 0;
    for ( pair<double, species_new*> &p:m_vpProducts  ){
        if ( iCount != m_vpProducts.size() - 1)
            if ( p.first != 1 )
                cout<< p.first << " " << p.second->getChemFormula() << " + " ;
            else
                cout<< p.second->getChemFormula() << " + ";
        else
            if ( p.first != 1 )
                cout<< p.first << " " << p.second->getChemFormula();
            else
                cout<< p.second->getChemFormula();

        iCount++;
    }

    cout<< endl;

}
