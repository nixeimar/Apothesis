#include "reaction_new.h"
#include "species_new.cpp"

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
