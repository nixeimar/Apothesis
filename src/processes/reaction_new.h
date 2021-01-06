#ifndef REACTION_NEW_H
#define REACTION_NEW_H

#include <iostream>
#include <string>
#include <vector>

#include "process_new.h"

using namespace std;

class species_new;

class reaction_new: public process_new
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
