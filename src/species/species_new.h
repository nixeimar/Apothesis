#ifndef SPECIES_NEW_H
#define SPECIES_NEW_H

#include <iostream>
#include <string>
#include <vector>
#include <string>

using namespace std;

class reaction;

class species_new
{
public:
    species_new();
    ~species_new();

    inline void setChemFormula( string chemForm ){ m_sChemForm = chemForm; }
    inline string getChemFormula(){ return m_sChemForm; }

    inline void setID( int id ){ m_iID = id; }
    inline int getID(){ return m_iID; }

    void addReaction( reaction* reaction );
    vector< reaction* > getReactions();

    inline void setMaxReacCoreff( double maxCoeff ){ m_dMaxCoeff = maxCoeff; }
    inline double getMaxReacCoreff(){ return m_dMaxCoeff; }

private:
    /// The chemical formula of this species
    string m_sChemForm;

    /// The id of this species
    int m_iID;

    /// The reactions that species participates in
    vector< reaction* > m_vReactions;

    /// This is the maximun number of the coefficient
    /// that this species participates in a reaction
    double m_dMaxCoeff;

};

#endif // SPECIES_NEW_H
