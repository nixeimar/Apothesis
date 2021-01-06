#ifndef ADSORPTION_NEW_H
#define ADSORPTION_NEW_H

#include <iostream>
#include <string>

#include "process_new.h"
#include "site.h"
#include "species_new.h"

using namespace std;

class adsorption_new: public process_new
{
public:
    adsorption_new();
    virtual ~adsorption_new();

    inline void setActivationEnergy( double nrg ){ m_dActNrg = nrg; }
    inline double getActivationEnergy(){ return m_dActNrg; }

    inline void setMolFrac( double val ){ m_dMolFrac = val; }
    inline double getMolFrac(){ return m_dMolFrac; }

    inline void setTargetSite( Site* site ){ m_Site = site;}
    inline Site* getTargetSite(){ return m_Site; }

    inline void setSpecies( species_new* s ){ m_Species = s; }
    inline species_new* getSpecies(){ return m_Species; }

    void perform() override;

private:
    ///The activation energy of the adsoprtion process
    double m_dActNrg;

    ///The mole fraction of the adsorption process
    double m_dMolFrac;

    ///The site that adsorption will be performed
    Site* m_Site;

    ///The species that this
    species_new* m_Species;
};

#endif // ADSORPTION_NEW_H
