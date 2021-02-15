#ifndef DESORPTION_NEW_H
#define DESORPTION_NEW_H

#include "process_new.h"

class Desorption_new:public Process_new
{
public:
    Desorption_new();
    ~Desorption_new();

    inline void setActivationEnergy( double nrg ){ m_dActNrg = nrg; }
    inline double getActivationEnergy(){ return m_dActNrg; }

    inline void setTargetSite( Site* site ){ m_Site = site;}
    inline Site* getTargetSite(){ return m_Site; }

    inline void setSpecies( species_new* s ){ m_Species = s; }
    inline species_new* getSpecies(){ return m_Species; }

    inline void setNumNeigh( int n ){ m_iNeigh = n; }

    double getProbability() override;

    void rules(set<string, std::any>) override {}
    void perform( int siteID  ) override;

private:
    ///The activation energy of the adsoprtion process
    double m_dActNrg;

    ///The site that adsorption will be performed
    Site* m_Site;

    ///The species that must be removed from the site
    species_new* m_Species;

    ///The neighbours of this diffusion process
    int m_iNeigh;
};

#endif // DESORPTION_NEW_H
