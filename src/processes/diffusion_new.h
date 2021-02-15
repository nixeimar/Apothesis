#ifndef DIFFUSION_NEW_H
#define DIFFUSION_NEW_H

#include <process_new.h>

class Diffusion_new:public Process_new
{
public:
    Diffusion_new();
    virtual ~Diffusion_new();

    inline void setActivationEnergy( double nrg ){ m_dActNrg = nrg; }
    inline double getActivationEnergy(){ return m_dActNrg; }

    inline void setOriginSite( Site* site ){ m_originSite = site;}
    inline Site* getOriginSite(){ return m_originSite; }

    inline void setTargetSite( Site* site ){ m_targetSite = site;}
    inline Site* getTargetSite(){ return m_targetSite; }

    inline void setSpecies( species_new* s ){ m_Species = s; }
    inline species_new* getSpecies(){ return m_Species; }

    inline void setNeigh(int n ){ m_iNeighNum = n; }


    double getProbability() override;
    void rules(set<string, std::any>) override {}
    void perform( int siteID ) override;

private:
    ///The activation energy of the adsoprtion process
    double m_dActNrg;

    ///The site to for the adsorption to be removed
    Site* m_originSite;

    ///The site that adsorption will be performed
    Site* m_targetSite;

    ///The species that must be removed from the site
    species_new* m_Species;

    /// The number of neighbours for calculating the probability
    int m_iNeighNum;
};

#endif // DIFFUSION_NEW_H
