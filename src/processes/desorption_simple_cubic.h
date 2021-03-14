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
#ifndef DESORPTIONSIMPLECUBIC_H
#define DESORPTIONSIMPLECUBIC_H

#include "process.h"

namespace MicroProcesses
{

class DesorptionSimpleCubic:public Process
{
public:
    DesorptionSimpleCubic();
    ~DesorptionSimpleCubic() override;

    inline void setActivationEnergy( double nrg ){ m_dActNrg = nrg; }
    inline double getActivationEnergy(){ return m_dActNrg; }

    inline void setTargetSite( Site* site ){ m_Site = site;}
    inline Site* getTargetSite(){ return m_Site; }

    inline void setSpecies( species_new* s ){ m_Species = s; }
    inline species_new* getSpecies(){ return m_Species; }

    inline void setNumNeigh( int n ){ m_iNeigh = n; }

    double getProbability() override;
    bool rules( Site* s) override;
    set<Site*> getAffectedSites() override;
    void perform( Site* ) override;

private:
    ///The activation energy of the adsoprtion process
    double m_dActNrg;

    ///The site that adsorption will be performed
    Site* m_Site;

    ///The species that must be removed from the site
    species_new* m_Species;

    ///A list holding all affected sites from this process
    set<Site*> m_seAffectedSites;

    ///The neighbours of this diffusion process
    int m_iNeigh;

    /// A member function to calculate the neighbors of a given site
    int mf_calculateNeighbors(Site*);

    REGISTER_PROCESS( DesorptionSimpleCubic )
};
}

#endif // DesorptionSimpleCubicSIMPLECUBIC_H
