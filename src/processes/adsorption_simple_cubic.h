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
#ifndef ADSORPTIONSIMPLECUBIC_H
#define ADSORPTIONSIMPLECUBIC_H

#include "process.h"

namespace MicroProcesses
{

class AdsorptionSimpleCubic: public Process
{
public:
    AdsorptionSimpleCubic();
    ~AdsorptionSimpleCubic() override;

    double getProbability() override;
    bool rules( Site* ) override;
    void perform( Site* ) override;

    inline void setActivationEnergy( double nrg ){ m_dActNrg = nrg; }
    inline double getActivationEnergy(){ return m_dActNrg; }

    inline void setMolFrac( double val ){ m_dMolFrac = val; }
    inline double getMolFrac(){ return m_dMolFrac; }

    inline void setTargetSite( Site* site ){ m_Site = site;}
    inline Site* getTargetSite(){ return m_Site; }

    inline void setSpecies( species_new* s ){ m_Species = s; }
    inline species_new* getSpecies(){ return m_Species; }

private:
    ///The activation energy of the adsoprtion process
    double m_dActNrg;

    ///The mole fraction of the AdsorptionSimpleCubic process
    double m_dMolFrac;

    ///The site that AdsorptionSimpleCubic will be performed
    Site* m_Site;

    ///The species that must adsopt
    species_new* m_Species;

    /// A member function to calculate the neighbors of a given site
    int mf_calculateNeighbors(Site*);

    REGISTER_PROCESS(AdsorptionSimpleCubic)
};
}

#endif // AdsorptionSimpleCubic_H
