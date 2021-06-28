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
#ifndef ADSORPTION2SSIMPLE_H
#define ADSORPTION2SSIMPLE_H

#include "process.h"
#include "algorithm"

namespace MicroProcesses
{

class AdsorptionFCC1102SSimple: public Process
{
public:
    AdsorptionFCC1102SSimple();
    ~AdsorptionFCC1102SSimple() override;

    double getProbability() override;

    /// Have at least one neighbour with below layer 4 neighs occupied and be less height than this site
    bool rules( Site* ) override;

    /// Simple adsorption and then 
    void perform( Site* ) override;

    inline void setActivationEnergy( double nrg ){ m_dActNrg = nrg; }
    inline double getActivationEnergy(){ return m_dActNrg; }

private:
    /// The site the this process is performed
    Site* m_psTargetSite;

    /// Its couple
    Site* m_psCoupledSite;

    /// The particles that have to be removed
    list<Site* > m_lToRemove;

    /// The number of neighbours that must have a site below it in order to be able to adsorb
    int m_iEnableNeighs;

    /// Finds the couple of
    Site* mf_findCouple();

    /// Check if it has a couple to be used else the site must be removed
    bool mf_hasCouple( Site* );

    /// Check if the sites below are complete
    bool mf_isLowLevelComplete( Site* s );

    ///The activation energy of the adsoprtion process
    double m_dActNrg;

    ///The mole fraction of the AdsorptionSimpleCubic process
    double m_dMolFrac;

    /// A member function to calculate the neighbors of a given site
    int mf_calculateNeighbors(Site*);

    /// Count the neights at the same level of a site
    void mf_setNeighsNum( Site* );

    REGISTER_PROCESS(AdsorptionFCC1102SSimple)
};
}

#endif // AdsorptionSimpleCubic_H
