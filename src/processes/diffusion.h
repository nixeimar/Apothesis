//============================================================================
//    Apothesis: A kinetic Monte Calro (KMC) code for deposition processes.
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

#ifndef DIFFUSIONSIMPLECUBIC_H
#define DIFFUSIONSIMPLECUBIC_H

#include <process.h>
#include <any>

#include "diffusion_types.h"
#include "diffusion_rules.h"
#include "diffusion_perform.h"


using namespace std;

namespace MicroProcesses
{

class Diffusion:public Process
{
public:
    Diffusion();
    ~Diffusion() override;

    bool rules( Site* ) override;
    void perform( Site* ) override;
    void init(vector<string> params) override;

    /// Sets the specific diffusion species label according to the input
    void setDiffused(string diffused){ m_sDiffused = diffused;}

    /// If keyrowd "all" is added then this is true
    inline void setAllNeighs( bool all ){  m_bAllNeihs = all; }

    /// A member function to calculate the neighbors of a given site
    int calculateNeighbors(Site*);

    /// Returns the diffusion rate
    inline double getDiffusionRate() { return m_dDiffusionRate; }

    /// Returns the activation energy if Arrhenius type
    double getActivationEnergy() {return m_dEd; }

    /// Returns the activation energy if Arrhenius type
    double getDifActivationEnergy() {return m_dEdm; }

    /// Returns the vibrational frequency if Arrhenius type
    double getVibrationalFrequency() {return m_dv0; }


protected:

    /// Pointers to functions in order to switch between different functions
    double (*m_fType)( Diffusion* );
    bool (*m_fRules)(Diffusion*, Site*);
    void (*m_fPerform)(Diffusion*, Site*);

private:

    bool mf_isInLowerStep( Site* s );
    bool mf_isInHigherStep( Site* s );

    /// If the keyword 'all' is used then the rule is based on the neighbours
    bool mf_allRule(Site* s);

    /// The label of the diffused species
    string m_sDiffused;

    /// True if the species that is to be diffused belongs to the growing film. Else it is false;
    bool m_isPartOfGrowth;

    /// If the user has "all" keyword this is set to true
    bool m_bAllNeihs;

    /// The diffusion rate
    double m_dDiffusionRate;

    /// The vibrational frequency  (if arrhenius)
    double m_dv0;

    /// The activation energy of the process (if arrhenius)
    double m_dEd;

    /// The diffuision activation energy of the process (if arrhenius)
    double m_dEdm;

    REGISTER_PROCESS( Diffusion )
};

}

#endif // DiffusionSimpleCubic_H
