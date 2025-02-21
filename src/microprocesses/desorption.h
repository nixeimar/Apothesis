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
#ifndef DESORPTIONSIMPLECUBIC_H
#define DESORPTIONSIMPLECUBIC_H

#include "process.h"

#include "desorption_types.h"
#include "desorption_rules.h"
#include "desorption_perform.h"

namespace MicroProcesses
{

class Desorption:public Process
{
public:
    Desorption();
    ~Desorption() override;

    bool rules( Site* s) override;
    void perform( Site* ) override;
    void init(vector<string> params) override;

    /// Sets the specific adsorption species label according to the input
    void setDesorbed(string desorbed){ m_sDesorbed = desorbed;}

    /// If keyrowd "all" is added then this is true
    inline void setAllNeighs( bool all ){  m_bAllNeihs = all; }

    /// A member function to calculate the neighbors of a given site
    int calculateNeighbors(Site*);

    /// The rate of desorption if contant type
    double getDesorptionRate() { return m_dDesorptionRate;}

    /// Returns the activation energy if Arrhenius type
    double getActivationEnergy() {return m_dEd; }

    /// Returns the vibrational frequency if Arrhenius type
    double getVibrationalFrequency() {return m_dv0; }

    /// The number of neighbors defined for this process
    double getNumNeighs(){ return m_iNumNeighs; }


protected: //pointers to functions

    /// Pointers to functions in order to switch between different functionalities
    double (*m_fType)(Desorption*);
    bool (*m_fRules)(Desorption*, Site*);
    void (*m_fPerform)(Desorption*, Site*);

private: //the data

    ///The site that adsorption will be performed
    Site* m_Site;

    /// The number of neighbours of this process
    int m_iNumNeighs;

    /// The species to be asdorbed
    string m_sDesorbed;

    /// If the user has "all" keyword this is set to true
    bool m_bAllNeihs;

    /// The desorption rate given as input from the user with the constant keyword
    double m_dDesorptionRate;

    /// The vibrational frequency  (if arrhenius)
    double m_dv0;

    /// The activation energy of the process (if arrhenius)
    double m_dEd;

    /// Checks if the site is in lower step (only for simple cubic lattice)
    bool isInLowerStep( Site* s );

    /// Checks if the site is in higher step (only for simple cubic lattice)
    bool isInHigherStep( Site* s );

    REGISTER_PROCESS(Desorption)
};
}

#endif // DesorptionSimpleCubicSIMPLECUBIC_H
