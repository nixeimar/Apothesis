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

namespace MicroProcesses
{

class Desorption:public Process
{
public:
    Desorption();
    ~Desorption() override;

    inline void setTargetSite( Site* site ){ m_Site = site;}
    inline Site* getTargetSite(){ return m_Site; }

    double getProbability() override;
    bool rules( Site* s) override;
    void perform( Site* ) override;
    void init(vector<string> params) override;

    void arrhenius( double, double, double, int);

    /// Sets the specific adsorption species label according to the input
    void setDesorbed(string desorbed){ m_sDesorbed = desorbed;}

    /// If keyrowd "all" is added then this is true
    inline void setAllNeighs( bool all ){  m_bAllNeihs = all; }

private:

    /// Pointers to functions in order to switch between different functions
    void (Desorption::*m_fType)();
    bool (Desorption::*m_fRules)(Site*);
    void (Desorption::*m_fPerform)(Site*);

    void mf_constantType();

    bool mf_isInLowerStep( Site* s );
    bool mf_isInHigherStep( Site* s );

    /// If the keyword 'all' is used then the rule is based on the neighbours
    bool mf_allRule(Site* s);

    /// Returns always true - this is actually as having uncoditional acceptance
    bool mf_basicRule(Site* s);

    /// For desorbing different species the site must be occupied
    bool mf_difSpeciesRule(Site* s);

    /// The process is PVD
    void mf_singleSpeciesSimpleDesorption(Site*);

    /// The process is CVD or ALD
    void mf_multiSpeciesSimpleDesorption(Site*);

    ///The site that adsorption will be performed
    Site* m_Site;

    /// A member function to calculate the neighbors of a given site
    int mf_calculateNeighbors(Site*);

    /// The number of neighbours of this process
    int m_iNumNeighs;

    bool mf_isPartOfGrowth();

    /// The species to be asdorbed
    string m_sDesorbed;

    /// If the user has "all" keyword this is set to true
    bool m_bAllNeihs;

    /// The desorption rate given as input from the user with the constant keyword
    double m_dDesorptionRate;

    REGISTER_PROCESS(Desorption)
};
}

#endif // DesorptionSimpleCubicSIMPLECUBIC_H
