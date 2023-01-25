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
#ifndef ADSORPTIONSIMPLECUBIC_H
#define ADSORPTIONSIMPLECUBIC_H

#include "process.h"

namespace MicroProcesses
{

class Adsorption: public Process
{
public:
    Adsorption();
    ~Adsorption() override;

    bool rules( Site* ) override;
    void perform( Site* ) override;
    void init( vector<string> params ) override;
    double getRateConstant() override;

    inline void setTargetSite( Site* site ){ m_Site = site;}
    inline Site* getTargetSite(){ return m_Site; }


    /// Sets the specific adsorption species label according to the input
    void setAdrorbed(string adsorbed){ m_sAdsorbed = adsorbed;}

    /// Set the number of sites that this adsorbed occupies.
    inline void setNumSites( int i ) { m_iNumSites = i;}

    /// Get the number of sites that this adsorbed occupies.
    inline int getNumSites() { return m_iNumSites;}

private:

    /// Pointers to functions in order to switch between different functions
    void (Adsorption::*m_fType)();
    bool (Adsorption::*m_fRules)(Site*);
    void (Adsorption::*m_fPerform)(Site*);

    /// The simple type for the adsorption process rate i.e.
    /// simple s0*f*P/(2*pi*MW*Ctot*kb*T) -> Sticking coefficient [-], f [-], C_tot [sites/m2], MW [kg/mol]
    void simpleType();

    /// The arrhenius type for the adsorption process rate i.e.
    /// arrhenius v0 A exp(-nE/kT), A = exp((E-Em)/kT) -> frequency v0 [-],  E (Joules), Em [Joules]
    void arrheniusType();

    /// Constant value for the adsorption process rate i.e.
    /// constant 1.0 [ML/s]
    void constantType();

    /// The process is PVD
    void signleSpeciesSimpleAdsorption(Site*);

    /// The process is PVD for multiple sites
    void signleSpeciesAdsorption(Site*);

    /// The process is CVD or ALD
    void multiSpeciesSimpleAdsorption(Site*);

    /// The process is CVD or ALD for multiple sites
    void multiSpeciesAdsorption(Site*);

    /// The uncoditional rule. The process is accepted without checked.
    bool uncoRule(Site* s);

    /// The basic rule for accepting this process.
    /// Check if the site is empty (i.e. the label is the same as the lattice species)
    /// then returns true (the processes can be performed).
    bool basicRule( Site* s);

    /// For adsorbing different species the sites must not be occupied (and TODO: the height must be the same)
    bool multiSpeciesRule(Site* s);

    /// For adsorbing different species in a single site must not be occupied (and TODO: the height must be the same)
    bool multiSpeciesSimpleRule(Site* s);

    /// Counts the vacants sites
    int countVacantSites( Site* s);

    /// Checks if the site is in lower step (only for simple cubic lattice)
    bool isInLowerStep( Site* s );

    /// Checks if the site is in higher step (only for simple cubic lattice)
    bool isInHigherStep( Site* s );

    ///The site that the process will be performed
    Site* m_Site;

    /// Calculates the neighbors of a given site
    int calculateNeighbors(Site*);

    /// For simple adsorption:
    ///The sticking coefficient [-]
    double m_dStick;

    ///The molar fraction [-]
    double m_dF;

    ///The concentration of sites [sites/m^2]
    double m_dCtot;

    ///The molecular weight of the species [kg/mol]
    double m_dMW;

    /// The species to be asdorbed
    string m_sAdsorbed;

    /// The adsorption rate given as input from the user with the constant keyword
    double m_dAdsorptionRate;

    REGISTER_PROCESS( Adsorption )
};
}

#endif // AdsorptionSimpleCubic_H
