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

    /// Sets the specific adsorption species label according to the input
    void setAdrorbed(string adsorbed){ m_sAdsorbed = adsorbed;}

    /// Set the number of sites that this adsorbed occupies.
    inline void setNumSites( int i ) { m_iNumSites = i;}

    /// Get the number of sites that this adsorbed occupies.
    inline int getNumSites() { return m_iNumSites;}

    /// Get the sticking coefficient value
    inline double getStickingCoef() { return m_dStick; }

    /// Get the molar fraction [-]
    inline double getMolarFraction() { return m_dF; }

    /// Get the concentration of sites [sites/m^2]
    inline double getSitesConc() { return m_dCtot; }

    /// Get the concentration of sites [sites/m^2]
    inline double getMolecularWeight() { return m_dMW; }

    /// Get the rate coefficient of this process. It must calculated by a "Type".
    inline void setRateCoefficient(double val ) { m_dRateConstant = val; }

    /// Get the adsorption rate given as input from the user with the constant keyword
    inline double getAdsorptionRate() { return m_dAdsorptionRate; }

    /// Calculates the neighbors of a given site - To be transerred to process?
    int calculateNeighbors(Site*);

    /// Counts the vacants sites - To be transerred to process?
    int countVacantSites( Site* s);

    /// Returns the label of the adsorbed species
    inline string getAdsorbedSpecies() { return m_sAdsorbed; }

private: //pointers to functions

    /// Pointers to functions in order to switch between different functionalities
    void (Adsorption::*m_fType)();
    bool (Adsorption::*m_fRules)(Site*);
    void (Adsorption::*m_fPerform)(Site*);

private: //types

    /// The simple type for the adsorption process rate i.e.
    /// simple s0*f*P/(2*pi*MW*Ctot*kb*T) -> Sticking coefficient [-], f [-], C_tot [sites/m2], MW [kg/mol]
    void simpleType();

    /// The arrhenius type for the adsorption process rate i.e.
    /// arrhenius v0 A exp(-nE/kT), A = exp((E-Em)/kT) -> frequency v0 [-],  E (Joules), Em [Joules]
    void arrheniusType();

    /// Constant value for the adsorption process rate i.e.
    /// constant 1.0 [ML/s]
    void constantType();

private: //rules

    /// The uncoditional rule. The process is accepted without checked.
    inline bool uncoRule(Site*) { return true; }

    /// The basic rule for accepting this process.
    /// Check if the site is empty (i.e. the label is the same as the lattice species)
    /// then returns true (the processes can be performed).
    bool basicRule( Site*);

    /// For adsorbing different species in a single site must not be occupied (and TODO: the height must be the same)
    bool multiSpeciesSimpleRule( Site*);

    /// For adsorbing different species the sites must not be occupied (and TODO: the height must be the same)
    bool multiSpeciesRule( Site*);

private: //perform

    /// The process is PVD
    void signleSpeciesSimpleAdsorption( Site* );

    /// The process is PVD for multiple sites
    void signleSpeciesAdsorption( Site*);

    /// The process is CVD or ALD
    void multiSpeciesSimpleAdsorption( Site*);

    /// The process is CVD or ALD for multiple sites
    void multiSpeciesAdsorption( Site*);

private: //data

    /// Checks if the site is in lower step (only for simple cubic lattice)
    bool isInLowerStep( Site* s );

    /// Checks if the site is in higher step (only for simple cubic lattice)
    bool isInHigherStep( Site* s );

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
