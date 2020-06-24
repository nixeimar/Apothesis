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

#ifndef SURFACE_REACTION_H
#define SURFACE_REACTION_H

#include "process.h"
#include "apothesis.h"
#include "site.h"

using namespace std; 
using namespace SurfaceTiles;

/* The surface reaction class. */
namespace MicroProcesses{

class SurfaceReaction: public Process
{
	public:
		/// Constructor
		SurfaceReaction
		(
			vector<string> species,
			vector<double> stoichiometry,
			double energy,
			double preExpFactor
		);
		
		/// Destructor
		virtual ~SurfaceReaction();
		
		void init();

		/// The name of the process.
		void setName(string s) { m_sName = s;}

		/// Compute the overall probabilities of this process and return it.
		double getProbability();

    	/// Returns the name of the process.
    	string getName();

    	/// Set the instance of Apothesis.
    	/// This allows to have access to all other functionalities of the KMC class.
    	void setInstance( Apothesis* apothesis ){ m_apothesis = apothesis; }

    	/// Set the lattice uppon which diffusion will be performed.
    	void activeSites( Lattice* );

    	/// The initial site that the diffusion will begin.
    	void selectSite();

    	/// Perform the process. This is what actually is called by the main KMC instance.
    	void perform();

    	/// A process map which is used between the different processes.
    	/// This need to be evaluated since it is not so handy.
    	void setProcessMap( map< Process*, list<Site* >* >* );

    	/// The list of active sites for diffusion.
    	list<Site* > getActiveList();

    	/// Here various tests should be putted in order to check for the validity of the process e.g.
    	/// the number of the particles in the active surface must be constant (mass is constant).
    	void test();

    	/// Returns true if the process can be performed in the site that callls it.
    	bool controltRules( Site* site );

	protected:

	    /// The kmc instance.
	    Apothesis* m_apothesis;

	    /// The name of the process
	    string m_sName;

	    /// Remove a site from a list
	    void mf_removeFromList();

	    /// Add a site to a list
	    void mf_addToList(Site* s);

	    /// Update the neighbour sites of this site. This is performed here since depending on the process
	    /// this changes. E.g. fotr adsortpion in a FCC lattice the first neighbors are different compared to the
	    /// adsroption in a BCC lattice.
	    void mf_updateNeighNum();

	    /// The site that diffusion is performed
	    Site* m_site;

	    /// The value of the probability of the process is stored here
	    double m_dProbability;

	    /// The number of neighs of this site
	    int m_iNeighNum;

	     /// Species that can diffuse
	    vector<string> m_reactionSpecies;

	    /// Stoichiometric coefficients of the reaction species.
		/// Reactants are negative, products are positive.
	    vector<double> m_stoichiometry;

	    /// Energy of the reaction
	    double m_energy;

		/// Pre-exponential factor of the reaction
		double m_preExpFactor;

	private:

	    /** The lattice of the process */
	    Lattice* m_pLattice;

	    /// The diffusion list which hold all the available sites for deposition
	    list<Site* > m_lAdsSites;

	    /** Pointer to the process map */
	    map< Process*, list<Site*>* >* m_pProcessMap;
};


}

#endif /// SURFACE_REACTION_H
