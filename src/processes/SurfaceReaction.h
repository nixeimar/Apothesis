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
#include "adsorption.h"
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
			Apothesis* apothesis,
			vector<Species*> species,
			vector<double> const& stoichiometry,
			double energy,
			double preExpFactor,
			bool immobilized
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
		
		/// Returns copy of m_stoichiometry
		const vector<double> getStoichiometry();

    	/// Set the instance of Apothesis.
    	/// This allows to have access to all other functionalities of the KMC class.
    	void setInstance( Apothesis* apothesis ){ m_apothesis = apothesis; }

    	/// Set the lattice uppon which diffusion will be performed.
    	void activeSites( Lattice* );

    	/// The initial site that the diffusion will begin.
    	void selectSite();

    	/// Perform the process. This is what actually is called by the main KMC instance.
    	void perform();

	 	/// The list of active sites for diffusion.
    	list<Site*> getActiveList();

    	/// Here various tests should be putted in order to check for the validity of the process e.g.
    	/// the number of the particles in the active surface must be constant (mass is constant).
    	void test();

    	/// Returns true if the process can be performed in the site that callls it.
    	bool controltRules( Site* site );

    	/// Compares the species on a site and checks to see if it can react
    	bool canReact(Site* site);

		/// To be deleted
		void setProcessMap(map< Process*, list<Site* >* >* procMap );

	protected:

	    /// The kmc instance.
	    Apothesis* m_apothesis;

	    /// The name of the process
	    string m_sName;

	    /// Remove a site from a list
	    void mf_removeFromList();

		/// Overloaded: remove specific site from list
		void mf_removeFromList(Site* s);

	    /// Add a site to a list
	    void mf_addToList(Site* s);

	    /// The value of the probability of the process is stored here
	    double m_dProbability;

	    /// The number of neighs of this site
	    int m_iNeighNum;

	    /// Reactants
	    vector<Species*> m_reactants;

		/// Products
		vector<Species*> m_products;

	    /// Stoichiometric coefficients of the reaction species.
		/// Reactants are negative, products are positive.
	    const vector<double> m_stoichiometry;

		/// Redundant but useful: reactant and product stoichiometry. All positive.
		vector<double> m_stoichReactants;

		vector<double> m_stoichProducts;

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

		/// can you further adsorb after required elements
		bool m_immobilized;

		/// Active sites
		int m_activeSites;
};


}

#endif /// SURFACE_REACTION_H
