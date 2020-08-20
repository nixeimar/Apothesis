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

#ifndef DESORPTION_H
#define DESORPTION_H

#include "process.h"

#include <iostream>
#include <list>
#include <math.h>

#include "apothesis.h"
#include "adsorption.h"
#include "diffusion.h"
#include "site.h"

using namespace std;
using namespace SurfaceTiles;

namespace MicroProcesses{

/** The adsorption class. Adsoprtion depending on the type of lattice can be performed in different sites.
For the simplest case e.g. BCC lattices all the sites are available for deposition. */

class Desorption: public Process
{
  public:
    /// Constructor
    Desorption
    (
      Apothesis* instance,
      string speciesName,
      Species* species,
      double energy,
      double frequency
    );

    /// Destructor
    virtual ~Desorption();

    /// The name of the process.
    void setName(string s){ m_sName = s;}

    //To be deleted
    void setInstance( Apothesis* apothesis ){}
    
    /// Compute the overall probabilities of this processus and return it.
    double getProbability();

    /// Returns the name of the process.
    string getName();

    /// Returns name of species
    const string getSpeciesName();

    /// Returns instance of species class
    const Species* getSpecies();

    /// Constructs the sites that adsorption can be performed.
    void activeSites( Lattice* );

    /// Slect the site that the adsorption will be performed from the available sites.
    void selectSite();

    void setProcessMap( map< Process*, list<Site* >* >* );

    /// Perform the process. This is what actually is called by the main KMC instance.
    void perform();

    /// The list of active sites for adsorption.
    list<Site*> getActiveList();

    /// Here various tests should be putted in order to check for the validity of the process e.g.
    /// the number of the particles in the active surface must be constant (mass is constant).
    void test();

    /// Returns true if the process can be performed in the site that callls it.
    bool controltRules( Site* site );

    /// Set associated adsorption pointer
    void setAdsorptionPointer(Adsorption* a);

    /// Set associated diffusion pointer
    void setDiffusionPointer(Diffusion* d);

    /// Set Diffusion switch
    void setDiffusion(bool canDiffuse);

    /// Add a site to a list
    void mf_addToList(Site* s);

    /// Update counter on number of sites with n neighbours
    void updateSiteCounter(int neighbours, bool addOrRemove);

    /// Update counts on neighbours
    void updateNeighbours(Site* s);

    /// Set site
    void setSite(Site* s);
    
  protected:
    /// The kmc instance.
    Apothesis* m_apothesis;

    /// The name of the process
    string m_sName;

    /// Remove a site from a list
    void mf_removeFromList();

    /// Update the neighbour sites of this site. This is performed here since depending on the process
    /// this changes. E.g. fotr adsortpion in a FCC lattice the first neighbors are different compared to the
    /// adsroption in a BCC lattice.
    void mf_updateNeighNum();

    /// The site that desorption is performed
    Site* m_site;

    /// The value of the probability of the process is stored here
    double m_dProbability;

    /// The number of neighs of this site
    int m_iNeighNum;

    /// Return names of desorption species
    const string getDesorptionSpecies();

    /// Return array of desorption energy
    const double getDesorptionEnergy();
   
    /// Return array of desorption frequency
    const double getDesorptionFrequency();

    /// Return pointer to corresponding adsorption class
    Adsorption* getAdsorption();

    /// Return pointer to corresponding diffusion class
    Diffusion* getDiffusion();

    bool m_canDiffuse;

    bool canDiffuse();

  private:
  
    /** The lattice of the process */
    Lattice* m_pLattice;

    /// The desorption list which hold all the available sites for deposition
    list<Site* > m_lDesSites;

    /// Species name that can desorb
    string m_desorptionSpeciesName;

    /// Species instance that can desorb
    Species* m_desorptionSpecies;

    /// Energy coefficients
    double m_desorptionEnergy;

    /// Frequency
    double m_desorptionFrequency;

    /// Pointer to associated adsorption class
    Adsorption* m_pAdsorption;

    /// Pointer to associated diffusion class
    Diffusion* m_pDiffusion;

    // Build probability table 
    vector<double> generateProbabilities();
    
    // Vector to hold the probabilities. Number of neighbour - 1 = index of list
    // TODO: How to initialize this as a const vector? The value should not change 
    vector<double> m_probabilities;

    // Number of sites with n number of neighbours
    vector<double> m_numNeighbours;

    // Maximum number of neighbours possible
    const int m_maxNeighbours;

};
}

#endif // ADSORPTION_H
