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

#ifndef ADSORPTION_H
#define ADSORPTION_H

#include "process.h"

#include <iostream>
#include <list>
#include <math.h>

#include "apothesis.h"
#include "desorption.h"
#include "diffusion.h"
#include "site.h"

using namespace std;
using namespace SurfaceTiles;

namespace MicroProcesses{

/** The adsorption class. Adsoprtion depending on the type of lattice can be performed in different sites.
For the simplest case e.g. BCC lattices all the sites are available for deposition. */

class Adsorption: public Process
{
  public:
    /// Constructor
    Adsorption
    (
      Apothesis* instance,
      string speciesName,
      Species* species,
      double stickingCoeffs,
      double massFraction
    );

    /// Destructor
    virtual ~Adsorption();

    /// The name of the process.
    void setName(string s){ m_sName = s;}

    /// Compute the overall probabilities of this processus and return it.
    double getProbability();

    /// Returns the name of the process.
    string getName();

    /// Returns name of spceies
    string getSpeciesName();

    /// Set the instance of the kmc.
    /// This allows to have access to all other functionalities of the KMC class.
    void setInstance( Apothesis* apothesis ){} //m_apothesis = apothesis; }

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

    /// Set associated desorption pointer
    void setDesorptionPointer(Desorption* d);

    /// Set Desorption switch
    void setDesorption(bool canDesorb);

    /// Set associated diffusion pointer
    void setDiffusionPointer(Diffusion* d);

    /// Set Diffusion switch
    void setDiffusion(bool canDiffuse);

    /// Set site
    void setSite(Site* s);

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
    int mf_updateNeighNum();

    const double getMassFraction();

    /// The site that adsorption is performed
    Site* m_site;

    /// The value of the probability of the process is stored here
    double m_dProbability;

    /// The number of neighs of this site
    int m_iNeighNum;

    /// Species that can adsorb
    string m_adsorptionSpeciesName;

    /// Instance of species class that can be adsorbed
    Species* m_adsorptionSpecies;

    /// Sticking coefficients
    double m_stickingCoeffs;
    
    /// Mass fractions
    double m_massfraction;

    /// Return pointer to corresponding desorption class
    Desorption* getDesorption();

    bool m_canDesorb;

    bool canDesorb();

    /// Return pointer to corresponding diffusion class
    Diffusion* getDiffusion();

    bool m_canDiffuse;
    
    bool canDiffuse();

  private:
  
    /** The lattice of the process */
    Lattice* m_pLattice;

    /// The adsorption list which hold all the available sites for adsorption
    list<Site* > m_lAdsSites;

    /// Pointer to associated desorption class
    Desorption* m_pDesorption;

    /// Pointer to associated diffusion class
    Diffusion* m_pDiffusion;


  };
}

#endif // ADSORPTION_H
