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

#ifndef DIFFUSION_H
#define DIFFUSION_H

#include "process.h"

/** The diffusion process. Performs the movement
 * of a particle to diffrent positions on the surface. */

namespace MicroProcesses{

class Diffusion: public Process
{
public:
    /// Constructor
    Diffusion();

    /// Destructor
    virtual ~Diffusion();

    /// Here everyhing needed for initialization of the process should be put
    void init();

    /// The name of the process.
    void setName(string s);

    /// Compute the overall probabilities of this process and return it.
    double getProbability();

    /// Returns the name of the process.
    string getName();

    /// Set the instance of Apothesis.
    /// This allows to have access to all other functionalities of the KMC class.
    void setInstance( Apothesis* apothesis ){ m_apothesis = apothesis; }

    /// Set the lattice uppon which adsorpiton will be performed.
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

  protected:
    /// The kmc instance.
    Apothesis* m_apothesis;

private:
    REGISTER_PROCESS( Diffusion)
};

}
#endif // DIFFUSION_H
