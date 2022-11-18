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

using namespace std;

namespace MicroProcesses
{

class Diffusion:public Process
{
public:
    Diffusion();
    ~Diffusion() override;

    inline void setOriginSite( Site* site ){ m_originSite = site;}
    inline Site* getOriginSite(){ return m_originSite; }

    inline void setTargetSite( Site* site ){ m_targetSite = site;}
    inline Site* getTargetSite(){ return m_targetSite; }

//    inline void setNeigh(int n ){ m_iNumNeighs = n; }

    double getProbability() override;
    bool rules( Site* ) override;
    void perform( Site* ) override;

    void init(vector<string> params) override;

    void arrhenius(double v0, double E, double Em, double T,  int n);

private:

    bool mf_isInLowerStep( Site* s );
    bool mf_isInHigherStep( Site* s );

    ///The site to for the adsorption to be removed
    Site* m_originSite;

    ///The site that adsorption will be performed
    Site* m_targetSite;

    /// A member function to calculate the neighbors of a given site
    int mf_calculateNeighbors(Site*);

    /// The number of neighbours for calculating the probability
    int m_iNumNeighs;

    REGISTER_PROCESS( Diffusion )
};

}

#endif // DiffusionSimpleCubic_H
