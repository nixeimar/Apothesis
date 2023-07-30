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

#include "diffusion.h"

namespace MicroProcesses
{

REGISTER_PROCESS_IMPL(Diffusion)

Diffusion::Diffusion(){}
Diffusion::~Diffusion(){}


void Diffusion::init(vector<string> params)
{
    //Here the params of this process are set and the probability is calcylated (either directly or though calling to a function.
    m_vParams = params;

    //In the first must always be the type
    m_sType = any_cast<string>(m_vParams[ 0 ]);
    if ( m_sType.compare("arrhenius") == 0 ){

        m_dv0 = stod(m_vParams[ 0 ]);
        m_dEd = stod(m_vParams[ 1 ]);
        m_dEdm = stod(m_vParams[ 2 ]);
        m_iNumNeighs = stoi( m_vParams[3] );

        m_fType = &arrheniusType;
    }
    else if ( m_sType.compare("constant") == 0 ){
        m_dDiffusionRate = stod(m_vParams[ 1 ]);

        m_fType = &constantType;
    }
    else {
        m_error->error_simple_msg("Not supported type of process -> " + m_sProcName + " | " + m_sType );
        EXIT
    }

    //Assign the type
    m_dRateConstant = (*m_fType)(this);

    m_isPartOfGrowth = isPartOfGrowth( m_sDiffused );

    //Select the rule for the diffusion process here
    if ( !m_isPartOfGrowth )
        if ( !m_bAllNeihs )
            m_fRules = &diffusionBasicRule;
        else
            m_fRules = &diffusionBasicAllRule;
    else
        m_fRules = &diffusionAllRule;

    //Check what process should be performed.
    //Desorption in PVD will lead to increasing the height of the site
    //Desorption in CVD will change the label of the site
    if ( !m_isPartOfGrowth )
        m_fPerform = &simpleDiffusion;
    else
        m_fPerform = &performPVD;
}

bool Diffusion::rules( Site* s)
{
    (*m_fRules)(this, s);
}

void Diffusion::perform( Site* s)
{
    m_seAffectedSites.clear();
    (*m_fPerform)(this, s);
}

int Diffusion::calculateSameNeighbors(Site* s){
    int neighs = 0;
    for ( Site* neigh:s->getNeighs() ) {
        if ( neigh->getLabel().compare(s->getLabel() ) == 0)
            neighs++;
    }
}

int Diffusion::calculateNeighbors(Site* s)
{
    int neighs = 0;

    if ( m_isPartOfGrowth ) {
        //If the particle of the lattice diffuses
        if ( m_pLattice->hasSteps() ){
            for ( Site* neigh:s->getNeighs() ) {
                if ( s->isLowerStep() && neigh->isHigherStep() ){
                    if ( neigh->getHeight() >= s->getHeight() + m_pLattice->getStepDiff() + 1 )
                        neighs++;
                }
                else if ( neigh->isLowerStep() && s->isHigherStep() ){
                    if ( neigh->getHeight() >= s->getHeight() - m_pLattice->getStepDiff() + 1 )
                        neighs++;
                }
                else {
                    if ( neigh->getHeight() >= s->getHeight() )
                        neighs++;
                }
            }

            s->setNeighsNum( neighs );
        }
        else {
            for ( Site* neigh:s->getNeighs() ) {
                if ( neigh->getHeight() >= s->getHeight() )
                    neighs++;
            }

            s->setNeighsNum( neighs );
        }
    } else {
        for ( Site* neigh:s->getNeighs() ) {
            if ( !neigh->isOccupied()  )
                neighs++;
        }
    }

    return neighs;
}

bool Diffusion::mf_isInLowerStep(Site* s)
{
    for (int j = 0; j < m_pLattice->getY(); j++)
        if ( s->getID() == m_pLattice->getSite( j, 0 )->getID() )
            return true;

    return false;
}

bool Diffusion::mf_isInHigherStep(Site* s)
{
    for (int j = 0; j < m_pLattice->getY(); j++){
        if ( s->getID() == m_pLattice->getSite( j, m_pLattice->getX() - 1 )->getID() ){
            return true;
        }
    }

    return false;
}

}
