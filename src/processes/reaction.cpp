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

#include "reaction.h"

Reaction::Reaction(){}
Reaction::~Reaction(){}

void Reaction::init(vector<string> params){
    //Here the params of this process are set and the probability is calcylated (either directly or though calling to a function.
    m_vParams = params;

    if ( m_mProducts.size() == 0 || m_mReactants.size() == 0){
        m_error->error_simple_msg("No reactants or products defined for reaction | " + m_sProcName );
        EXIT
    }

    //In the first must always be the type
    m_sType = any_cast<string>(m_vParams[ 0 ]);

    if ( m_sType.compare("arrhenius") == 0 ){
        arrheniusType( stod(m_vParams[ 1 ]), stod(m_vParams[ 2 ]), m_pUtilParams->getTemperature() );
    }
    else if (m_sType.compare("constant") == 0){
        m_dReactionRate = stod( m_vParams[1] );
        constantType();
    }
    else {
        m_error->error_simple_msg("Not supported type of process -> " + m_sProcName + " | " + m_sType );
        EXIT
    }

    //Create the rule for the adsoprtion process.
/*    if ( m_bAllNeihs &&  isPartOfGrowth( m_sDesorbed ) )
        m_fRules = &Desorption::mf_allRule;
    else if ( !m_bAllNeihs &&  isPartOfGrowth( m_sDesorbed ) )
        m_fRules = &Desorption::mf_basicRule;
    else
        m_fRules = &Desorption::mf_difSpeciesRule;*/
}


void Reaction::constantType(){}


void Reaction::arrheniusType(double v0, double Ed, double T)
{
    double k = m_pUtilParams->dkBoltz;
    Ed = Ed/m_pUtilParams->dAvogadroNum;

    m_dProb = v0*exp(-Ed/(k*T));
}

bool Reaction::rules(Site *s)
{
    if ( !s->isOccupied() ) return false;
    if ( !isReactant( s ) )
        return false;
    else {

        // First search for the site
        int stoichio = m_mReactants[ s->getLabel() ];

        if ( stoichio == 1) {
            //Search for the other sites
            for ( pair<string, int> r:m_mReactants ){

                if ( r.first.compare( s->getLabel() ) != 0 ){

                    // First search for the site
                    double stoichio = r.second;

                    double iCount = 0;
                    for (Site* neigh:s->getNeighs() ){
                        if ( neigh->isOccupied() ) {
                            if ( neigh->getLabel().compare( r.first ) == 0 ) {
                                iCount++;

                                if ( iCount == stoichio )
                                    break;
                            }
                        }
                    }

                    if ( iCount != stoichio )
                        return false;
                }
            }
        }
    }

    return true;
}

bool Reaction::isReactant(Site* s){
    for ( pair<string, double> p:m_mReactants ) {
        if ( p.first.compare( s->getLabel() ) == 0 )
            return true;
    }

    return false;
}

void Reaction::perform(Site *s)
{


}

double Reaction::getProbability(){
    return 0.;
}

