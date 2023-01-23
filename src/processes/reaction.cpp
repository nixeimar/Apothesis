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

Reaction::Reaction(): m_bLeadsToGrowth(false){}
Reaction::~Reaction(){}

void Reaction::init(vector<string> params){
    //Here the params of this process are set and the probability is calcylated (either directly or though calling to a function.
    m_vParams = params;

    if ( m_vProducts.size() == 0 || m_vReactants.size() == 0 || m_mProducts.size() == 0 || m_mReactants.size() == 0){
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

    //Holds species info in the case of 1-1 reaction e.g. A+B->AB
    m_idReacting.resize( m_pLattice->getSize() );

    vector<string> gSpecies = m_pUtilParams->getGrowthSpecies();

    buildTransformationMatrix();

    if ( gSpecies.size() <= 0 )
        m_bLeadsToGrowth = false;
    else {
        for ( string rs:m_vReactants) {
            for ( string gs:gSpecies ){
                if ( gs.compare( rs) == 0 ) {
                    m_bLeadsToGrowth = true;
                    break;
                }
            }
        }
    }

    if ( !m_bLeadsToGrowth ) {
        m_fRules = &Reaction::simpleRule;
        m_fPerform = &Reaction::catalysis;
    }
    else {

        if ( allReactCoeffOne() && m_vReactants.size() == 2 && m_vProducts.size() <= 2  ){
            m_fRules = &Reaction::oneOneRule;
            m_fPerform = &Reaction::oneOneReaction;
        }
    }
}

void Reaction::buildTransformationMatrix(){
    int iCount = 0;
    for (string r:m_vReactants ) {
        if ( iCount < m_vProducts.size() ) {
            m_mTransformationMatrix[ r ] = m_vProducts[ iCount ];
        }
        else
            m_mTransformationMatrix[ r ] = "";

        iCount++;
    }
}

bool Reaction::allReactCoeffOne(){
    for ( int i:m_vCoefReactants )
        if ( i != 1 )
            return false;

    return true;
}

void Reaction::constantType(){
    m_dProb = m_dReactionRate; //*m_pLattice->getSize();
}

void Reaction::arrheniusType(double v0, double Ed, double T)
{
    double k = m_pUtilParams->dkBoltz;
    Ed = Ed/m_pUtilParams->dAvogadroNum;

    m_dProb = v0*exp(-Ed/(k*T));
}

void Reaction::oneOneReaction( Site* s){

}

//Reaction A* + B* -> AB*
bool Reaction::oneOneRule( Site* s){
    if ( !s->isOccupied() ) return false;
    if ( !isReactant( s ) )
        return false;

    return true;
}

bool Reaction::simpleRule(Site* s){
    if ( !s->isOccupied() ) return false;

    if ( !isReactant( s ) ) return false;
    else {
        //Search for the other sites
        for ( pair<string, int> r:m_mReactants ){
            if ( r.first.compare( s->getLabel() ) != 0 ){
                for (Site* neigh:s->getNeighs() ){
                    if ( neigh->getLabel().compare( r.first ) == 0 ) {
                        m_idReacting[ s->getID() ].insert( neigh->getID() );
                        m_idReacting[ neigh->getID() ].insert( s->getID() );
                    }
                    else {
                        m_idReacting[ s->getID() ].erase( neigh->getID() );
                        m_idReacting[ neigh->getID() ].erase( s->getID() );
                    }
                }
            }
        }
    }

    if (  m_idReacting[ s->getID() ].size() != 0 )
        return true;
    return false;
}

bool Reaction::rules(Site *s)
{
    (this->*m_fRules)(s);
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
    (this->*m_fPerform)(s);
}

void Reaction::catalysis(Site *s){

    if ( s->getID() == 88 )
        cout << endl;

    string test = s->getLabel();

    s->setOccupied(false);
    s->setLabel( s->getBelowLabel() );

    m_seAffectedSites.insert( s );
    for ( Site* neigh:s->getNeighs() )
        m_seAffectedSites.insert( neigh );

    int lucky = m_pRandomGen->getIntRandom(0, m_idReacting[ s->getID() ].size() - 1 );

    int x = *std::next(m_idReacting[ s->getID() ].begin(), lucky );

    Site* otherSite =  m_pLattice->getSite( x );

    if ( !isReactant(otherSite )){

        cout << test << " " << otherSite->getLabel() << endl;
        cout << s->getID() << " " << otherSite->getID() << endl;

        cout << "Oh! Fuck!" << endl;

        otherSite->setOccupied( false );
        otherSite->setLabel( "X" );

        std::string name="SurfaceSpecies_" + std::to_string( 11111111 ) + ".dat";
        std::ofstream file(name);

        file << "Time (s): " << time << endl;

        for (int i = 0; i < m_pLattice->getY(); i++){
            for (int j = 0; j < m_pLattice->getX(); j++)
                file << m_pLattice->getSite( i*m_pLattice->getX() + j )->getLabel() << " " ;

            file << endl;
        }

        EXIT;
    }

    otherSite->setOccupied( false );
    otherSite->setLabel( otherSite->getBelowLabel() );

    m_seAffectedSites.insert( otherSite );
    for ( Site* neigh:otherSite->getNeighs() )
        m_seAffectedSites.insert( neigh );

    m_idReacting[  s->getID() ].erase( otherSite->getID() );
    m_idReacting[ otherSite->getID() ].erase( s->getID() );
}

double Reaction::getProbability(){ return m_dProb; }

