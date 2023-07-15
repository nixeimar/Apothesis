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
        m_iNumNeighs = stoi( m_vParams[3] );
        arrhenius( stod(m_vParams[ 1 ]), stod(m_vParams[ 2 ]), stod(m_vParams[ 3 ]), m_pUtilParams->getTemperature(), m_iNumNeighs+1 );
    }
    else {
        m_error->error_simple_msg("Not supported type of process -> " + m_sProcName + " | " + m_sType );
        EXIT
    }

    //Create the rule for the adsoprtion process.
    if ( m_bAllNeihs )
        m_fRules = &Diffusion::mf_allRule;
    else
        m_fRules = &Diffusion::mf_basicRule;

    //Check what process should be performed.
    //Desorption in PVD will lead to increasing the height of the site
    //Desorption in CVD will change the label of the site
    if ( isPartOfGrowth() )
        m_fPerform = &Diffusion::mf_performPVD;
    else
        m_fPerform = &Diffusion::mf_performCVDALD;

    cout << endl;
}

bool Diffusion::isPartOfGrowth(){
    if (std::find(m_pUtilParams->getGrowthSpecies().begin(), m_pUtilParams->getGrowthSpecies().end(), m_sDiffused ) != m_pUtilParams->getGrowthSpecies().end())
        return true;

    return false;
}

void Diffusion::arrhenius(double v0, double E, double Em, double T,  int n)
{
    /*--- Taken from  Lam and Vlachos (2000)PHYSICAL REVIEW B, VOLUME 64, 035401 - DOI: 10.1103/PhysRevB.64.035401 ---*/
    /*    double Na = 6.0221417930e+23;				// Avogadro's number [1/mol]
    double P = 101325;					// [Pa]
    double T = any_cast<double>(m_vParams[0]); //500;						// [K]
    double k = any_cast<double>(m_vParams[1]); // 1.3806503e-23;			// Boltzmann's constant [j/K]
    double s0 = 0.1;
    double C_tot = 1.0e+19;				// [sites/m^2] Vlachos code says [moles sites/m^2]
    double E_d = any_cast<double>(m_vParams[2]); //(7.14e+4)/Na;			// [j]
    double E = 71128/Na;   //(7.14e+4)/Na;			// [j] -> 17 kcal
    double m = 32e-3/Na;				// [kg]
    double E_m = any_cast<double>(m_vParams[3]); //(4.28e+4)/Na;			// [j]
    double k_d = 1.0e+13;				// [s^-1]
    double y = 2.0e-3;					// Mole fraction of the precursor on the wafer
    /*--------------------------------------------------*/

    //   double v0 = k_d; //*exp(-E/(k*T));
    //   double A = exp( (E_d-E_m)/(k*T) );

    //--------------------- Transitions probability ----------------------------------------//
    //  return 0;// A*v0*exp( -(double)any_cast<int>(m_mParams["neighs"])*E/(k*T) );
    //----------------------------------------------------------------------------------------//

    double k = m_pUtilParams->dkBoltz;
    E = E/m_pUtilParams->dAvogadroNum;
    Em = Em/m_pUtilParams->dAvogadroNum;
    double A = exp(E-Em)/(k*T);
    m_dProb = v0*A*exp(-(double)n*E/(k*T));
}

bool Diffusion::mf_allRule(Site* s){
    if ( s->getNeighsNum() == m_iNumNeighs )
        return true;
    return false;}

bool Diffusion::mf_basicRule(Site* s){
    return true;
}

void Diffusion::mf_performCVDALD( Site* s){
}


void Diffusion::mf_diffusionSingleAtom(Site* s){

    mf_calculateNeighbors(s);
    vector<Site*>vacantSites;

    Site* neighWest=s->Site::getNeighPosition(Site::WEST);
    if(!neighWest->isOccupied()){
        vacantSites.push_back(neighWest);
    }

    Site* neighEast=s->Site::getNeighPosition(Site::EAST);
     if(!neighWest->isOccupied()){
        vacantSites.push_back(neighWest);
    }

    Site* neighNorth=s->Site::getNeighPosition(Site::NORTH);
     if(!neighWest->isOccupied()){
        vacantSites.push_back(neighWest);
    }

    Site* neighSouth=s->Site::getNeighPosition(Site::SOUTH);
     if(!neighWest->isOccupied()){
        vacantSites.push_back(neighWest);
    }
    
    if(vacantSites.size()==1){

     Site* newAdsorbSite=vacantSites[0];

    }
    else if(vacantSites.size()==2){

    Site* newAdsorbSite;
    if ( m_pRandomGen )
        newAdsorbSite= vacantSites[( m_pRandomGen->getIntRandom(0,1))];
    else{
        cout << "The random generator has not been defined." << endl;
        EXIT
    }

    }
    else if(vacantSites.size()==3){
           Site* newAdsorbSite;
    if ( m_pRandomGen )
        newAdsorbSite= vacantSites[( m_pRandomGen->getIntRandom(0,2))];
    else{
        cout << "The random generator has not been defined." << endl;
        EXIT
    }

    }
    else if(vacantSites.size()==4){
           Site* newAdsorbSite;
    if ( m_pRandomGen )
        newAdsorbSite= vacantSites[( m_pRandomGen->getIntRandom(0,3))];
    else{
        cout << "The random generator has not been defined." << endl;
        EXIT
    }

    }
    else{
        //doNothing
    }

}

void Diffusion::mf_performPVD( Site* s){
    //----- This is desorption ------------------------------------------------------------->
    s->decreaseHeight( 1 );
    mf_calculateNeighbors( s ) ;
    m_seAffectedSites.insert( s );
    for ( Site* neigh:s->getNeighs() ) {
        mf_calculateNeighbors( neigh );
        m_seAffectedSites.insert( neigh );

        for ( Site* firstNeigh:neigh->getNeighs() ){
            firstNeigh->setNeighsNum( mf_calculateNeighbors( firstNeigh ) );
            m_seAffectedSites.insert( firstNeigh );
        }
    }
    //--------------------------------------------------------------------------------------<

    // Random pick a site to re-adsorpt
    Site* adsorbSite;
    if ( m_pRandomGen )
        adsorbSite = s->getNeighs().at( m_pRandomGen->getIntRandom(0, m_iNumNeighs-2) );
    else{
        cout << "The random generator has not been defined." << endl;
        EXIT
    }

    //----- This is adsoprtion ------------------------------------------------------------->
    s->increaseHeight( 1 );
    mf_calculateNeighbors( s );
    m_seAffectedSites.insert( s ) ;

    for ( Site* neigh:s->getNeighs() ) {
        mf_calculateNeighbors( neigh );
        m_seAffectedSites.insert( neigh ) ;
    }
    //--------------------------------------------------------------------------------------<
}

void Diffusion::perform( Site* s)
{
    (this->*m_fPerform)(s);
}

int Diffusion::mf_calculateNeighbors(Site* s)
{
    int neighs = 0;

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

bool Diffusion::rules( Site* s)
{
    (this->*m_fRules)(s);
}

double Diffusion::getRateConstant(){ m_dProb; }

}
