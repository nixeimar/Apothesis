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

#include "BCC.h"
#include "read.h"

#include <map>

BCC::BCC(Apothesis *apothesis) : Lattice(apothesis), m_iMinNeigs(1)
{
    ;
}

BCC::BCC(Apothesis *apothesis, bool step, vector<int> stepInfo) : Lattice(apothesis),
    m_hasSteps(step),
    m_stepInfo(stepInfo)
{
    ;
}

void BCC::buildSteps(int iSize, int jSize, int kSize)
{
    //e.g. Step 20 1 0

    if ( m_vSites.size()%iSize != 0){
        cout << "Cannot create  stepped surface because it cannot be divided exaclty." << endl;
        EXIT;
    }

    int iStep = 0;
    int iHeight = 0;

    //Steps in x-direction
    for (int i = 0; i < m_iSizeX; i++) {
        if ( iStep == iSize ) {
            iStep = 0;
            iHeight += jSize;
        }

        for (int j = 0; j < m_iSizeY; j++)
            getSite(j, i)->setHeight( m_iHeight + iHeight );

        iStep++;
    }

//    print();

    for ( int i = 0; i< m_vSites.size(); i++)
        calculateNeighNum( i );

  //  printNeighNum();
}

void BCC::setInitialHeight(int height) { m_iHeight = height; }

void BCC::build()
{
    if (m_Type == NONE)
    {
        cout << "Not supported lattice type" << endl;
        EXIT;
    }

    if (m_iSizeX == 0 || m_iSizeY == 0)
    {
        m_errorHandler->error_simple_msg("The lattice size cannot be zero in either dimension.");
        EXIT;
    }

    if (m_iHeight < 5)
    {
        m_errorHandler->warningSimple_msg("The lattice initial height is too small.Consider revising.");
    }

    // The sites of the lattice.
    m_vSites.resize( getSize() );
    for (int i = 0; i < m_vSites.size(); i++)
        m_vSites[i] = new Site();

    //This is OK
    for (int i = 0; i < m_iSizeX; i++)
    {
        for (int j = i * m_iSizeY; j < (m_iSizeY + i * m_iSizeY); j++)
        {
            m_vSites[j]->setID(j);
            m_vSites[j]->setHeight(m_iHeight - 1);
        }
    }

    mf_neigh();
}

BCC::~BCC()
{
    for (int i = 0; i < getSize(); i++)
        delete m_vSites[i];
}

void BCC::setSteps(bool hasSteps)
{
    m_bHasSteps = hasSteps;
}

void BCC::setStepInfo(int sizeX, int sizeY, int sizeZ)
{
    m_iStepX = sizeX;
    m_iStepY = sizeY;
    m_iStepZ = sizeZ;
}

void BCC::mf_buildSteps()
{

    if (m_iSizeX % m_iStepX != 0)
    {
        m_errorHandler->error_simple_msg("ERROR: The number of steps you provided doesn't conform with the lattice size ");
        exit(0);
    }
    if (m_iStepY != 0) // Be sure that we do have steps. If indi_y = 0 (1 0 0) then we have an initial flat surface
    {
        unsigned int steps = m_iSizeX / m_iStepX;
        for (unsigned int step = 1; step < steps; step++)
            for (unsigned int i = step * m_iStepX; i < (step + 1) * m_iStepX; i++)
                for (unsigned int j = 0; j < m_iSizeY; j++)
                    m_vSites[i * m_iStepY + j] += m_iStepY * step;
        //(*mesh)[i][j] += m_iStepY * step;///
        cout << "Number of steps:" << steps << endl;
    }
}

void BCC::mf_neigh()
{
    /* All except the boundaries */
    for (int i = 1; i < m_iSizeY - 1; i++){
        for (int j = 1; j < m_iSizeX - 1; j++){
            m_vSites[i * m_iSizeX + j]->setNeigh(m_vSites[(i - 1) * m_iSizeX + j]);
            m_vSites[i * m_iSizeX + j]->setNeigh(m_vSites[(i + 1) * m_iSizeX + j]);
            m_vSites[i * m_iSizeX + j]->setNeigh(m_vSites[i * m_iSizeX + j + 1]);
            m_vSites[i * m_iSizeX + j]->setNeigh(m_vSites[i * m_iSizeX + j - 1]);
        }
    }

    int iFirstCorner = 0;
    int iSecondCorner = m_iSizeX - 1;
    int iThirdCorner = m_iSizeX * m_iSizeY - m_iSizeX;
    int iForthCorner = m_iSizeX * m_iSizeY - 1;

    /*First row */
    for (int j = iFirstCorner; j <= iSecondCorner; j++){
        if (j != 0 && j != m_iSizeX - 1){
            m_vSites[j]->setNeigh(m_vSites[j - 1]);
            m_vSites[j]->setNeigh(m_vSites[j + 1]);
            m_vSites[j]->setNeigh(m_vSites[j + m_iSizeX]);
            m_vSites[j]->setNeigh(m_vSites[iThirdCorner + j]);
        }
        else if (j == iFirstCorner){
            m_vSites[j]->setNeigh(m_vSites[iSecondCorner]);
            m_vSites[j]->setNeigh(m_vSites[1]);
            m_vSites[j]->setNeigh(m_vSites[iSecondCorner + 1]);
            m_vSites[j]->setNeigh(m_vSites[iThirdCorner]);
        }
        else if (j == iSecondCorner){
            m_vSites[j]->setNeigh(m_vSites[j - 1]);
            m_vSites[j]->setNeigh(m_vSites[0]);
            m_vSites[j]->setNeigh(m_vSites[2 * m_iSizeX - 1]);
            m_vSites[j]->setNeigh(m_vSites[iForthCorner]);
        }
    }

    /*Last row */
    int iPos = 1;
    for (int j = iThirdCorner; j <= iForthCorner; j++){
        if (j != iThirdCorner && j != iForthCorner){
            m_vSites[j]->setNeigh(m_vSites[j - 1]);
            m_vSites[j]->setNeigh(m_vSites[j + 1]);
            m_vSites[j]->setNeigh(m_vSites[iFirstCorner + iPos]);
            m_vSites[j]->setNeigh(m_vSites[j - m_iSizeX]);
            iPos++;
        }
        else if (j == iThirdCorner){
            m_vSites[j]->setNeigh(m_vSites[iForthCorner]);
            m_vSites[j]->setNeigh(m_vSites[iThirdCorner + 1]);
            m_vSites[j]->setNeigh(m_vSites[iFirstCorner]);
            m_vSites[j]->setNeigh(m_vSites[iThirdCorner - m_iSizeX]);
        }
        else if (j == iForthCorner){
            m_vSites[j]->setNeigh(m_vSites[iForthCorner - 1]);
            m_vSites[j]->setNeigh(m_vSites[iThirdCorner]);
            m_vSites[j]->setNeigh(m_vSites[iSecondCorner]);
            m_vSites[j]->setNeigh(m_vSites[iThirdCorner - 1]);
        }
    }

    /* First column */
    for (int j = iFirstCorner + m_iSizeX; j < iThirdCorner; j += m_iSizeX){
        m_vSites[j]->setNeigh(m_vSites[j + m_iSizeX - 1]);
        m_vSites[j]->setNeigh(m_vSites[j + 1]);
        m_vSites[j]->setNeigh(m_vSites[j + m_iSizeX]);
        m_vSites[j]->setNeigh(m_vSites[j - m_iSizeX]);
    }

    /* Last column */
    for (int j = iSecondCorner + m_iSizeX; j < iForthCorner; j += m_iSizeX){
        m_vSites[j]->setNeigh(m_vSites[j - 1]);
        m_vSites[j]->setNeigh(m_vSites[j - m_iSizeX + 1]);
        m_vSites[j]->setNeigh(m_vSites[j + m_iSizeX]);
        m_vSites[j]->setNeigh(m_vSites[j - m_iSizeX]);
    }

    /*	int iCount = 0;
    int pos = 0;
    while (iCount < 100) {
        cout << "Enter pos to print neighbours: ";
        cin >> pos;
        cout << m_vSites[pos]->getID() << ": " << endl;
//		for (int i = 0; i < 4; i++) {
            cout << "WEST: " << m_vSites[pos]->getNeighPosition( Site::WEST )->getID()  << endl;
            cout << "EAST: " << m_vSites[pos]->getNeighPosition(Site::EAST)->getID() << endl;
            cout << "NORTH: " << m_vSites[pos]->getNeighPosition(Site::NORTH)->getID() << endl;
            cout << "SOUTH: " << m_vSites[pos]->getNeighPosition(Site::SOUTH)->getID() << endl;
            //		}
    }*/
}

//Site* BCC::getSite(int id) { return m_vSites[ id ]; }

void BCC::check()
{
    cout << "Checking lattice..." << endl;

    int test = 2;
    cout << test << ": ";
    cout << "W:" << getSite(test)->getNeighPosition(Site::WEST)->getID() << " ";\
    cout << "E:" << getSite(test)->getNeighPosition(Site::EAST)->getID() << " ";
    cout << "N:" << getSite(test)->getNeighPosition(Site::NORTH)->getID() << " ";
    cout << "S:" << getSite(test)->getNeighPosition(Site::SOUTH)->getID() << endl;
}


int BCC::calculateNeighNum( int id )
{
    int neighs = 1;
    for ( Site* s:m_vSites[ id ]->getNeighs() ) {
        if ( s->getHeight() >= m_vSites[ id ]->getHeight() )
            neighs++;
    }

    // THIS IS BAD! REFACTOR ....
    m_vSites[ id ]->setNeighsNum( neighs );
    return neighs;
}

// We have to define the processes e.g. Adsorption Simple, Asorption MultipleSites, Diffusion 1s Neighbors,
// Desorption 1st Neighbors, Reaction Decomposition etc...

/*void BCC::adsorp( int siteID, species_new* chemSpec )
{
    // >--------  For Lam & Vlachos (2000) ------------------------------------//
    //Remove site and its neihbors from its previous position in diffusion and desorption classes
    for ( auto &p:*m_pProcMap ) {
        if ( !IO::contains( p.first, "Adsoprtion" ) ){
            p.second.erase( siteID );

            for ( Site* s:m_vSites[ siteID ]->getNeighs() ) {
                p.second.erase( s->getID() );
                for ( Site* firstNeigh:s->getNeighs() )
                    p.second.erase( firstNeigh->getID() );
            }
        }
    }

    //For PVD results
    m_vSites[ siteID ]->increaseHeight( 1 );
    m_iSiteNeighsNum = calculateNeighNum( siteID );

    string strProc = "Desorption " + to_string( m_iSiteNeighsNum ) + "N";
    m_pProcMap->at( strProc ).insert( siteID );

    strProc = "Diffusion " + to_string( m_iSiteNeighsNum ) + "N";
    m_pProcMap->at( strProc ).insert( siteID );

    for ( Site* s:m_vSites[ siteID ]->getNeighs() ) {
        m_iSiteNeighsNum = calculateNeighNum( s->getID() );
        string strProc = "Desorption " + to_string( m_iSiteNeighsNum ) + "N";
        m_pProcMap->at( strProc ).insert( s->getID() );
        strProc = "Diffusion " + to_string( m_iSiteNeighsNum ) + "N";
        m_pProcMap->at( strProc ).insert( s->getID()  );

        for ( Site* firstNeigh:s->getNeighs() ){
            m_iSiteNeighsNum = calculateNeighNum( firstNeigh->getID() );
            string strProc = "Desorption " + to_string( m_iSiteNeighsNum ) + "N";
            m_pProcMap->at( strProc ).insert( firstNeigh->getID() );
            strProc = "Diffusion " + to_string( m_iSiteNeighsNum ) + "N";
            m_pProcMap->at( strProc ).insert( firstNeigh->getID()  );
        }
    }
    // ---------  For Lam & Vlachos (2000) ------------------------------------<//


    //1.Add to this site species the new adroped species
    //m_vSites[ siteID ]->getReactSpecies()[ chemSpec->getID() ]++;

    //2.Check if this site can react (according to the reactio matrix) and put it in the appropriate process map classes.
    //  2a. If a reaction with the max coeff has been completed. Remove it from adsoprtion and wait to react.

    //3.Check if this site (according to the reactio matrix) can accept diffused species from the neighbour sites or be a donor.
    //  3a. Do that for the neighbour sites and remove/add them accordingly in the process map.
}

void BCC::desorp(int siteID, species_new *chemSpecies)
{
    // <--------  For Lam & Vlachos (2000) ------------------------------------//
    //Remove site and its neihbors from its previous position in diffusion and desorption classes
    for ( auto &p:*m_pProcMap ) {
        if ( !IO::contains( p.first, "Adsoprtion" ) ){
            p.second.erase( siteID );

            for ( Site* s:m_vSites[ siteID ]->getNeighs() ) {
                p.second.erase( s->getID() );
                for ( Site* firstNeigh:s->getNeighs() )
                    p.second.erase( firstNeigh->getID() );
            }
        }
    }

    //For PVD results
    m_vSites[ siteID ]->decreaseHeight( 1 );
    m_iSiteNeighsNum = calculateNeighNum( siteID );

    string strProc = "Desorption " + to_string( m_iSiteNeighsNum ) + "N";
    m_pProcMap->at( strProc ).insert( siteID );

    strProc = "Diffusion " + to_string( m_iSiteNeighsNum ) + "N";
    m_pProcMap->at( strProc ).insert( siteID );

    for ( Site* s:m_vSites[ siteID ]->getNeighs() ) {
        m_iSiteNeighsNum = calculateNeighNum( s->getID() );
        string strProc = "Desorption " + to_string( m_iSiteNeighsNum ) + "N";
        m_pProcMap->at( strProc ).insert( s->getID() );
        strProc = "Diffusion " + to_string( m_iSiteNeighsNum ) + "N";
        m_pProcMap->at( strProc ).insert( s->getID()  );

        for ( Site* firstNeigh:s->getNeighs() ){
            m_iSiteNeighsNum = calculateNeighNum( firstNeigh->getID() );
            string strProc = "Desorption " + to_string( m_iSiteNeighsNum ) + "N";
            m_pProcMap->at( strProc ).insert( firstNeigh->getID() );
            strProc = "Diffusion " + to_string( m_iSiteNeighsNum ) + "N";
            m_pProcMap->at( strProc ).insert( firstNeigh->getID()  );

        }
    }
    // ---------  For Lam & Vlachos (2000) ------------------------------------>//


    // <--------  For Lam & Vlachos (2000) ------------------------------------//
    //Remove site from its previous positin in diffusion and desorption classes
/*    for ( auto &p:*m_pProcMap ) {
        if ( !IO::contains( p.first, "Adsoprtion" ) )
            p.second.erase( siteID );
    }

    this->print();

    //For PVD results
    m_vSites[ siteID ]->decreaseHeight();

    //This is for the basic site
    int neighs = 1;
    for ( Site* s:m_vSites[ siteID ]->getNeighs() ){
        if ( s->getHeight() <= m_vSites[ siteID ]->getHeight() )
            m_iMinNeigs++;
        else {
            for ( auto &p:*m_pProcMap ) {
                if ( !IO::contains( p.first, "Adsoprtion" ) )
                    p.second.erase( s->getID() );
            }

            s->decreaseNeighsNum();

            string strProcN = "Desorption " + to_string( s->getNeighboursNum() ) + "N";
            m_pProcMap->at( strProcN ).insert( s->getID() );

            strProcN = "Diffusion " + to_string( s->getNeighboursNum() ) + "N";
            m_pProcMap->at( strProcN ).insert( s->getID() );
        }
    }
    m_vSites[ siteID ]->setNeighboursNum( m_iMinNeigs );

    string strProc = "Desorption " + to_string( m_iMinNeigs ) + "N";
    m_pProcMap->at( strProc ).insert( siteID );

    strProc = "Diffusion " + to_string( m_iMinNeigs ) + "N";
    m_pProcMap->at( strProc ).insert( siteID );*/

// ---------  For Lam & Vlachos (2000) ------------------------------------>//

/*}

void BCC::react(int siteID)
{
    //1.React means increase the sites height by one.

    //2.Remove every species from this site (that means that the site is occupied by a lattice site).

    //3.Update neighbours and process map according to the new height as we would do in PVD.
}*/
