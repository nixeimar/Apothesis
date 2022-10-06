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

#include <map>

BCC::BCC(Apothesis *apothesis) : Lattice(apothesis), m_iMinNeigs(1)
{
    ;
}

BCC::BCC(Apothesis *apothesis, bool step, vector<int> stepInfo) : Lattice(apothesis),
    m_bHasSteps(step),
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

    int h = getSite( m_iSizeY-1, m_iSizeX-1)->getHeight() ;
    m_iStepDiff = abs( getSite( m_iSizeY-1, m_iSizeX-1)->getHeight() - getSite( 0, 0 )->getHeight() ) + 1;

//    printNeighNum();

    for (int j = 0; j < m_iSizeY; j++){
        getSite( j, 0 )->setLowerStep( true );
        getSite( j, m_iSizeX - 1  )->setHigherStep( true );
    }
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

void BCC::mf_neigh()
{
    /* All except the boundaries */
    for (int i = 1; i < m_iSizeY - 1; i++)
    {
        for (int j = 1; j < m_iSizeX - 1; j++)
        {
            m_vSites[i * m_iSizeX + j]->setNeigh(m_vSites[(i - 1) * m_iSizeX + j]);
            m_vSites[i * m_iSizeX + j]->setNeighPosition(m_vSites[(i - 1) * m_iSizeX + j], Site::NORTH);

            m_vSites[i * m_iSizeX + j]->setNeigh(m_vSites[(i + 1) * m_iSizeX + j]);
            m_vSites[i * m_iSizeX + j]->setNeighPosition(m_vSites[(i + 1) * m_iSizeX + j], Site::SOUTH);

            m_vSites[i * m_iSizeX + j]->setNeigh(m_vSites[i * m_iSizeX + j + 1]);
            m_vSites[i * m_iSizeX + j]->setNeighPosition(m_vSites[i * m_iSizeX + j + 1], Site::EAST);

            m_vSites[i * m_iSizeX + j]->setNeigh(m_vSites[i * m_iSizeX + j - 1]);
            m_vSites[i * m_iSizeX + j]->setNeighPosition(m_vSites[i * m_iSizeX + j - 1], Site::WEST);
        }
    }

    int iFirstCorner = 0;
    int iSecondCorner = m_iSizeX - 1;
    int iThirdCorner = m_iSizeX * m_iSizeY - m_iSizeX;
    int iForthCorner = m_iSizeX * m_iSizeY - 1;

    /*First row */
    for (int j = iFirstCorner; j <= iSecondCorner; j++)
    {
        if (j != 0 && j != m_iSizeX - 1)
        {
            m_vSites[j]->setNeigh(m_vSites[j - 1]);
            m_vSites[j]->setNeighPosition(m_vSites[j - 1], Site::WEST);

            m_vSites[j]->setNeigh(m_vSites[j + 1]);
            m_vSites[j]->setNeighPosition(m_vSites[j + 1], Site::EAST);

            m_vSites[j]->setNeigh(m_vSites[j + m_iSizeX]);
            m_vSites[j]->setNeighPosition(m_vSites[j + m_iSizeX], Site::SOUTH);

            m_vSites[j]->setNeigh(m_vSites[iThirdCorner + j]);
            m_vSites[j]->setNeighPosition(m_vSites[iThirdCorner + j], Site::NORTH);
        }
        else if (j == iFirstCorner)
        {
            m_vSites[j]->setNeigh(m_vSites[iSecondCorner]);
            m_vSites[j]->setNeighPosition(m_vSites[iSecondCorner], Site::WEST);

            m_vSites[j]->setNeigh(m_vSites[1]);
            m_vSites[j]->setNeighPosition(m_vSites[1], Site::EAST);

            m_vSites[j]->setNeigh(m_vSites[iSecondCorner + 1]);
            m_vSites[j]->setNeighPosition(m_vSites[iSecondCorner + 1], Site::SOUTH);

            m_vSites[j]->setNeigh(m_vSites[iThirdCorner]);
            m_vSites[j]->setNeighPosition(m_vSites[iThirdCorner], Site::NORTH);
        }
        else if (j == iSecondCorner)
        {
            m_vSites[j]->setNeigh(m_vSites[j - 1]);
            m_vSites[j]->setNeighPosition(m_vSites[j - 1], Site::WEST);

            m_vSites[j]->setNeigh(m_vSites[0]);
            m_vSites[j]->setNeighPosition(m_vSites[0], Site::EAST);

            m_vSites[j]->setNeigh(m_vSites[2 * m_iSizeX - 1]);
            m_vSites[j]->setNeighPosition(m_vSites[2 * m_iSizeX - 1], Site::SOUTH);

            m_vSites[j]->setNeigh(m_vSites[iForthCorner]);
            m_vSites[j]->setNeighPosition(m_vSites[iForthCorner], Site::NORTH);
        }
    }

    /*Last row */
    int iPos = 1;
    for (int j = iThirdCorner; j <= iForthCorner; j++)
    {
        if (j != iThirdCorner && j != iForthCorner)
        {
            m_vSites[j]->setNeigh(m_vSites[j - 1]);
            m_vSites[j]->setNeighPosition(m_vSites[j - 1], Site::WEST);

            m_vSites[j]->setNeigh(m_vSites[j + 1]);
            m_vSites[j]->setNeighPosition(m_vSites[j + 1], Site::EAST);

            m_vSites[j]->setNeigh(m_vSites[iFirstCorner + iPos]);
            m_vSites[j]->setNeighPosition(m_vSites[iFirstCorner + iPos], Site::SOUTH);

            m_vSites[j]->setNeigh(m_vSites[j - m_iSizeX]);
            m_vSites[j]->setNeighPosition(m_vSites[j - m_iSizeX], Site::NORTH);
            iPos++;
        }
        else if (j == iThirdCorner)
        {
            m_vSites[j]->setNeigh(m_vSites[iForthCorner]);
            m_vSites[j]->setNeighPosition(m_vSites[iForthCorner], Site::WEST);

            m_vSites[j]->setNeigh(m_vSites[iThirdCorner + 1]);
            m_vSites[j]->setNeighPosition(m_vSites[iThirdCorner + 1], Site::EAST);

            m_vSites[j]->setNeigh(m_vSites[iFirstCorner]);
            m_vSites[j]->setNeighPosition(m_vSites[iFirstCorner], Site::SOUTH);

            m_vSites[j]->setNeigh(m_vSites[iThirdCorner - m_iSizeX]);
            m_vSites[j]->setNeighPosition(m_vSites[iThirdCorner - m_iSizeX], Site::NORTH);
        }
        else if (j == iForthCorner)
        {
            m_vSites[j]->setNeigh(m_vSites[iForthCorner - 1]);
            m_vSites[j]->setNeighPosition(m_vSites[iForthCorner - 1], Site::WEST);

            m_vSites[j]->setNeigh(m_vSites[iThirdCorner]);
            m_vSites[j]->setNeighPosition(m_vSites[iThirdCorner], Site::EAST);

            m_vSites[j]->setNeigh(m_vSites[iSecondCorner]);
            m_vSites[j]->setNeighPosition(m_vSites[iSecondCorner], Site::SOUTH);

            m_vSites[j]->setNeigh(m_vSites[iThirdCorner - 1]);
            m_vSites[j]->setNeighPosition(m_vSites[iThirdCorner - 1], Site::NORTH);
        }
    }

    /* First column */
    for (int j = iFirstCorner + m_iSizeX; j < iThirdCorner; j += m_iSizeX)
    {
        m_vSites[j]->setNeigh(m_vSites[j + m_iSizeX - 1]);
        m_vSites[j]->setNeighPosition(m_vSites[j + m_iSizeX - 1], Site::WEST);

        m_vSites[j]->setNeigh(m_vSites[j + 1]);
        m_vSites[j]->setNeighPosition(m_vSites[j + 1], Site::EAST);

        m_vSites[j]->setNeigh(m_vSites[j + m_iSizeX]);
        m_vSites[j]->setNeighPosition(m_vSites[j + m_iSizeX], Site::SOUTH);

        m_vSites[j]->setNeigh(m_vSites[j - m_iSizeX]);
        m_vSites[j]->setNeighPosition(m_vSites[j - m_iSizeX], Site::NORTH);
    }

    /* Last column */
    for (int j = iSecondCorner + m_iSizeX; j < iForthCorner; j += m_iSizeX)
    {
        m_vSites[j]->setNeigh(m_vSites[j - 1]);
        m_vSites[j]->setNeighPosition(m_vSites[j - 1], Site::WEST);

        m_vSites[j]->setNeigh(m_vSites[j - m_iSizeX + 1]);
        m_vSites[j]->setNeighPosition(m_vSites[j - m_iSizeX + 1], Site::EAST);

        m_vSites[j]->setNeigh(m_vSites[j + m_iSizeX]);
        m_vSites[j]->setNeighPosition(m_vSites[j + m_iSizeX], Site::SOUTH);

        m_vSites[j]->setNeigh(m_vSites[j - m_iSizeX]);
        m_vSites[j]->setNeighPosition(m_vSites[j - m_iSizeX], Site::NORTH);
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


void BCC::writeLatticeHeights( double time, int step )
{
    std::ofstream file("Lattice_" + to_string(step) );

    for (int i = 0; i < m_iSizeY; i++){
        for (int j = 0; j < m_iSizeX; j++)
            file << m_vSites[ i*m_iSizeX + j ]->getHeight() << " ";
        file  << endl;
    }


    for (int i = 0; i < m_iSizeY; i++){
        for (int j = 0; j < m_iSizeX; j++)
            file << m_vSites[ i*m_iSizeX + j ]->getLabel() + to_string( m_vSites[ i*m_iSizeX + j ]->getID() )  << "\t" << "( " << m_vSites[ i*m_iSizeX + j ]->getHeight() << " ) " ;
        file  << endl;
    }

    file.close();
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
