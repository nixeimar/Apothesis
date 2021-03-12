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

#include "lattice.h"
#include "read.h"

Lattice::Lattice(Apothesis *apothesis) : Pointers(apothesis)
{

}

void Lattice::setType(string sType)
{
    if (sType == "FCC")
        m_Type = FCC;
    else if (sType == "BCC")
        m_Type = BCC;
    else
        m_Type = NONE;
}

void Lattice::setX(int x) { m_iSizeX = x; }

void Lattice::setY(int y) { m_iSizeY = y; }

void Lattice::setInitialHeight(int height) { m_iHeight = height; }

Lattice::~Lattice()
{
}

vector<Site *> Lattice::getSites()
{
    return m_vSites;
}

Lattice::Type Lattice::getType()
{
    switch (m_Type)
    {
    case FCC:
        return FCC;
    case BCC:
        return BCC;
    default:
        return NONE;
    }
}

Site* Lattice::getSite(int id) { return m_vSites[id]; }

Site* Lattice::getSite(int i, int j)
{
    if ( i >= m_iSizeX ){
        m_errorHandler->error_simple_msg("ERROR: The site's index exceeds the size of the x-dimension of the lattice. ");
        exit(1);
    }

    if ( j >= m_iSizeY ){
        m_errorHandler->error_simple_msg("ERROR: The site's index exceeds the size of the y-dimension of the lattice. ");
        exit(1);
    }

    return m_vSites[ i*m_iSizeY + j ];
}

void Lattice::setProcMap( map<string, set< int > >* procMap )
{
    if ( procMap )
        m_pProcMap = procMap;
    else {
        m_errorHandler->warningSimple_msg("Process map could not be created.");
        exit(-1);
    }
}

void Lattice::print()
{
    for (int i = 0; i < m_iSizeY; i++){
        for (int j = 0; j < m_iSizeX; j++)
            cout << m_vSites[ i*m_iSizeX + j ]->getID() << "( " << m_vSites[ i*m_iSizeX + j ]->getHeight() << " )" ;

        cout  << endl;
    }
}


void Lattice::printNeighs( int ID )
{
    cout << "======= Printing neigbors ============ " << endl;

    if ( ID < getSize() ){
        cout << "Level 0 neighs: ";
        for ( int i = 0; i< m_vSites[ ID ]->get1stNeihbors()[ 0 ].size(); i++ )
            cout << m_vSites[ ID ]->get1stNeihbors()[ 0 ].at( i )->getID() << " ";

        cout << endl;

        cout << "Level -1 neighs: ";
        for ( int i = 0; i< m_vSites[ ID ]->get1stNeihbors()[ -1 ].size(); i++ )
            cout << m_vSites[ ID ]->get1stNeihbors()[ -1 ].at( i )->getID() << " ";

        cout << endl;

        cout << "======= end printing neigbors ============ " << endl;

    }
    else {
        cout << "Cannot print neighs because ID exceeds the available number of sites." << endl;
        EXIT;
    }

}
