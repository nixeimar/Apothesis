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

#include "lattice.h"

Lattice::Lattice(Apothesis *apothesis) : Pointers(apothesis),m_iStepDiff(0)
{

}

void Lattice::setType(string sType)
{

    m_sType = sType;

    if (sType == "FCC")
        m_Type = FCC;
    else if (sType == "SimpleCubic")
        m_Type = SimpleCubic;
    else if (sType == "HCP")
        m_Type = HCP;
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

string Lattice::getTypeAsString(){ return  m_sType; }

Lattice::Type Lattice::getType()
{
    switch (m_Type)
    {
    case FCC:
        return FCC;
    case SimpleCubic:
        return SimpleCubic;
    case HCP:
        return HCP;
    default:
        return NONE;
    }
}

void Lattice::buildSteps(){;}


Site* Lattice::getSite(int id) { return m_vSites[id]; }

Site* Lattice::getSite(int i, int j)
{
    return m_vSites[ i*m_iSizeX + j ];
}

void Lattice::print()
{
    for (int i = 0; i < m_iSizeY; i++){
        for (int j = 0; j < m_iSizeX; j++)
            cout << m_vSites[ i*m_iSizeX + j ]->getLabel() + to_string( m_vSites[ i*m_iSizeX + j ]->getID() )  << "\t" << "( " << m_vSites[ i*m_iSizeX + j ]->getHeight() << " ) " ;
        cout  << endl;
    }
}

void Lattice::printNeighNum()
{
    for (int i = 0; i < m_iSizeY; i++){
        for (int j = 0; j < m_iSizeX; j++)
            cout << m_vSites[ i*m_iSizeX + j ]->getID() << "( " << m_vSites[ i*m_iSizeX + j ]->getNeighsNum() << " )" ;

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
        EXIT
    }

}

void Lattice::writeXYZ( string filename ){;}
void Lattice::writeLatticeHeights( double, int ){;}

void Lattice::printInfo() {

    cout << endl;
    cout << "--- start info lattice parameteres ----- " << endl;
    cout << "---------------------------------------- " << endl;
    cout << "Type: "; cout << getType() << endl;
    cout << "Size X: "; cout << getX() << endl;
    cout << "Size Y: "; cout << getY() << endl;
    cout << "Lattice species: "; cout << getLabels() << endl;

    if ( hasSteps() ) {
        cout << "Number of steps: "; cout << getNumSteps() << endl;
        cout << "Step height: "; cout << getStepHeight() << endl;
    }
    cout << "--- end lattice parameters info -------- " << endl;
    cout << endl;
}

