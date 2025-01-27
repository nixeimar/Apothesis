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

#include "SimpleCubic.h"

#include <map>

SimpleCubic::SimpleCubic(Apothesis *apothesis) : Lattice(apothesis), m_iMinNeigs(1)
{
    m_Type = Lattice::SimpleCubic;
}

void SimpleCubic::buildSteps()
{    
    int iPerStep = 0;
    if ( m_vSites.size()%m_iNumSteps != 0){
        cout << "Cannot create  stepped surface because it cannot be divided exaclty." << endl;
        EXIT
    }
    else {
        iPerStep = m_iSizeX/m_iNumSteps;
    }

    int iStep = 0;
    int iHeight = 0;

    //Steps in x-direction
    for (int i = 0; i < m_iSizeX; i++) {
        if ( iStep == iPerStep ) {
            iStep = 0;
            iHeight += m_iStepHeight;
        }

        for (int j = 0; j < m_iSizeY; j++)
            getSite(j, i)->setHeight( m_iHeight + iHeight );

        iStep++;
    }

    for ( int i = 0; i< m_vSites.size(); i++)
        calculateNeighNum( i );

    int h = getSite( m_iSizeY-1, m_iSizeX-1)->getHeight() ;
    m_iStepDiff = abs( getSite( m_iSizeY-1, m_iSizeX-1)->getHeight() - getSite( 0, 0 )->getHeight() ) + 1;

    for (int j = 0; j < m_iSizeY; j++){
        getSite( j, 0 )->setLowerStep( true );
        getSite( j, m_iSizeX - 1  )->setHigherStep( true );
    }
}

void SimpleCubic::readHeightsFromFile() {

    ifstream heightFile("heights.dat");

    if (heightFile.good()) {

        vector<vector<int>> heights;
        string line;

        // Read the file line by line
        while (std::getline(heightFile, line)) {

            if (line.empty() ) continue;

            std::vector<int> row;
            std::istringstream iss(line);
            int num;
            // Extract integers from the line and add them to the row
            while (iss >> num) {
                row.push_back(num);
            }
            // Add the row to the matrix
            heights.push_back(row);
        }

        int icount = 0;

        for (int i = 0; i < m_iSizeX; ++i) {
            for (int j = 0; j < m_iSizeY; ++j) {
                icount = i * m_iSizeY + j;
                m_vSites[icount]->setID(icount);
                m_vSites[icount]->setHeight( heights[ i ][ j ] );
            }
        }

    }
    heightFile.close();
}

void SimpleCubic::readSpeciesFromFile(){

    ifstream speciesFile("species.dat");

    if (speciesFile.good()) {

        vector<vector<string>> species;
        string line;

        // Read the file line by line
        while (std::getline(speciesFile, line)) {

            if (line.empty() ) continue;

            line = trim(line);

            std::vector<string> row;
            std::istringstream iss(line);
            string s;
            // Extract strings from the line and add them to the row
            while (iss >> s ) {
                row.push_back(s);
            }
            // Add the row to the matrix
            species.push_back(row);
        }

        int icount = 0;

        for (int i = 0; i < m_iSizeX; ++i) {
            for (int j = 0; j < m_iSizeY; ++j) {
                icount = i * m_iSizeY + j;
                m_vSites[icount]->setLabel( species[ i ][ j ] );

                if ( species[ i ][ j ].find("*") != std::string::npos)
                    m_vSites[ icount ]->setOccupied( true);
            }
        }
    }
    speciesFile.close();
}

void SimpleCubic::build()
{
    if (m_Type == NONE)
    {
        cout << "Not supported lattice type" << endl;
        EXIT
    }

    if (m_iSizeX == 0 || m_iSizeY == 0)
    {
        m_errorHandler->error_simple_msg("The lattice size cannot be zero in either dimension.");
        EXIT
    }


    mf_neigh();
}

SimpleCubic::~SimpleCubic()
{
    for (int i = 0; i < getSize(); i++)
        delete m_vSites[i];
}

void SimpleCubic::setSteps(bool hasSteps)
{
    m_bHasSteps = hasSteps;
}

void SimpleCubic::mf_neigh()
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

void SimpleCubic::check()
{
    cout << "Checking lattice..." << endl;

    int test = 2;
    cout << test << ": ";
    cout << "W:" << getSite(test)->getNeighPosition(Site::WEST)->getID() << " ";\
        cout << "E:" << getSite(test)->getNeighPosition(Site::EAST)->getID() << " ";
    cout << "N:" << getSite(test)->getNeighPosition(Site::NORTH)->getID() << " ";
    cout << "S:" << getSite(test)->getNeighPosition(Site::SOUTH)->getID() << endl;
}


void SimpleCubic::writeLatticeHeights( double time, int step )
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

int SimpleCubic::calculateNeighNum( int id )
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

unordered_map<string, double> SimpleCubic::computeCoverages( vector<string> species ) {
    for ( string name:species){
        m_mCoverages[ name ] = 0.;

        int iCount = 0;
        for ( int i =0; i< getSize(); i++){
            if ( m_vSites[ i ]->getLabel().compare( name ) == 0 )
                iCount++;
        }

        m_mCoverages[ name ] = (double)iCount/getSize();
    }

    return m_mCoverages;
}

