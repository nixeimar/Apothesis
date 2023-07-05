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

#include "HCP.h"

#include <map>
#include <unordered_map>

HCP::HCP(Apothesis *apothesis) : Lattice(apothesis), m_iMinNeigs(1)
{
    m_Type = Lattice::HCP;
}

void HCP::setInitialHeight(int height) { m_iHeight = height; }

void HCP::build()
{

    //CG: HPC resembles FCC but without the neighbors reside in difference height.
    //Example of a 6x6 lattice. The x-dim of the HCP lattice must be even number.

    // 0       2       4
    //     1       3       5
    // 6       8      10
    //     7       9      11
    // 12     14      16
    //    13      15      17
    // 18     20      22
    //    19      21      23
    // 24     26      28
    //    25      27      29
    // 30     32      34
    //    31      33      35

    // 0       2       4       6
    //     1       3       5       7
    // 8       10      12      14
    //     9       11      13      15
    // 16      18      20      22
    //     17      19      21      23
    // 24      26      28      30
    //     25      27      29       31
    // 32      34      36      38
    //     33      35      37      39


    // 0       2       4       6
    //     1       3       5       7
    // 8       10      12      14
    //     9       11      13      15
    // 16      18      20      22
    //     17      19      21      23
    // 24      26      28      30
    //     25      27      29      31
    // 32      34      36      38
    //     33      35      37      39
    // 40      42      44      46
    //     41      43      45      47
    // 48      50      52       54
    //     49      51      53      55
    // 56      58      60       62
    //     57      59      61      63


    // Example:
    // The neighbors for 14 are: 7 (West Up), 13 (West down), 9 (East up), 15 (East down), 8 (North), 20 (South)

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

    if (m_iSizeX%2 != 0 )
    {
        m_errorHandler->error_simple_msg("In HPC lattices, the x-dimension of the lattice must be an even number.");
        EXIT
    }


    if (m_iHeight < 5)
    {
        m_errorHandler->warningSimple_msg("The lattice initial height is too small.Consider revising.");
    }

    // The sites of the lattice.
    m_vSites.resize(getSize());
    for (int i = 0; i < m_vSites.size(); i++)
        m_vSites[i] = new Site();

    // This is OK
    for (int i = 0; i < m_iSizeX; i++)
    {
        for (int j = i * m_iSizeY; j < (m_iSizeY + i * m_iSizeY); j++)
        {
            m_vSites[j]->setID(j);
            m_vSites[j]->setHeight(m_iHeight);
        }
    }

    mf_neigh();

    for (Site* s:m_vSites[ 39 ]->getNeighs() ){
        cout << s->getID() << endl;
    }

    // Here we set the label of the species
    for (int i = 0; i < m_iSizeY; i++)
    {
        for (int j = 0; j < m_iSizeX; j++)
            m_vSites[i * m_iSizeX + j]->setLabel(m_sLabel);
    }
}

HCP::~HCP()
{
    for (int i = 0; i < getSize(); i++)
        delete m_vSites[i];
}

void HCP::mf_neigh()
{
    /* All except the boundaries */
    int iPos = 0;
    /* All except the boundaries */
    for ( int i = 1; i < m_iSizeY - 1; i++ ){
        for (int j = 2; j < m_iSizeX -2; j++) {
            iPos  = i*m_iSizeX + j;

            //Same level
            m_vSites[ iPos ]->setNeigh( m_vSites[ iPos - m_iSizeX ] ); //North
            m_vSites[ iPos ]->setNeighPosition( m_vSites[ iPos - m_iSizeX ],  Site::NORTH );

            m_vSites[ iPos ]->setNeigh( m_vSites[ iPos + m_iSizeX ]); //South
            m_vSites[ iPos ]->setNeighPosition( m_vSites[ iPos + m_iSizeX ],  Site::SOUTH);

            //The lattice x-dim is even
            if ( iPos%2 == 0){
                m_vSites[ iPos ]->setNeigh( m_vSites[ iPos - m_iSizeX + 1 ] );
                m_vSites[ iPos ]->setNeighPosition( m_vSites[ iPos - m_iSizeX + 1 ],  Site::EAST_UP );

                m_vSites[ iPos ]->setNeigh( m_vSites[ iPos - m_iSizeX - 1 ] );
                m_vSites[ iPos ]->setNeighPosition( m_vSites[ iPos - m_iSizeX - 1 ],  Site::WEST_UP );

                m_vSites[ iPos ]->setNeigh( m_vSites[ iPos - 1 ] );
                m_vSites[ iPos ]->setNeighPosition( m_vSites[ iPos - 1 ],  Site::WEST_DOWN );
                m_vSites[ iPos ]->setNeigh( m_vSites[ iPos + 1 ] );
                m_vSites[ iPos ]->setNeighPosition( m_vSites[ iPos + 1 ],  Site::EAST_DOWN );
            }
            else {
                m_vSites[ iPos ]->setNeigh( m_vSites[ iPos + 1 ] );
                m_vSites[ iPos ]->setNeighPosition( m_vSites[ iPos +1 ],  Site::EAST_UP );
                m_vSites[ iPos ]->setNeigh( m_vSites[ iPos - 1 ] );
                m_vSites[ iPos ]->setNeighPosition( m_vSites[ iPos - 1 ],  Site::WEST_UP );
                m_vSites[ iPos ]->setNeigh( m_vSites[ iPos + m_iSizeX + 1 ] );
                m_vSites[ iPos ]->setNeighPosition( m_vSites[ iPos + m_iSizeX  + 1 ],  Site::EAST_DOWN );
                m_vSites[ iPos ]->setNeigh( m_vSites[ iPos + m_iSizeX  - 1 ] );
                m_vSites[ iPos ]->setNeighPosition( m_vSites[ iPos + m_iSizeX  - 1 ],  Site::WEST_DOWN );
            }

        }
    }

    //Mark the corners. Helps to define the rest of the neighours
    // N1, N5 ---- N6, N2
    // |                |
    // |                |
    // N3, N8 ---- N7, N4

    // In the 6x6 example lattice (see above):
    // N1 = 0, N2 = 5, N3 = 30, N4 = 35
    // N5 = 1, N6 = 4, N7 = 34, N8 = 31

    int N1 = 0, N2 = m_iSizeX - 1, N3 = getSize() - m_iSizeX, N4 = getSize() - 1;
    int N5 = 1, N6 = N2 - 1, N7 = N4 - 1, N8 = N3 + 1;

    //N1
    m_vSites[ N1 ]->setNeigh( m_vSites[ m_iSizeX ] );
    m_vSites[ N1 ]->setNeighPosition( m_vSites[ iPos + m_iSizeX ],  Site::SOUTH );
    m_vSites[ N1 ]->setNeigh( m_vSites[ getSize() - m_iSizeX ] );
    m_vSites[ N1 ]->setNeighPosition( m_vSites[ iPos + m_iSizeX ],  Site::NORTH );
    m_vSites[ N1 ]->setNeigh( m_vSites[ N1 + 1 ] );
    m_vSites[ N1 ]->setNeighPosition( m_vSites[ N1 + 1 ],  Site::EAST_DOWN );
    m_vSites[ N1 ]->setNeigh( m_vSites[ N3 + 1 ] );
    m_vSites[ N1 ]->setNeighPosition( m_vSites[ N3 + 1 ],  Site::EAST_UP );
    m_vSites[ N1 ]->setNeigh( m_vSites[ N2 ] );
    m_vSites[ N1 ]->setNeighPosition( m_vSites[ N2 ],  Site::WEST_DOWN );
    m_vSites[ N1 ]->setNeigh( m_vSites[ N4 ] );
    m_vSites[ N1 ]->setNeighPosition( m_vSites[ N4 ],  Site::WEST_UP );

    //N2
    m_vSites[ N2 ]->setNeigh( m_vSites[ N4 ] );
    m_vSites[ N2 ]->setNeighPosition( m_vSites[ N4 ],  Site::NORTH );
    m_vSites[ N2 ]->setNeigh( m_vSites[ N2 + m_iSizeX ] );
    m_vSites[ N2 ]->setNeighPosition( m_vSites[ N2 + m_iSizeX ],  Site::SOUTH );
    m_vSites[ N2 ]->setNeigh( m_vSites[ N2 - 1 ] );
    m_vSites[ N2 ]->setNeighPosition( m_vSites[ N2 - 1 ],  Site::WEST_UP );
    m_vSites[ N2 ]->setNeigh( m_vSites[ N2 + m_iSizeX - 1 ] );
    m_vSites[ N2 ]->setNeighPosition( m_vSites[ N2 + m_iSizeX - 1 ],  Site::WEST_DOWN );
    m_vSites[ N2 ]->setNeigh( m_vSites[ m_iSizeX ] );
    m_vSites[ N2 ]->setNeighPosition( m_vSites[ m_iSizeX ],  Site::EAST_DOWN );
    m_vSites[ N2 ]->setNeigh( m_vSites[ N1 ] );
    m_vSites[ N2 ]->setNeighPosition( m_vSites[ N1 ],  Site::EAST_UP );

    //N3
    m_vSites[ N3 ]->setNeigh( m_vSites[ N1 ] );
    m_vSites[ N3 ]->setNeighPosition( m_vSites[ N1 ],  Site::SOUTH );
    m_vSites[ N3 ]->setNeigh( m_vSites[ N3 - m_iSizeX ] );
    m_vSites[ N3 ]->setNeighPosition( m_vSites[ N3 - m_iSizeX ],  Site::NORTH );
    m_vSites[ N3 ]->setNeigh( m_vSites[ N3 - m_iSizeX + 1 ] );
    m_vSites[ N3 ]->setNeighPosition( m_vSites[ N3 - m_iSizeX + 1 ],  Site::EAST_UP );
    m_vSites[ N3 ]->setNeigh( m_vSites[ N3 + 1 ] );
    m_vSites[ N3 ]->setNeighPosition( m_vSites[ N3 + 1 ],  Site::EAST_DOWN );
    m_vSites[ N3 ]->setNeigh( m_vSites[ N4 ] );
    m_vSites[ N3 ]->setNeighPosition( m_vSites[ N4 ],  Site::WEST_DOWN );
    m_vSites[ N3 ]->setNeigh( m_vSites[ N3 - 1 ] );
    m_vSites[ N3 ]->setNeighPosition( m_vSites[ N3 - 1 ],  Site::WEST_UP );

    //N4
    m_vSites[ N4 ]->setNeigh( m_vSites[ N2 ] );
    m_vSites[ N4 ]->setNeighPosition( m_vSites[ N2 ],  Site::SOUTH );
    m_vSites[ N4 ]->setNeigh( m_vSites[ N4 - m_iSizeX ] );
    m_vSites[ N4 ]->setNeighPosition( m_vSites[ N4 - m_iSizeX  ],  Site::NORTH );
    m_vSites[ N4 ]->setNeigh( m_vSites[ N3 ] );
    m_vSites[ N4 ]->setNeighPosition( m_vSites[ N3 ],  Site::WEST_UP );
    m_vSites[ N4 ]->setNeigh( m_vSites[ N6 ] );
    m_vSites[ N4 ]->setNeighPosition( m_vSites[ N6 ],  Site::WEST_DOWN );
    m_vSites[ N4 ]->setNeigh( m_vSites[ N4 - 1 ] );
    m_vSites[ N4 ]->setNeighPosition( m_vSites[ N4 - 1 ],  Site::EAST_UP );
    m_vSites[ N4 ]->setNeigh( m_vSites[ N1 ] );
    m_vSites[ N4 ]->setNeighPosition( m_vSites[ N1 ],  Site::EAST_DOWN );

    //N5
    m_vSites[ N5 ]->setNeigh( m_vSites[ N3 + 1 ] );
    m_vSites[ N5 ]->setNeighPosition( m_vSites[ N3 + 1 ],  Site::NORTH );
    m_vSites[ N5 ]->setNeigh( m_vSites[ N5 + m_iSizeX ] );
    m_vSites[ N5 ]->setNeighPosition( m_vSites[ N5 + m_iSizeX  ],  Site::SOUTH );
    m_vSites[ N5 ]->setNeigh( m_vSites[ N5 - 1 ] );
    m_vSites[ N5 ]->setNeighPosition( m_vSites[ N5 - 1 ],  Site::WEST_UP );
    m_vSites[ N5 ]->setNeigh( m_vSites[ N2 + 1 ] );
    m_vSites[ N5 ]->setNeighPosition( m_vSites[ N2 + 1 ],  Site::WEST_DOWN );
    m_vSites[ N5 ]->setNeigh( m_vSites[ N5 + 1 ] );
    m_vSites[ N5 ]->setNeighPosition( m_vSites[ N5 + 1 ],  Site::EAST_UP );
    m_vSites[ N5 ]->setNeigh( m_vSites[ N5 + m_iSizeX + 1] );
    m_vSites[ N5 ]->setNeighPosition( m_vSites[ N5 + m_iSizeX + 1 ],  Site::EAST_DOWN );

    //N6
    m_vSites[ N6 ]->setNeigh( m_vSites[ N7 ] );
    m_vSites[ N6 ]->setNeighPosition( m_vSites[ N7 ],  Site::NORTH );
    m_vSites[ N6 ]->setNeigh( m_vSites[ N6 + m_iSizeX ] );
    m_vSites[ N6 ]->setNeighPosition( m_vSites[ N6 + m_iSizeX  ],  Site::SOUTH );
    m_vSites[ N6 ]->setNeigh( m_vSites[ N7 - 1 ] );
    m_vSites[ N6 ]->setNeighPosition( m_vSites[ N7 - 1 ],  Site::WEST_UP );
    m_vSites[ N6 ]->setNeigh( m_vSites[ N6 - 1 ] );
    m_vSites[ N6 ]->setNeighPosition( m_vSites[ N6 - 1 ],  Site::WEST_DOWN );
    m_vSites[ N6 ]->setNeigh( m_vSites[ N4 ] );
    m_vSites[ N6 ]->setNeighPosition( m_vSites[ N4 ],  Site::EAST_UP );
    m_vSites[ N6 ]->setNeigh( m_vSites[ N2 ] );
    m_vSites[ N6 ]->setNeighPosition( m_vSites[ N2 ],  Site::EAST_DOWN );

    //N7
    m_vSites[ N7 ]->setNeigh( m_vSites[ N7 - m_iSizeX ] );
    m_vSites[ N7 ]->setNeighPosition( m_vSites[ N7 - m_iSizeX ],  Site::NORTH );
    m_vSites[ N7 ]->setNeigh( m_vSites[ N6 ] );
    m_vSites[ N7 ]->setNeighPosition( m_vSites[ N6  ],  Site::SOUTH );
    m_vSites[ N7 ]->setNeigh( m_vSites[ N7 - m_iSizeX - 1 ] );
    m_vSites[ N7 ]->setNeighPosition( m_vSites[ N7 - m_iSizeX - 1 ],  Site::WEST_UP );
    m_vSites[ N7 ]->setNeigh( m_vSites[ N7 - 1 ] );
    m_vSites[ N7 ]->setNeighPosition( m_vSites[ N7 - 1 ],  Site::WEST_DOWN );
    m_vSites[ N7 ]->setNeigh( m_vSites[ N7 - m_iSizeX + 1  ] );
    m_vSites[ N7 ]->setNeighPosition( m_vSites[ N7 - m_iSizeX + 1  ],  Site::EAST_UP );
    m_vSites[ N7 ]->setNeigh( m_vSites[ N7 + 1 ] );
    m_vSites[ N7 ]->setNeighPosition( m_vSites[ N7 + 1 ],  Site::EAST_DOWN );

    //N8
    m_vSites[ N8 ]->setNeigh( m_vSites[ N8 - m_iSizeX ] );
    m_vSites[ N8 ]->setNeighPosition( m_vSites[ N8 - m_iSizeX ],  Site::NORTH );
    m_vSites[ N8 ]->setNeigh( m_vSites[ N5 ] );
    m_vSites[ N8 ]->setNeighPosition( m_vSites[ N5  ],  Site::SOUTH );
    m_vSites[ N8 ]->setNeigh( m_vSites[ N8 - 1 ] );
    m_vSites[ N8 ]->setNeighPosition( m_vSites[ N8 - 1 ],  Site::WEST_UP );
    m_vSites[ N8 ]->setNeigh( m_vSites[ N1 ] );
    m_vSites[ N8 ]->setNeighPosition( m_vSites[ N1 ],  Site::WEST_DOWN );
    m_vSites[ N8 ]->setNeigh( m_vSites[ N8 + 1  ] );
    m_vSites[ N8 ]->setNeighPosition( m_vSites[ N8 + 1  ],  Site::EAST_UP );
    m_vSites[ N8 ]->setNeigh( m_vSites[ N5 + 1 ] );
    m_vSites[ N8 ]->setNeighPosition( m_vSites[ N5 + 1 ],  Site::EAST_DOWN );


    //first line
    for (int i = 2; i < m_iSizeX-2; i++){
        m_vSites[ i ]->setNeigh( m_vSites[ i + m_iSizeX ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeX ],  Site::SOUTH );
        m_vSites[ i ]->setNeigh( m_vSites[ N3 + i ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ N3 + i ],  Site::NORTH );

        if ( i%2 == 0){
            m_vSites[ i ]->setNeigh( m_vSites[ i - 1 ] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ i - 1 ],  Site::WEST_DOWN );
            m_vSites[ i ]->setNeigh( m_vSites[ i + 1 ] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ i + 1 ],  Site::EAST_DOWN );
            m_vSites[ i ]->setNeigh( m_vSites[ N3 + i - 1] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ N3 + i - 1 ],  Site::WEST_UP );
            m_vSites[ i ]->setNeigh( m_vSites[ N3 + i + 1] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ N3 + i + 1 ],  Site::EAST_UP );
        }
        else {
            m_vSites[ i ]->setNeigh( m_vSites[ i - 1 ] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ i - 1 ],  Site::WEST_UP );
            m_vSites[ i ]->setNeigh( m_vSites[ i + 1 ] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ i + 1 ],  Site::EAST_UP );
            m_vSites[ i ]->setNeigh( m_vSites[ i + m_iSizeX - 1 ] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeX - 1  ],  Site::WEST_DOWN );
            m_vSites[ i ]->setNeigh( m_vSites[ i + m_iSizeX + 1] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeX + 1 ],  Site::EAST_UP );
        }
    }

    //last line
    for (int i = (N8+1), iCount = 2; i < N7;  i++, iCount++){
        m_vSites[ i ]->setNeigh( m_vSites[ N1 + iCount ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ N1 + iCount  ],  Site::NORTH );
        m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeX ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeX ],  Site::SOUTH );

        if ( i%2 == 0){
            m_vSites[ i ]->setNeigh( m_vSites[ i - 1 ] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ i - 1 ],  Site::WEST_DOWN );
            m_vSites[ i ]->setNeigh( m_vSites[ i + 1 ] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ i + 1 ],  Site::EAST_DOWN );
            m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeX - 1 ] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeX - 1 ],  Site::WEST_UP );
            m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeX + 1 ] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeX + 1 ],  Site::EAST_UP );
        }
        else {
            m_vSites[ i ]->setNeigh( m_vSites[ i - 1 ] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ i - 1 ],  Site::WEST_UP );
            m_vSites[ i ]->setNeigh( m_vSites[ i + 1 ] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ i + 1 ],  Site::EAST_UP );
            m_vSites[ i ]->setNeigh( m_vSites[ iCount + 1 ] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ iCount + 1  ],  Site::WEST_DOWN );
            m_vSites[ i ]->setNeigh( m_vSites[ iCount - 1 ] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ iCount - 1 ],  Site::EAST_UP );
        }
    }

    //First column
    for (int i = ( N1 + m_iSizeX ); i < N3;  i += m_iSizeX ){
        m_vSites[ i ]->setNeigh( m_vSites[ i + m_iSizeX ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeX  ],  Site::SOUTH );
        m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeX ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeX ],  Site::NORTH );
        m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeX + 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeX + 1 ],  Site::EAST_UP );
        m_vSites[ i ]->setNeigh( m_vSites[ i + 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + 1 ],  Site::EAST_DOWN );
        m_vSites[ i ]->setNeigh( m_vSites[ i  - 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - 1 ],  Site::WEST_UP );
        m_vSites[ i ]->setNeigh( m_vSites[ i + m_iSizeX - 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeX - 1 ],  Site::WEST_DOWN );
    }

    //Second column
    for (int i = ( N5 + m_iSizeX ); i < N8;  i += m_iSizeX ){
        m_vSites[ i ]->setNeigh( m_vSites[ i + m_iSizeX ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeX  ],  Site::SOUTH );
        m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeX ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeX ],  Site::NORTH );
        m_vSites[ i ]->setNeigh( m_vSites[ i + 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + 1 ],  Site::EAST_UP );
        m_vSites[ i ]->setNeigh( m_vSites[ i + m_iSizeX + 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeX + 1 ],  Site::EAST_DOWN );
        m_vSites[ i ]->setNeigh( m_vSites[ i - 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - 1 ],  Site::WEST_UP );
        m_vSites[ i ]->setNeigh( m_vSites[ i + m_iSizeX - 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeX - 1 ],  Site::WEST_DOWN );
    }


    //Last column
    for ( int i = (N2 + m_iSizeX); i < (getSize() - m_iSizeX); i+=m_iSizeX ){
        m_vSites[ i ]->setNeigh( m_vSites[ i + m_iSizeX ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeX  ],  Site::SOUTH );
        m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeX ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeX ],  Site::NORTH );
        m_vSites[ i ]->setNeigh( m_vSites[ i - 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - 1 ],  Site::WEST_UP );
        m_vSites[ i ]->setNeigh( m_vSites[ i + m_iSizeX - 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeX - 1 ],  Site::WEST_DOWN );
        m_vSites[ i ]->setNeigh( m_vSites[ i + 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + 1 ],  Site::EAST_UP );
        m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeX + 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeX + 1 ],  Site::EAST_DOWN );
    }

    //Last column - 1
    for ( int i = (N6 + m_iSizeX); i < N7; i+=m_iSizeX ){
        m_vSites[ i ]->setNeigh( m_vSites[ i + m_iSizeX ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeX  ],  Site::SOUTH );
        m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeX ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeX ],  Site::NORTH );
        m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeX - 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeX - 1 ],  Site::WEST_UP );
        m_vSites[ i ]->setNeigh( m_vSites[ i - 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - 1 ],  Site::WEST_DOWN );
        m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeX + 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeX + 1 ],  Site::EAST_UP );
        m_vSites[ i ]->setNeigh( m_vSites[ i + 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + 1 ],  Site::EAST_DOWN );
    }

    //A very simple test that all neighs have been defined at least in number
    for ( int i = 0; i < getSize(); i++){
        if ( m_vSites[ i ]->getNeighs().size() != 6 ){
            cout << "Check neighs in site: " << i << " " << m_vSites[i]->getNeighs().size() << endl;
            EXIT
        }
    }

}

void HCP::check()
{
    cout << "Checking lattice..." << endl;

    int test = 2;
    cout << test << ": ";
    cout << "W:" << getSite(test)->getNeighPosition(Site::WEST)->getID() << " ";
    cout << "E:" << getSite(test)->getNeighPosition(Site::EAST)->getID() << " ";
    cout << "N:" << getSite(test)->getNeighPosition(Site::NORTH)->getID() << " ";
    cout << "S:" << getSite(test)->getNeighPosition(Site::SOUTH)->getID() << endl;
}

int HCP::calculateNeighNum(int id)
{
    int neighs = 1;
    for (Site *s : m_vSites[id]->getNeighs())
    {
        if (s->getHeight() >= m_vSites[id]->getHeight())
            neighs++;
    }

    // THIS IS BAD! REFACTOR ....
    m_vSites[id]->setNeighsNum(neighs);
    return neighs;
}

unordered_map<string, double> HCP::computeCoverages( vector<string> species)
{
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

