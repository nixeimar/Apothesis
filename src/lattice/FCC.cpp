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

#include "FCC.h"
#include "read.h"

FCC::FCC(Apothesis* apothesis):Lattice(apothesis)
{ 
    ;
}

void FCC::setInitialHeight( int  height ) { m_iHeight = height; }

void FCC::build()
{
    // The sites of the lattice.
    m_vSites.resize( getSize() );
    for ( int i = 0; i < m_vSites.size(); i++)
        m_vSites[ i ] = new Site();

    if ( m_sOrient == "111"){
        if ( m_iSizeX%2 != 0 || m_iSizeX%2 != 0){
            cout << "Error:The size of the lattice must be an even number in each direction." << endl;
            for ( int i = 0; i < m_vSites.size(); i++)
                delete m_vSites[ i ];
            EXIT;
        }

        if ( m_iSizeX == 0 || m_iSizeY == 0) {
            m_errorHandler->error_simple_msg("The lattice size cannot be zero in either dimension.");
            EXIT;
        }

        if ( m_iHeight < 5) {
            m_errorHandler->warningSimple_msg("The lattice initial height is too small.Consider revising.");
        }

        for (int i = 0; i < m_iSizeX; i++){
            if (i%2 == 0)
                for (int j = i*m_iSizeY; j < (m_iSizeY + i*m_iSizeY); j++)
                    if (j%2 == 0){
                        m_vSites[ j]->setID( j);
                        m_vSites[ j]->setHeight( m_iHeight -1 );
                    }
                    else{
                        m_vSites[ j]->setID( j);
                        m_vSites[ j]->setHeight( m_iHeight );
                    }
            else
                for (int j = i*m_iSizeY; j < (m_iSizeY + i*m_iSizeY); j++)
                    if (j%2 == 0){
                        m_vSites[ j]->setID( j);
                        m_vSites[ j]->setHeight( m_iHeight );
                    }
                    else{
                        m_vSites[ j]->setID( j);
                        m_vSites[ j]->setHeight( m_iHeight - 1);
                    }
        }

        mf_neigh_111();
    }
    else if ( m_sOrient == "110" ){

        if ( m_iSizeX%2 != 0 ){
            cout << "Error:The size of the lattice in the x-direciton must be an even number." << endl;
            for ( int i = 0; i < m_vSites.size(); i++)
                delete m_vSites[ i ];
            EXIT;
        }


        int iPos = 0;
        for (int i = 0; i< m_iSizeY; i++ ){
            for (int j = 0; j< m_iSizeX; j++ ){
                iPos = i*m_iSizeX + j;

                if ( j%2 != 0 ){
                    m_vSites[ iPos ]->setID( iPos );
                    m_vSites[ iPos ]->setHeight( m_iHeight - 1 );
                }
                else {
                    m_vSites[ iPos ]->setID( iPos );
                    m_vSites[ iPos ]->setHeight( m_iHeight );
                }
            }
        }

        mf_neigh_110();
    }
    else if ( m_sOrient == "100") {  }

}

FCC::~FCC()
{
    for ( int i = 0; i<getSize(); i++)
        delete m_vSites[i];
}

void FCC::mf_neigh_100()
{

}

void FCC::mf_neigh_110()
{

    //Nikos: FCC(110) must be an even number because there are differences in the level(-1,1) neighbors.

    //Example of a 6x6 lattice. Note that it the end row must always be a height below the frist one in order to have periodic conditions
    //The level (-1, 1) neighbors are different for the odd and even columns.

    //Example of a 6x6 lattice. The even columns are a height level above the odd units.

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

    // Example:
    // The neighbors for 14 in level 0 (same height) are: 12 (West), 16 (East), 8 (North), 20 (South)
    //             level -1, 1 (a level below/above) are: 7 (West Up), 13 (West down), 9 (East up), 15 (East down)

    // The neighbors for 15 in level 0 (same height) are: 13 (West), 17 (East), 9 (North), 21 (South)
    //             level -1, 1 (a level below/above) are: 14 (West Up), 20 (West down), 16 (East up), 22 (East down)

    int iPos = 0;
    /* All except the boundaries */
    for ( int i = 1; i < m_iSizeY - 1; i++ ){
        for (int j = 2; j < m_iSizeX -2; j++) {
            iPos  = i*m_iSizeX + j;

            //Same level
            m_vSites[ iPos ]->set1stNeibors(0, m_vSites[ iPos - m_iSizeX ] ); //North
            m_vSites[ iPos ]->setNeighPosition( m_vSites[ iPos - m_iSizeX ],  Site::NORTH );

            m_vSites[ iPos ]->set1stNeibors(0, m_vSites[ iPos + m_iSizeX ]); //South
            m_vSites[ iPos ]->setNeighPosition( m_vSites[ iPos + m_iSizeX ],  Site::SOUTH);

            m_vSites[ iPos ]->set1stNeibors(0, m_vSites[ iPos - 2 ]); //West
            m_vSites[ iPos ]->setNeighPosition( m_vSites[ iPos + m_iSizeX ],  Site::WEST);

            m_vSites[ iPos ]->set1stNeibors(0, m_vSites[ iPos + 2 ]); //East
            m_vSites[ iPos ]->setNeighPosition( m_vSites[ iPos + m_iSizeX ],  Site::WEST);

            //Next|Previous level
            if ( iPos%2 == 0){
                m_vSites[ iPos ]->set1stNeibors(-1, m_vSites[ iPos - m_iSizeX + 1 ] );
                m_vSites[ iPos ]->setNeighPosition( m_vSites[ iPos - m_iSizeX + 1 ],  Site::EAST_UP );

                m_vSites[ iPos ]->set1stNeibors(-1, m_vSites[ iPos - m_iSizeX - 1 ] );
                m_vSites[ iPos ]->setNeighPosition( m_vSites[ iPos - m_iSizeX - 1 ],  Site::WEST_UP );

                m_vSites[ iPos ]->set1stNeibors(-1, m_vSites[ iPos - 1 ] );
                m_vSites[ iPos ]->setNeighPosition( m_vSites[ iPos - 1 ],  Site::WEST_DOWN );
                m_vSites[ iPos ]->set1stNeibors(-1, m_vSites[ iPos + 1 ] );
                m_vSites[ iPos ]->setNeighPosition( m_vSites[ iPos + 1 ],  Site::EAST_DOWN );
            }
            else {
                m_vSites[ iPos ]->set1stNeibors(-1, m_vSites[ iPos + 1 ] );
                m_vSites[ iPos ]->setNeighPosition( m_vSites[ iPos +1 ],  Site::EAST_UP );
                m_vSites[ iPos ]->set1stNeibors(-1, m_vSites[ iPos - 1 ] );
                m_vSites[ iPos ]->setNeighPosition( m_vSites[ iPos - 1 ],  Site::WEST_UP );
                m_vSites[ iPos ]->set1stNeibors(-1, m_vSites[ iPos + m_iSizeX + 1 ] );
                m_vSites[ iPos ]->setNeighPosition( m_vSites[ iPos + m_iSizeX  + 1 ],  Site::EAST_DOWN );
                m_vSites[ iPos ]->set1stNeibors(-1, m_vSites[ iPos + m_iSizeX  - 1 ] );
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
    m_vSites[ N1 ]->set1stNeibors(0, m_vSites[ 2 ] );
    m_vSites[ N1 ]->setNeighPosition( m_vSites[ 2 ],  Site::EAST );
    m_vSites[ N1 ]->set1stNeibors(0, m_vSites[ N2 - 1 ] );
    m_vSites[ N1 ]->setNeighPosition( m_vSites[ N2 - 1 ],  Site::WEST );
    m_vSites[ N1 ]->set1stNeibors(0, m_vSites[ m_iSizeX ] );
    m_vSites[ N1 ]->setNeighPosition( m_vSites[ iPos + m_iSizeX ],  Site::SOUTH );
    m_vSites[ N1 ]->set1stNeibors(0, m_vSites[ getSize() - m_iSizeX ] );
    m_vSites[ N1 ]->setNeighPosition( m_vSites[ iPos + m_iSizeX ],  Site::NORTH );
    m_vSites[ N1 ]->set1stNeibors( -1, m_vSites[ N1 + 1 ] );
    m_vSites[ N1 ]->setNeighPosition( m_vSites[ N1 + 1 ],  Site::EAST_DOWN );
    m_vSites[ N1 ]->set1stNeibors( -1, m_vSites[ N3 + 1 ] );
    m_vSites[ N1 ]->setNeighPosition( m_vSites[ N3 + 1 ],  Site::EAST_UP );
    m_vSites[ N1 ]->set1stNeibors( -1, m_vSites[ N2 ] );
    m_vSites[ N1 ]->setNeighPosition( m_vSites[ N2 ],  Site::WEST_DOWN );
    m_vSites[ N1 ]->set1stNeibors( -1, m_vSites[ N4 ] );
    m_vSites[ N1 ]->setNeighPosition( m_vSites[ N4 ],  Site::WEST_UP );

    //N2
    m_vSites[ N2 ]->set1stNeibors(0, m_vSites[ N2 - 2 ] );
    m_vSites[ N2 ]->setNeighPosition( m_vSites[ N2 - 2 ],  Site::WEST );
    m_vSites[ N2 ]->set1stNeibors(0, m_vSites[ N1 + 1 ] );
    m_vSites[ N2 ]->setNeighPosition( m_vSites[ N1 + 1 ],  Site::EAST );
    m_vSites[ N2 ]->set1stNeibors(0, m_vSites[ N4 ] );
    m_vSites[ N2 ]->setNeighPosition( m_vSites[ N4 ],  Site::NORTH );
    m_vSites[ N2 ]->set1stNeibors(0, m_vSites[ N2 + m_iSizeX ] );
    m_vSites[ N2 ]->setNeighPosition( m_vSites[ N2 + m_iSizeX ],  Site::SOUTH );
    m_vSites[ N2 ]->set1stNeibors( -1, m_vSites[ N2 - 1 ] );
    m_vSites[ N2 ]->setNeighPosition( m_vSites[ N2 - 1 ],  Site::WEST_UP );
    m_vSites[ N2 ]->set1stNeibors( -1, m_vSites[ N2 + m_iSizeX - 1 ] );
    m_vSites[ N2 ]->setNeighPosition( m_vSites[ N2 + m_iSizeX - 1 ],  Site::WEST_DOWN );
    m_vSites[ N2 ]->set1stNeibors( -1, m_vSites[ m_iSizeX ] );
    m_vSites[ N2 ]->setNeighPosition( m_vSites[ m_iSizeX ],  Site::EAST_DOWN );
    m_vSites[ N2 ]->set1stNeibors( -1, m_vSites[ N1 ] );
    m_vSites[ N2 ]->setNeighPosition( m_vSites[ N1 ],  Site::EAST_UP );

    //N3
    m_vSites[ N3 ]->set1stNeibors(0, m_vSites[ N3 + 2 ] );
    m_vSites[ N3 ]->setNeighPosition( m_vSites[ N3 + 2 ],  Site::EAST );
    m_vSites[ N3 ]->set1stNeibors(0, m_vSites[ N4 - 1 ] );
    m_vSites[ N3 ]->setNeighPosition( m_vSites[ N4 - 1 ],  Site::WEST );
    m_vSites[ N3 ]->set1stNeibors(0, m_vSites[ N1 ] );
    m_vSites[ N3 ]->setNeighPosition( m_vSites[ N1 ],  Site::SOUTH );
    m_vSites[ N3 ]->set1stNeibors(0, m_vSites[ N3 - m_iSizeX ] );
    m_vSites[ N3 ]->setNeighPosition( m_vSites[ N3 - m_iSizeX ],  Site::NORTH );
    m_vSites[ N3 ]->set1stNeibors( -1, m_vSites[ N3 - m_iSizeX + 1 ] );
    m_vSites[ N3 ]->setNeighPosition( m_vSites[ N3 - m_iSizeX + 1 ],  Site::EAST_UP );
    m_vSites[ N3 ]->set1stNeibors( -1, m_vSites[ N3 + 1 ] );
    m_vSites[ N3 ]->setNeighPosition( m_vSites[ N3 + 1 ],  Site::EAST_DOWN );
    m_vSites[ N3 ]->set1stNeibors( -1, m_vSites[ N4 ] );
    m_vSites[ N3 ]->setNeighPosition( m_vSites[ N4 ],  Site::WEST_DOWN );
    m_vSites[ N3 ]->set1stNeibors( -1, m_vSites[ N3 - 1 ] );
    m_vSites[ N3 ]->setNeighPosition( m_vSites[ N3 - 1 ],  Site::WEST_UP );

    //N4
    m_vSites[ N4 ]->set1stNeibors(0, m_vSites[ N3 + 1 ] );
    m_vSites[ N4 ]->setNeighPosition( m_vSites[ N3 + 1 ],  Site::EAST );
    m_vSites[ N4 ]->set1stNeibors(0, m_vSites[ N4 - 2 ] );
    m_vSites[ N4 ]->setNeighPosition( m_vSites[ N4 - 2 ],  Site::WEST );
    m_vSites[ N4 ]->set1stNeibors(0, m_vSites[ N2 ] );
    m_vSites[ N4 ]->setNeighPosition( m_vSites[ N2 ],  Site::SOUTH );
    m_vSites[ N4 ]->set1stNeibors(0, m_vSites[ N4 - m_iSizeX ] );
    m_vSites[ N4 ]->setNeighPosition( m_vSites[ N4 - m_iSizeX  ],  Site::NORTH );
    m_vSites[ N4 ]->set1stNeibors( -1, m_vSites[ N3 ] );
    m_vSites[ N4 ]->setNeighPosition( m_vSites[ N3 ],  Site::WEST_UP );
    m_vSites[ N4 ]->set1stNeibors( -1, m_vSites[ N6 ] );
    m_vSites[ N4 ]->setNeighPosition( m_vSites[ N6 ],  Site::WEST_DOWN );
    m_vSites[ N4 ]->set1stNeibors( -1, m_vSites[ N4 - 1 ] );
    m_vSites[ N4 ]->setNeighPosition( m_vSites[ N4 - 1 ],  Site::EAST_UP );
    m_vSites[ N4 ]->set1stNeibors( -1, m_vSites[ N1 ] );
    m_vSites[ N4 ]->setNeighPosition( m_vSites[ N1 ],  Site::EAST_DOWN );

    //N5
    m_vSites[ N5 ]->set1stNeibors(0, m_vSites[ N5 + 2 ] );
    m_vSites[ N5 ]->setNeighPosition( m_vSites[ N5 + 2 ],  Site::EAST );
    m_vSites[ N5 ]->set1stNeibors(0, m_vSites[ N4 - 2 ] );
    m_vSites[ N5 ]->setNeighPosition( m_vSites[ N4 - 2 ],  Site::WEST );
    m_vSites[ N5 ]->set1stNeibors(0, m_vSites[ N3 + 1 ] );
    m_vSites[ N5 ]->setNeighPosition( m_vSites[ N3 + 1 ],  Site::NORTH );
    m_vSites[ N5 ]->set1stNeibors(0, m_vSites[ N5 + m_iSizeX ] );
    m_vSites[ N5 ]->setNeighPosition( m_vSites[ N5 + m_iSizeX  ],  Site::SOUTH );
    m_vSites[ N5 ]->set1stNeibors( -1, m_vSites[ N5 - 1 ] );
    m_vSites[ N5 ]->setNeighPosition( m_vSites[ N5 - 1 ],  Site::WEST_UP );
    m_vSites[ N5 ]->set1stNeibors( -1, m_vSites[ N2 + 1 ] );
    m_vSites[ N5 ]->setNeighPosition( m_vSites[ N2 + 1 ],  Site::WEST_DOWN );
    m_vSites[ N5 ]->set1stNeibors( -1, m_vSites[ N5 + 1 ] );
    m_vSites[ N5 ]->setNeighPosition( m_vSites[ N5 + 1 ],  Site::EAST_UP );
    m_vSites[ N5 ]->set1stNeibors( -1, m_vSites[ N5 + m_iSizeX + 1] );
    m_vSites[ N5 ]->setNeighPosition( m_vSites[ N5 + m_iSizeX + 1 ],  Site::EAST_DOWN );

    //N6
    m_vSites[ N6 ]->set1stNeibors(0, m_vSites[ N1 ] );
    m_vSites[ N6 ]->setNeighPosition( m_vSites[ N1 ],  Site::EAST );
    m_vSites[ N6 ]->set1stNeibors(0, m_vSites[ N6 - 2 ] );
    m_vSites[ N6 ]->setNeighPosition( m_vSites[ N6 - 2 ],  Site::WEST );
    m_vSites[ N6 ]->set1stNeibors(0, m_vSites[ N7 ] );
    m_vSites[ N6 ]->setNeighPosition( m_vSites[ N7 ],  Site::NORTH );
    m_vSites[ N6 ]->set1stNeibors(0, m_vSites[ N6 + m_iSizeX ] );
    m_vSites[ N6 ]->setNeighPosition( m_vSites[ N6 + m_iSizeX  ],  Site::SOUTH );
    m_vSites[ N6 ]->set1stNeibors( -1, m_vSites[ N7 - 1 ] );
    m_vSites[ N6 ]->setNeighPosition( m_vSites[ N7 - 1 ],  Site::WEST_UP );
    m_vSites[ N6 ]->set1stNeibors( -1, m_vSites[ N6 - 1 ] );
    m_vSites[ N6 ]->setNeighPosition( m_vSites[ N6 - 1 ],  Site::WEST_DOWN );
    m_vSites[ N6 ]->set1stNeibors( -1, m_vSites[ N4 ] );
    m_vSites[ N6 ]->setNeighPosition( m_vSites[ N4 ],  Site::EAST_UP );
    m_vSites[ N6 ]->set1stNeibors( -1, m_vSites[ N2 ] );
    m_vSites[ N6 ]->setNeighPosition( m_vSites[ N2 ],  Site::EAST_DOWN );

    //N7
    m_vSites[ N7 ]->set1stNeibors(0, m_vSites[ N3 ] );
    m_vSites[ N7 ]->setNeighPosition( m_vSites[ N3 ],  Site::EAST );
    m_vSites[ N7 ]->set1stNeibors(0, m_vSites[ N7 - 2 ] );
    m_vSites[ N7 ]->setNeighPosition( m_vSites[ N7 - 2 ],  Site::WEST );
    m_vSites[ N7 ]->set1stNeibors(0, m_vSites[ N7 - m_iSizeX ] );
    m_vSites[ N7 ]->setNeighPosition( m_vSites[ N7 - m_iSizeX ],  Site::NORTH );
    m_vSites[ N7 ]->set1stNeibors(0, m_vSites[ N6 ] );
    m_vSites[ N7 ]->setNeighPosition( m_vSites[ N6  ],  Site::SOUTH );
    m_vSites[ N7 ]->set1stNeibors( -1, m_vSites[ N7 - m_iSizeX - 1 ] );
    m_vSites[ N7 ]->setNeighPosition( m_vSites[ N7 - m_iSizeX - 1 ],  Site::WEST_UP );
    m_vSites[ N7 ]->set1stNeibors( -1, m_vSites[ N7 - 1 ] );
    m_vSites[ N7 ]->setNeighPosition( m_vSites[ N7 - 1 ],  Site::WEST_DOWN );
    m_vSites[ N7 ]->set1stNeibors( -1, m_vSites[ N7 - m_iSizeX + 1  ] );
    m_vSites[ N7 ]->setNeighPosition( m_vSites[ N7 - m_iSizeX + 1  ],  Site::EAST_UP );
    m_vSites[ N7 ]->set1stNeibors( -1, m_vSites[ N7 + 1 ] );
    m_vSites[ N7 ]->setNeighPosition( m_vSites[ N7 + 1 ],  Site::EAST_DOWN );

    //N8
    m_vSites[ N8 ]->set1stNeibors(0, m_vSites[ N8 + 2 ] );
    m_vSites[ N8 ]->setNeighPosition( m_vSites[ N8 + 2 ],  Site::EAST );
    m_vSites[ N8 ]->set1stNeibors(0, m_vSites[ N4 ] );
    m_vSites[ N8 ]->setNeighPosition( m_vSites[ N4 ],  Site::WEST );
    m_vSites[ N8 ]->set1stNeibors(0, m_vSites[ N8 - m_iSizeX ] );
    m_vSites[ N8 ]->setNeighPosition( m_vSites[ N8 - m_iSizeX ],  Site::NORTH );
    m_vSites[ N8 ]->set1stNeibors(0, m_vSites[ N5 ] );
    m_vSites[ N8 ]->setNeighPosition( m_vSites[ N5  ],  Site::SOUTH );
    m_vSites[ N8 ]->set1stNeibors( -1, m_vSites[ N8 - 1 ] );
    m_vSites[ N8 ]->setNeighPosition( m_vSites[ N8 - 1 ],  Site::WEST_UP );
    m_vSites[ N8 ]->set1stNeibors( -1, m_vSites[ N1 ] );
    m_vSites[ N8 ]->setNeighPosition( m_vSites[ N1 ],  Site::WEST_DOWN );
    m_vSites[ N8 ]->set1stNeibors( -1, m_vSites[ N8 + 1  ] );
    m_vSites[ N8 ]->setNeighPosition( m_vSites[ N8 + 1  ],  Site::EAST_UP );
    m_vSites[ N8 ]->set1stNeibors( -1, m_vSites[ N5 + 1 ] );
    m_vSites[ N8 ]->setNeighPosition( m_vSites[ N5 + 1 ],  Site::EAST_DOWN );


    //first line
    for (int i = 2; i < m_iSizeX-2; i++){
        m_vSites[ i ]->set1stNeibors( 0, m_vSites[ i - 2 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - 2 ],  Site::WEST );
        m_vSites[ i ]->set1stNeibors( 0, m_vSites[ i + 2 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + 2 ],  Site::EAST );
        m_vSites[ i ]->set1stNeibors( 0, m_vSites[ i + m_iSizeX ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeX ],  Site::SOUTH );
        m_vSites[ i ]->set1stNeibors( 0, m_vSites[ N3 + i ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ N3 + i ],  Site::NORTH );

        if ( i%2 == 0){
            m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i - 1 ] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ i - 1 ],  Site::WEST_DOWN );
            m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i + 1 ] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ i + 1 ],  Site::EAST_DOWN );
            m_vSites[ i ]->set1stNeibors( -1, m_vSites[ N3 + i - 1] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ N3 + i - 1 ],  Site::WEST_UP );
            m_vSites[ i ]->set1stNeibors( -1, m_vSites[ N3 + i + 1] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ N3 + i + 1 ],  Site::EAST_UP );
        }
        else {
            m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i - 1 ] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ i - 1 ],  Site::WEST_UP );
            m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i + 1 ] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ i + 1 ],  Site::EAST_UP );
            m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i + m_iSizeX - 1 ] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeX - 1  ],  Site::WEST_DOWN );
            m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i + m_iSizeX + 1] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeX + 1 ],  Site::EAST_UP );
        }
    }

    //last line
    for (int i = (N8+1), iCount = 2; i < N7;  i++, iCount++){
        m_vSites[ i ]->set1stNeibors( 0, m_vSites[ i - 2 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - 2 ],  Site::WEST );
        m_vSites[ i ]->set1stNeibors( 0, m_vSites[ i + 2 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + 2 ],  Site::EAST );
        m_vSites[ i ]->set1stNeibors( 0, m_vSites[ N1 + iCount ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ N1 + iCount  ],  Site::SOUTH );
        m_vSites[ i ]->set1stNeibors( 0, m_vSites[ i - m_iSizeX ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeX ],  Site::NORTH );

        if ( i%2 == 0){
            m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i - 1 ] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ i - 1 ],  Site::WEST_DOWN );
            m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i + 1 ] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ i + 1 ],  Site::EAST_DOWN );
            m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i - m_iSizeX - 1 ] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeX - 1 ],  Site::WEST_UP );
            m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i - m_iSizeX + 1 ] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeX + 1 ],  Site::EAST_UP );
        }
        else {
            m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i - 1 ] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ i - 1 ],  Site::WEST_UP );
            m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i + 1 ] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ i + 1 ],  Site::EAST_UP );
            m_vSites[ i ]->set1stNeibors( -1, m_vSites[ iCount + 1 ] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ iCount + 1  ],  Site::WEST_DOWN );
            m_vSites[ i ]->set1stNeibors( -1, m_vSites[ iCount - 1 ] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ iCount - 1 ],  Site::EAST_UP );
        }
    }


    //First column
    for (int i = ( N1 + m_iSizeX ); i < N3;  i += m_iSizeX ){
        m_vSites[ i ]->set1stNeibors( 0, m_vSites[ i + m_iSizeX - 2 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeX - 2 ],  Site::WEST );
        m_vSites[ i ]->set1stNeibors( 0, m_vSites[ i + 2 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + 2 ],  Site::EAST );
        m_vSites[ i ]->set1stNeibors( 0, m_vSites[ i + m_iSizeX ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeX  ],  Site::SOUTH );
        m_vSites[ i ]->set1stNeibors( 0, m_vSites[ i - m_iSizeX ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeX ],  Site::NORTH );
        m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i - m_iSizeX + 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeX + 1 ],  Site::EAST_UP );
        m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i + 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + 1 ],  Site::EAST_DOWN );
        m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i  - 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - 1 ],  Site::WEST_UP );
        m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i + m_iSizeX - 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeX - 1 ],  Site::WEST_DOWN );
    }

    //Second column
    for (int i = ( N5 + m_iSizeX ); i < N8;  i += m_iSizeX ){
        m_vSites[ i ]->set1stNeibors( 0, m_vSites[ i + m_iSizeX  - 2 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeX  - 2 ],  Site::WEST );
        m_vSites[ i ]->set1stNeibors( 0, m_vSites[ i + 2 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + 2 ],  Site::EAST );
        m_vSites[ i ]->set1stNeibors( 0, m_vSites[ i + m_iSizeX ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeX  ],  Site::SOUTH );
        m_vSites[ i ]->set1stNeibors( 0, m_vSites[ i - m_iSizeX ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeX ],  Site::NORTH );
        m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i + 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + 1 ],  Site::EAST_UP );
        m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i + m_iSizeX + 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeX + 1 ],  Site::EAST_DOWN );
        m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i - 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - 1 ],  Site::WEST_UP );
        m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i + m_iSizeX - 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeX - 1 ],  Site::WEST_DOWN );
    }


    //Last column
    for ( int i = (N2 + m_iSizeX); i < (getSize() - m_iSizeX); i+=m_iSizeX ){
        m_vSites[ i ]->set1stNeibors( 0, m_vSites[ i - 2 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - 2 ],  Site::WEST );
        m_vSites[ i ]->set1stNeibors( 0, m_vSites[ i - m_iSizeX + 2 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeX + 2 ],  Site::EAST );
        m_vSites[ i ]->set1stNeibors( 0, m_vSites[ i + m_iSizeX ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeX  ],  Site::SOUTH );
        m_vSites[ i ]->set1stNeibors( 0, m_vSites[ i - m_iSizeX ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeX ],  Site::NORTH );
        m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i - 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - 1 ],  Site::WEST_UP );
        m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i + m_iSizeX - 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeX - 1 ],  Site::WEST_DOWN );
        m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i + 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + 1 ],  Site::EAST_UP );
        m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i - m_iSizeX + 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeX + 1 ],  Site::EAST_DOWN );
     }

    //Last column - 1
    for ( int i = (N6 + m_iSizeX); i < N7; i+=m_iSizeX ){
        m_vSites[ i ]->set1stNeibors( 0, m_vSites[ i - 2 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - 2 ],  Site::WEST );
        m_vSites[ i ]->set1stNeibors( 0, m_vSites[ i - m_iSizeX + 2 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeX + 2 ],  Site::EAST );
        m_vSites[ i ]->set1stNeibors( 0, m_vSites[ i + m_iSizeX ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeX  ],  Site::SOUTH );
        m_vSites[ i ]->set1stNeibors( 0, m_vSites[ i - m_iSizeX ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeX ],  Site::NORTH );
        m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i - m_iSizeX - 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeX - 1 ],  Site::WEST_UP );
        m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i - 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - 1 ],  Site::WEST_DOWN );
        m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i - m_iSizeX + 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeX + 1 ],  Site::EAST_UP );
        m_vSites[ i ]->set1stNeibors( -1, m_vSites[ i + 1 ] );
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + 1 ],  Site::EAST_DOWN );
    }

    //A very simple test that all neighs have been defined at least in number
    for ( int i = 0; i < getSize(); i++){
        if ( m_vSites[ i ]->get1stNeihbors()[ 0 ].size() != 4 ){
            cout << "Check level 0 neighs in site: " << i << endl;
            EXIT;
        }

        if ( m_vSites[ i ]->get1stNeihbors()[ 0 ].size() != 4 ){
            cout << "Check level 0 neighs in site: " << i << endl;
            EXIT;
        }
    }
}

void FCC::mf_neigh_111()
{
    /* All except the boundaries */
    for ( int i = 2; i < m_iSizeX - 2; i++ ){
        for (int j = (i*m_iSizeY + 2); j < (i*m_iSizeY + m_iSizeY - 2); j++ ){
            m_vSites[ j]->setNeigh( m_vSites[ j - 2]);
            m_vSites[ j]->setNeighPosition( m_vSites[ j - 2], Site::EAST);

            m_vSites[ j]->setNeigh( m_vSites[ j + 2]);
            m_vSites[ j]->setNeighPosition( m_vSites[ j + 2 ], Site::WEST);

            m_vSites[ j]->setNeigh( m_vSites[ j + m_iSizeY - 1]);
            m_vSites[ j]->setNeighPosition( m_vSites[ j + m_iSizeY - 1 ], Site::EAST_DOWN);

            m_vSites[ j]->setNeigh( m_vSites[ j + m_iSizeY + 1 ]);
            m_vSites[ j]->setNeighPosition( m_vSites[  j + m_iSizeY + 1 ], Site::WEST_DOWN);

            m_vSites[ j]->setNeigh( m_vSites[ j - m_iSizeY - 1]);
            m_vSites[ j]->setNeighPosition( m_vSites[ j - m_iSizeY - 1 ], Site::EAST_UP);

            m_vSites[ j]->setNeigh( m_vSites[ j - m_iSizeY + 1]);
            m_vSites[ j]->setNeighPosition( m_vSites[ j - m_iSizeY + 1 ], Site::WEST_UP);

            m_vSites[ j]->storeActivationSite( m_vSites[ j - 1 ], Site::ACTV_EAST );
            m_vSites[ j]->storeActivationSite( m_vSites[ j + 1], Site::ACTV_WEST );

            m_vSites[ j]->storeActivationSite( m_vSites[ j + m_iSizeY ], Site::ACTV_SOUTH );
            m_vSites[ j]->storeActivationSite( m_vSites[ j - m_iSizeY ], Site::ACTV_NORTH );

            m_vSites[ j]->setNeighPosition( m_vSites[ j - 2*m_iSizeY  ], Site::NORTH);
            m_vSites[ j]->setNeighPosition( m_vSites[ j + 2*m_iSizeY ], Site::SOUTH);
        }
    }

    /* First row */
    for ( int i = 0; i < m_iSizeY; i++){
        /* First site */
        if ( i == 0) {
            m_vSites[ i ]->setNeigh( m_vSites[ 2]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ 2], Site::EAST);

            m_vSites[ i ]->setNeigh( m_vSites[ m_iSizeY-2]);
            m_vSites[ i ]->setNeighPosition( m_vSites[ m_iSizeY-2], Site::WEST);

            m_vSites[ i ]->setNeigh( m_vSites[ m_iSizeY+1]);
            m_vSites[ i ]->setNeighPosition( m_vSites[ m_iSizeY+1], Site::EAST_DOWN);

            m_vSites[ i ]->setNeigh( m_vSites[ 2*m_iSizeY-1]);
            m_vSites[ i ]->setNeighPosition( m_vSites[  2*m_iSizeY-1], Site::WEST_DOWN);

            m_vSites[ i ]->setNeigh( m_vSites[ getSize()-m_iSizeY+1]);
            m_vSites[ i ]->setNeighPosition( m_vSites[  getSize()-m_iSizeY+1], Site::EAST_UP);

            m_vSites[ i ]->setNeigh( m_vSites[ getSize()-1]);
            m_vSites[ i ]->setNeighPosition( m_vSites[  getSize()-1], Site::WEST_UP);

            m_vSites[ i ]->setNeighPosition( m_vSites[  getSize() - 2*m_iSizeY], Site::NORTH);
            m_vSites[ i ]->setNeighPosition( m_vSites[  2*m_iSizeY ], Site::SOUTH);

            m_vSites[ i ]->storeActivationSite( m_vSites[ m_iSizeY - 1], Site::ACTV_EAST );
            m_vSites[ i ]->storeActivationSite( m_vSites[ 1], Site::ACTV_WEST );
            m_vSites[ i ]->storeActivationSite( m_vSites[ m_iSizeY ], Site::ACTV_SOUTH );
            m_vSites[ i ]->storeActivationSite( m_vSites[ getSize() - m_iSizeY ], Site::ACTV_NORTH );
        }
        /* Second site */
        else if ( i == 1){
            m_vSites[ i ]->setNeigh( m_vSites[ 3]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ 3], Site::WEST);

            m_vSites[ i ]->setNeigh( m_vSites[ m_iSizeY - 1]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ m_iSizeY - 1], Site::EAST);

            m_vSites[ i ]->setNeigh( m_vSites[ m_iSizeY]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ m_iSizeY], Site::EAST_DOWN);

            m_vSites[ i ]->setNeigh( m_vSites[ m_iSizeY + 2]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ m_iSizeY + 2], Site::WEST_DOWN);

            m_vSites[ i ]->setNeigh( m_vSites[ getSize() - m_iSizeY]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ getSize() - m_iSizeY], Site::EAST_UP);

            m_vSites[ i ]->setNeigh( m_vSites[ getSize() -m_iSizeY + 2]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ getSize() - m_iSizeY + 2], Site::WEST_UP);

            m_vSites[ i ]->setNeighPosition(m_vSites[ getSize() - 2*m_iSizeY + 1], Site::NORTH);
            m_vSites[ i ]->setNeighPosition(m_vSites[ 2*m_iSizeY + 1 ], Site::SOUTH);

            m_vSites[ i ]->storeActivationSite( m_vSites[ 0], Site::ACTV_EAST );
            m_vSites[ i ]->storeActivationSite( m_vSites[ 2], Site::ACTV_WEST );
            m_vSites[ i ]->storeActivationSite( m_vSites[ m_iSizeY + 1 ], Site::ACTV_SOUTH );
            m_vSites[ i ]->storeActivationSite( m_vSites[ getSize() - m_iSizeY + 1 ], Site::ACTV_NORTH );
        }
        /* Before last site */
        else if (i == m_iSizeY - 2){
            m_vSites[ i ]->setNeigh( m_vSites[ 0]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ 0], Site::WEST);

            m_vSites[ i ]->setNeigh( m_vSites[ m_iSizeY-4]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ m_iSizeY-4], Site::EAST);

            m_vSites[ i ]->setNeigh( m_vSites[ 2*m_iSizeY -3]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ 2*m_iSizeY -3], Site::EAST_DOWN);

            m_vSites[ i ]->setNeigh( m_vSites[ 2*m_iSizeY -1]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ 2*m_iSizeY -1], Site::WEST_DOWN);

            m_vSites[ i ]->setNeigh( m_vSites[ getSize() -1]);
            m_vSites[ i ]->setNeighPosition(m_vSites[  getSize() -1], Site::WEST_UP);

            m_vSites[ i ]->setNeigh( m_vSites[ getSize() -3]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ getSize() -3], Site::EAST_UP);

            m_vSites[ i ]->setNeighPosition(m_vSites[ getSize() - m_iSizeY - 2], Site::NORTH);
            m_vSites[ i ]->setNeighPosition(m_vSites[ 3*m_iSizeY - 2], Site::SOUTH);

            m_vSites[ i ]->storeActivationSite( m_vSites[ i - 1], Site::ACTV_EAST );
            m_vSites[ i ]->storeActivationSite( m_vSites[ i + 1], Site::ACTV_WEST );
            m_vSites[ i ]->storeActivationSite( m_vSites[ m_iSizeY + i ], Site::ACTV_SOUTH );
            m_vSites[ i ]->storeActivationSite( m_vSites[ getSize() - 2 ], Site::ACTV_NORTH );
        }
        /* End site */
        else if( i == m_iSizeY-1){
            m_vSites[ i ]->setNeigh( m_vSites[ 1]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ 1], Site::WEST);

            m_vSites[ i ]->setNeigh( m_vSites[ m_iSizeY-3]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ m_iSizeY-3], Site::EAST);

            m_vSites[ i ]->setNeigh( m_vSites[ m_iSizeY]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ m_iSizeY], Site::WEST_DOWN);

            m_vSites[ i ]->setNeigh( m_vSites[ 2*m_iSizeY - 2]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ 2*m_iSizeY - 2], Site::EAST_DOWN);

            m_vSites[ i ]->setNeigh( m_vSites[ getSize()-m_iSizeY]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ getSize()-m_iSizeY], Site::WEST_UP);

            m_vSites[ i ]->setNeigh( m_vSites[ getSize()-2]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ getSize()-2], Site::EAST_UP);

            m_vSites[ i ]->setNeighPosition(m_vSites[ getSize() - m_iSizeY - 1], Site::NORTH);
            m_vSites[ i ]->setNeighPosition(m_vSites[ 3*m_iSizeY - 1], Site::SOUTH);

            m_vSites[ i ]->storeActivationSite( m_vSites[ i - 1], Site::ACTV_EAST );
            m_vSites[ i ]->storeActivationSite( m_vSites[ 0], Site::ACTV_WEST );
            m_vSites[ i ]->storeActivationSite( m_vSites[ m_iSizeY + i ], Site::ACTV_SOUTH );
            m_vSites[ i ]->storeActivationSite( m_vSites[ getSize() - 1 ], Site::ACTV_NORTH );
        }
        else {
            m_vSites[ i ]->setNeigh( m_vSites[ i+m_iSizeY-1]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ i + m_iSizeY - 1], Site::EAST_DOWN);

            m_vSites[ i ]->setNeigh( m_vSites[ i+m_iSizeY+1]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ i+m_iSizeY+1], Site::WEST_DOWN);

            m_vSites[ i ]->setNeigh( m_vSites[ i-2]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ i-2], Site::EAST);

            m_vSites[ i ]->setNeigh( m_vSites[ i+2]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ i+2], Site::WEST);

            m_vSites[ i ]->setNeigh( m_vSites[ i+getSize()-m_iSizeY-1]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ i+getSize()-m_iSizeY-1], Site::EAST_UP);

            m_vSites[ i ]->setNeigh( m_vSites[ i+getSize()-m_iSizeY+1]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ i+getSize()-m_iSizeY+1], Site::WEST_UP);

            m_vSites[ i ]->setNeighPosition(m_vSites[ getSize() - 2*m_iSizeY + i], Site::NORTH);
            m_vSites[ i ]->setNeighPosition(m_vSites[ i + 2*m_iSizeY ], Site::SOUTH);

            m_vSites[ i ]->storeActivationSite( m_vSites[ i - 1], Site::ACTV_EAST );
            m_vSites[ i ]->storeActivationSite( m_vSites[ i + 1], Site::ACTV_WEST );
            m_vSites[ i ]->storeActivationSite( m_vSites[ m_iSizeY + i ], Site::ACTV_SOUTH );
            m_vSites[ i ]->storeActivationSite( m_vSites[ getSize() - m_iSizeY + i ], Site::ACTV_NORTH );
        }
    }

    /* Last row */
    for ( int i = (getSize() - m_iSizeY); i < getSize(); i++ ){
        if ( i == ( getSize() - m_iSizeY ) ){ //bottom left site
            m_vSites[ i ]->setNeigh( m_vSites[ m_iSizeY - 1] );
            m_vSites[ i ]->setNeighPosition( m_vSites[ m_iSizeY - 1 ], Site::EAST_DOWN);

            m_vSites[ i ]->setNeigh( m_vSites[ 1]);
            m_vSites[ i ]->setNeighPosition( m_vSites[ 1], Site::WEST_DOWN);

            m_vSites[ i ]->setNeigh( m_vSites[ i + 2]);
            m_vSites[ i ]->setNeighPosition( m_vSites[ i + 2], Site::WEST);

            m_vSites[ i ]->setNeigh( m_vSites[ getSize() - 2 ]);
            m_vSites[ i ]->setNeighPosition( m_vSites[ getSize() - 2 ], Site::EAST);

            m_vSites[ i ]->setNeigh( m_vSites[ i - 1]);
            m_vSites[ i ]->setNeighPosition( m_vSites[ i - 1], Site::EAST_UP);

            m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeY + 1]);
            m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeY + 1], Site::WEST_UP);

            m_vSites[ i ]->setNeighPosition( m_vSites[ i - 2*m_iSizeY ], Site::NORTH );
            m_vSites[ i ]->setNeighPosition( m_vSites[ m_iSizeY ], Site::SOUTH );

            m_vSites[ i ]->storeActivationSite( m_vSites[ getSize() - 1 ], Site::ACTV_EAST );
            m_vSites[ i ]->storeActivationSite( m_vSites[ i + 1], Site::ACTV_WEST );
            m_vSites[ i ]->storeActivationSite( m_vSites[ 0 ], Site::ACTV_SOUTH );
            m_vSites[ i ]->storeActivationSite( m_vSites[ i - m_iSizeY ], Site::ACTV_NORTH );
        }
        else if ( i == ( getSize() - 1 ) ){ //bottom right site
            m_vSites[ i ]->setNeigh( m_vSites[  m_iSizeY - 2 ] );
            m_vSites[ i ]->setNeighPosition( m_vSites[  m_iSizeY - 2 ], Site::EAST_DOWN);

            m_vSites[ i ]->setNeigh( m_vSites[ 0]);
            m_vSites[ i ]->setNeighPosition( m_vSites[ 0], Site::WEST_DOWN);

            m_vSites[ i ]->setNeigh( m_vSites[ getSize() - m_iSizeY + 1]);
            m_vSites[ i ]->setNeighPosition( m_vSites[ getSize() - m_iSizeY + 1], Site::WEST);

            m_vSites[ i ]->setNeigh( m_vSites[ i - 2 ]);
            m_vSites[ i ]->setNeighPosition( m_vSites[ i - 2 ], Site::EAST);

            m_vSites[ i ]->setNeigh( m_vSites[ getSize() - m_iSizeY - 2]);
            m_vSites[ i ]->setNeighPosition( m_vSites[ getSize() - m_iSizeY - 2], Site::EAST_UP);

            m_vSites[ i ]->setNeigh( m_vSites[ getSize() - 2*m_iSizeY ]);
            m_vSites[ i ]->setNeighPosition( m_vSites[ getSize() - 2*m_iSizeY ], Site::WEST_UP);

            m_vSites[ i ]->setNeighPosition( m_vSites[ getSize() - 2*m_iSizeY - 1 ], Site::NORTH );
            m_vSites[ i ]->setNeighPosition( m_vSites[ 2*m_iSizeY - 1  ], Site::SOUTH );

            m_vSites[ i ]->storeActivationSite( m_vSites[ i - 1 ], Site::ACTV_EAST );
            m_vSites[ i ]->storeActivationSite( m_vSites[ getSize() - m_iSizeY ], Site::ACTV_WEST );
            m_vSites[ i ]->storeActivationSite( m_vSites[ m_iSizeY  - 1 ], Site::ACTV_SOUTH );
            m_vSites[ i ]->storeActivationSite( m_vSites[ getSize() - m_iSizeY - 1 ], Site::ACTV_NORTH );
        }
        else{
            m_vSites[ i ]->setNeighPosition( m_vSites[ i - 2*m_iSizeY ], Site::NORTH );
            m_vSites[ i ]->setNeighPosition( m_vSites[ i - (m_iSizeX-2)*m_iSizeY ], Site::SOUTH );

            m_vSites[ i ]->storeActivationSite( m_vSites[ i - 1 ], Site::ACTV_EAST );
            m_vSites[ i ]->storeActivationSite( m_vSites[ i + 1], Site::ACTV_WEST );
            m_vSites[ i ]->storeActivationSite( m_vSites[ i - (m_iSizeX-1)*m_iSizeY ], Site::ACTV_SOUTH );
            m_vSites[ i ]->storeActivationSite( m_vSites[ i - m_iSizeY ], Site::ACTV_NORTH );

            if ( i == (getSize() - m_iSizeY + 1 ) ){
                m_vSites[ i ]->setNeigh( m_vSites[ getSize() - 1 ]);
                m_vSites[ i ]->setNeighPosition( m_vSites[ getSize() - 1 ], Site::EAST);

                m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeY - 1]);
                m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeY - 1], Site::EAST_UP);

                m_vSites[ i ]->setNeigh( m_vSites[ 0]);
                m_vSites[ i ]->setNeighPosition( m_vSites[0], Site::EAST_DOWN);

                m_vSites[ i ]->setNeigh( m_vSites[ i + 2 ]);
                m_vSites[ i ]->setNeighPosition( m_vSites[ i + 2 ], Site::WEST);

                m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeY + 1 ]);
                m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeY + 1 ], Site::WEST_UP);

                m_vSites[ i ]->setNeigh( m_vSites[ i - (m_iSizeX-1)*m_iSizeY + 1]);
                m_vSites[ i ]->setNeighPosition( m_vSites[ i - (m_iSizeX-1)*m_iSizeY + 1 ], Site::WEST_DOWN);
            }
            else if ( i == getSize() - 2 ){
                m_vSites[ i ]->setNeigh( m_vSites[ getSize() - m_iSizeY ]);
                m_vSites[ i ]->setNeighPosition( m_vSites[ getSize() - m_iSizeY ], Site::WEST);

                m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeY + 1]);
                m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeY + 1  ], Site::WEST_UP);

                m_vSites[ i ]->setNeigh( m_vSites[ i - (m_iSizeX-1)*m_iSizeY + 1 ]);
                m_vSites[ i ]->setNeighPosition( m_vSites[ i - (m_iSizeX-1)*m_iSizeY + 1 ], Site::WEST_DOWN);

                m_vSites[ i ]->setNeigh( m_vSites[ i - 2 ]);
                m_vSites[ i ]->setNeighPosition( m_vSites[ i - 2 ], Site::EAST);

                m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeY - 1]);
                m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeY - 1], Site::EAST_UP);

                m_vSites[ i ]->setNeigh( m_vSites[  i - (m_iSizeX-1)*m_iSizeY - 1]);
                m_vSites[ i ]->setNeighPosition( m_vSites[ i - (m_iSizeX-1)*m_iSizeY - 1 ], Site::EAST_DOWN);
            }
            else{
                m_vSites[ i ]->setNeigh( m_vSites[ i + 2 ]);
                m_vSites[ i ]->setNeighPosition( m_vSites[ i + 2 ], Site::WEST);

                m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeY + 1]);
                m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeY + 1 ], Site::WEST_UP);

                m_vSites[ i ]->setNeigh( m_vSites[ i - (m_iSizeX-1)*m_iSizeY + 1]);
                m_vSites[ i ]->setNeighPosition( m_vSites[ i - (m_iSizeX-1)*m_iSizeY + 1 ], Site::WEST_DOWN);

                m_vSites[ i ]->setNeigh( m_vSites[ i - 2 ]);
                m_vSites[ i ]->setNeighPosition( m_vSites[ i - 2 ], Site::EAST);

                m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeY - 1]);
                m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeY - 1], Site::EAST_UP);

                m_vSites[ i ]->setNeigh( m_vSites[  i - (m_iSizeX-1)*m_iSizeY - 1]);
                m_vSites[ i ]->setNeighPosition( m_vSites[ i - (m_iSizeX-1)*m_iSizeY - 1 ], Site::EAST_DOWN);
            }
        }
    }

    /* First East column */
    for (int i = m_iSizeY; i < (getSize() - m_iSizeY); i+= m_iSizeY ){
        m_vSites[ i ]->setNeigh( m_vSites[ i + m_iSizeY - 2]);
        m_vSites[ i ]->setNeighPosition(m_vSites[ i + m_iSizeY - 2], Site::EAST);

        m_vSites[ i ]->setNeigh( m_vSites[ i + 2]);
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + 2], Site::WEST);

        m_vSites[ i ]->setNeigh( m_vSites[ i + 2*m_iSizeY - 1]);
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + 2*m_iSizeY - 1 ], Site::EAST_DOWN);

        m_vSites[ i ]->setNeigh( m_vSites[ i + m_iSizeY + 1 ]);
        m_vSites[ i ]->setNeighPosition( m_vSites[  i + m_iSizeY + 1 ], Site::WEST_DOWN);

        m_vSites[ i ]->setNeigh( m_vSites[ i - 1]);
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - 1 ], Site::EAST_UP);

        m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeY + 1]);
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeY + 1 ], Site::WEST_UP);

        m_vSites[ i ]->storeActivationSite( m_vSites[ i + m_iSizeY - 1 ], Site::ACTV_EAST );
        m_vSites[ i ]->storeActivationSite( m_vSites[ i + 1], Site::ACTV_WEST );
        m_vSites[ i ]->storeActivationSite( m_vSites[ i + m_iSizeY ], Site::ACTV_SOUTH );
        m_vSites[ i ]->storeActivationSite( m_vSites[ i - m_iSizeY ], Site::ACTV_NORTH );

        if ( i == m_iSizeY ){
            m_vSites[ i ]->setNeighPosition( m_vSites[ getSize() - m_iSizeY ], Site::NORTH);
            m_vSites[ i ]->setNeighPosition( m_vSites[ i + 2*m_iSizeY ], Site::SOUTH);
        }
        else if ( i == getSize() - 2*m_iSizeY ){
            m_vSites[ i ]->setNeighPosition( m_vSites[ i - 2*m_iSizeY ], Site::NORTH);
            m_vSites[ i ]->setNeighPosition( m_vSites[ 0 ], Site::SOUTH);
        }
        else{
            m_vSites[ i ]->setNeighPosition( m_vSites[ i - 2*m_iSizeY ], Site::NORTH);
            m_vSites[ i ]->setNeighPosition( m_vSites[ i + 2*m_iSizeY ], Site::SOUTH);
        }
    }

    /* First West column side */
    for (int i = (2*m_iSizeY -1); i < (getSize() - m_iSizeY); i+= m_iSizeY ){
        m_vSites[ i ]->setNeigh( m_vSites[ i - 2]);
        m_vSites[ i ]->setNeighPosition(m_vSites[ i - 2], Site::EAST);

        m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeY +  2]);
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeY +  2 ], Site::WEST);

        m_vSites[ i ]->setNeigh( m_vSites[ i + m_iSizeY - 1]);
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeY - 1 ], Site::EAST_DOWN);

        m_vSites[ i ]->setNeigh( m_vSites[ i + 1 ]);
        m_vSites[ i ]->setNeighPosition( m_vSites[  i + 1 ], Site::WEST_DOWN);

        m_vSites[ i ]->setNeigh( m_vSites[ i -m_iSizeY - 1]);
        m_vSites[ i ]->setNeighPosition( m_vSites[ i -m_iSizeY - 1 ], Site::EAST_UP);

        m_vSites[ i ]->setNeigh( m_vSites[ i - 2*m_iSizeY + 1]);
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - 2*m_iSizeY + 1 ], Site::WEST_UP);

        m_vSites[ i ]->storeActivationSite( m_vSites[ i - 1 ], Site::ACTV_EAST );
        m_vSites[ i ]->storeActivationSite( m_vSites[ i - m_iSizeY + 1], Site::ACTV_WEST );
        m_vSites[ i ]->storeActivationSite( m_vSites[ i + m_iSizeY ], Site::ACTV_SOUTH );
        m_vSites[ i ]->storeActivationSite( m_vSites[ i - m_iSizeY ], Site::ACTV_NORTH );

        if ( i == 2*m_iSizeY - 1 ){
            m_vSites[ i ]->setNeighPosition( m_vSites[ getSize() - 1 ], Site::NORTH);
            m_vSites[ i ]->setNeighPosition( m_vSites[ i + 2*m_iSizeY ], Site::SOUTH);
        }
        else if ( i == getSize() - m_iSizeY - 1 ){
            m_vSites[ i ]->setNeighPosition( m_vSites[ i - 2*m_iSizeY ], Site::NORTH);
            m_vSites[ i ]->setNeighPosition( m_vSites[ m_iSizeY - 1 ], Site::SOUTH);
        }
        else{
            m_vSites[ i ]->setNeighPosition( m_vSites[ i - 2*m_iSizeY ], Site::NORTH);
            m_vSites[ i ]->setNeighPosition( m_vSites[ i + 2*m_iSizeY ], Site::SOUTH);
        }
    }

    /* Second row */
    for ( int i = (m_iSizeY + 1); i < (2*m_iSizeY - 1); i++){
        m_vSites[ i ]->setNeigh( m_vSites[ i + m_iSizeY - 1]);
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeY - 1 ], Site::EAST_DOWN);

        m_vSites[ i ]->setNeigh( m_vSites[ i + m_iSizeY + 1 ]);
        m_vSites[ i ]->setNeighPosition( m_vSites[  i + m_iSizeY + 1 ], Site::WEST_DOWN);

        m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeY - 1]);
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeY - 1 ], Site::EAST_UP);

        m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeY + 1]);
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeY + 1 ], Site::WEST_UP);

        m_vSites[ i ]->storeActivationSite( m_vSites[ i - 1 ], Site::ACTV_EAST );
        m_vSites[ i ]->storeActivationSite( m_vSites[ i + 1], Site::ACTV_WEST );
        m_vSites[ i ]->storeActivationSite( m_vSites[ i + m_iSizeY ], Site::ACTV_SOUTH );
        m_vSites[ i ]->storeActivationSite( m_vSites[ i - m_iSizeY ], Site::ACTV_NORTH );

        m_vSites[ i ]->setNeighPosition( m_vSites[ i + (m_iSizeX - 2)*m_iSizeY ], Site::NORTH);
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + 2*m_iSizeY ], Site::SOUTH);

        if ( i == (m_iSizeY + 1)  ){
            m_vSites[ i ]->setNeigh( m_vSites[ i + m_iSizeY - 2]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ i + m_iSizeY - 2], Site::EAST);

            m_vSites[ i ]->setNeigh( m_vSites[ i + 2]);
            m_vSites[ i ]->setNeighPosition( m_vSites[ i + 2 ], Site::WEST);
        }
        else if ( i == 2*m_iSizeY - 2 ){
            m_vSites[ i ]->setNeigh( m_vSites[ i - 2]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ i - 2], Site::EAST);

            m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeY + 2]);
            m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeY + 2 ], Site::WEST);
        }
        else{
            m_vSites[ i ]->setNeigh( m_vSites[ i - 2]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ i - 2], Site::EAST);

            m_vSites[ i ]->setNeigh( m_vSites[ i + 2]);
            m_vSites[ i ]->setNeighPosition( m_vSites[ i + 2 ], Site::WEST);
        }
    }

    /* Second before end row */
    for ( int i = ( getSize() - 2*m_iSizeY + 1 ); i < (getSize() - m_iSizeY - 1); i++){
        m_vSites[ i ]->setNeigh( m_vSites[ i + m_iSizeY - 1]);
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeY - 1 ], Site::EAST_DOWN);

        m_vSites[ i ]->setNeigh( m_vSites[ i + m_iSizeY + 1 ]);
        m_vSites[ i ]->setNeighPosition( m_vSites[  i + m_iSizeY + 1 ], Site::WEST_DOWN);

        m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeY - 1]);
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeY - 1 ], Site::EAST_UP);

        m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeY + 1]);
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeY + 1 ], Site::WEST_UP);

        m_vSites[ i ]->storeActivationSite( m_vSites[ i - 1 ], Site::ACTV_EAST );
        m_vSites[ i ]->storeActivationSite( m_vSites[ i + 1], Site::ACTV_WEST );
        m_vSites[ i ]->storeActivationSite( m_vSites[ i + m_iSizeY ], Site::ACTV_SOUTH );
        m_vSites[ i ]->storeActivationSite( m_vSites[ i - m_iSizeY ], Site::ACTV_NORTH );

        m_vSites[ i ]->setNeighPosition( m_vSites[ i - 2*m_iSizeY  ], Site::NORTH);
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - (m_iSizeX - 2)*m_iSizeY ], Site::SOUTH);

        if ( i == ( getSize() - 2*m_iSizeY + 1 )  ){
            m_vSites[ i ]->setNeigh( m_vSites[ i + m_iSizeY - 2]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ i + m_iSizeY - 2], Site::EAST);

            m_vSites[ i ]->setNeigh( m_vSites[ i + 2]);
            m_vSites[ i ]->setNeighPosition( m_vSites[ i + 2 ], Site::WEST);
        }
        else if ( i == getSize() - m_iSizeY - 2 ){
            m_vSites[ i ]->setNeigh( m_vSites[ i - 2]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ i - 2], Site::EAST);

            m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeY + 2]);
            m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeY + 2 ], Site::WEST);
        }
        else{
            m_vSites[ i ]->setNeigh( m_vSites[ i - 2]);
            m_vSites[ i ]->setNeighPosition(m_vSites[ i - 2], Site::EAST);

            m_vSites[ i ]->setNeigh( m_vSites[ i + 2]);
            m_vSites[ i ]->setNeighPosition( m_vSites[ i + 2 ], Site::WEST);
        }
    }

    /* Second east column */
    /* Since we have defined the corners in the "second row" and in the "second before end row" no extra action is needed - No ifs...*/
    for ( int i = ( 2*m_iSizeY + 1 ); i < ( (m_iSizeX-3)*m_iSizeY + 2 ); i+=m_iSizeY){
        m_vSites[ i ]->setNeigh( m_vSites[ i + m_iSizeY - 2]);
        m_vSites[ i ]->setNeighPosition(m_vSites[ i + m_iSizeY - 2], Site::EAST);

        m_vSites[ i ]->setNeigh( m_vSites[ i + 2]);
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + 2 ], Site::WEST);

        m_vSites[ i ]->setNeigh( m_vSites[ i + m_iSizeY - 1]);
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeY - 1 ], Site::EAST_DOWN);

        m_vSites[ i ]->setNeigh( m_vSites[ i + m_iSizeY + 1 ]);
        m_vSites[ i ]->setNeighPosition( m_vSites[  i + m_iSizeY + 1 ], Site::WEST_DOWN);

        m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeY - 1]);
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeY - 1 ], Site::EAST_UP);

        m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeY + 1]);
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeY + 1 ], Site::WEST_UP);

        m_vSites[ i ]->storeActivationSite( m_vSites[ i - 1 ], Site::ACTV_EAST );
        m_vSites[ i ]->storeActivationSite( m_vSites[ i + 1], Site::ACTV_WEST );
        m_vSites[ i ]->storeActivationSite( m_vSites[ i + m_iSizeY ], Site::ACTV_SOUTH );
        m_vSites[ i ]->storeActivationSite( m_vSites[ i - m_iSizeY ], Site::ACTV_NORTH );

        m_vSites[ i ]->setNeighPosition( m_vSites[ i - 2*m_iSizeY  ], Site::NORTH);
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + 2*m_iSizeY ], Site::SOUTH);
    }

    /* Second west column */
    /* Since we have defined the corners in the "second row" and in the "second before end row" no extra action is needed - No ifs...*/
    for ( int i = (3*m_iSizeY - 2) ; i < ( (m_iSizeX-1)*m_iSizeY - 2 ); i+=m_iSizeY){
        m_vSites[ i ]->setNeigh( m_vSites[ i - 2]);
        m_vSites[ i ]->setNeighPosition(m_vSites[ i - 2], Site::EAST);

        m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeY + 2]);
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeY + 2 ], Site::WEST);

        m_vSites[ i ]->setNeigh( m_vSites[ i + m_iSizeY - 1]);
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + m_iSizeY - 1 ], Site::EAST_DOWN);

        m_vSites[ i ]->setNeigh( m_vSites[ i + m_iSizeY + 1 ]);
        m_vSites[ i ]->setNeighPosition( m_vSites[  i + m_iSizeY + 1 ], Site::WEST_DOWN);

        m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeY - 1]);
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeY - 1 ], Site::EAST_UP);

        m_vSites[ i ]->setNeigh( m_vSites[ i - m_iSizeY + 1]);
        m_vSites[ i ]->setNeighPosition( m_vSites[ i - m_iSizeY + 1 ], Site::WEST_UP);

        m_vSites[ i ]->storeActivationSite( m_vSites[ i - 1 ], Site::ACTV_EAST );
        m_vSites[ i ]->storeActivationSite( m_vSites[ i + 1], Site::ACTV_WEST );
        m_vSites[ i ]->storeActivationSite( m_vSites[ i + m_iSizeY ], Site::ACTV_SOUTH );
        m_vSites[ i ]->storeActivationSite( m_vSites[ i - m_iSizeY ], Site::ACTV_NORTH );

        m_vSites[ i ]->setNeighPosition( m_vSites[ i - 2*m_iSizeY  ], Site::NORTH);
        m_vSites[ i ]->setNeighPosition( m_vSites[ i + 2*m_iSizeY ], Site::SOUTH);
    }
}

void FCC::check()
{
    int k = 0;

    cout << "Checking lattice..." << endl;

    int test = 2;
    cout << test << ": ";
    cout << "W:" <<  getSite( test )->getNeighPosition( Site::WEST)->getID() << " ";
    cout << "Wu:" << getSite( test )->getNeighPosition( Site::WEST_UP)->getID() << " ";
    cout << "WD:" << getSite( test )->getNeighPosition( Site::WEST_DOWN)->getID() << " ";
    cout << "E:" << getSite( test )->getNeighPosition( Site::EAST)->getID() << " ";
    cout << "EU:" << getSite( test )->getNeighPosition( Site::EAST_UP)->getID() << " ";
    cout << "ED:" << getSite( test )->getNeighPosition( Site::EAST_DOWN)->getID() << " ";
    cout << "N:" <<getSite( test )->getNeighPosition( Site::NORTH)->getID() << " ";
    cout << "S:" <<getSite( test )->getNeighPosition( Site::SOUTH)->getID() << endl;

    cout << "Activation: " << endl;

    cout << "N:" <<  getSite( test )->getActivationSite(Site::ACTV_NORTH)->getID() << " ";
    cout << "S:" << getSite( test )->getActivationSite(Site::ACTV_SOUTH)->getID() << " ";
    cout << "E:" << getSite( test )->getActivationSite(Site::ACTV_EAST)->getID() << " ";
    cout << "W:" << getSite( test )->getActivationSite(Site::ACTV_WEST)->getID() << endl;
}
