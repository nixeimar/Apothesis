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

FCC::FCC(Apothesis *apothesis, bool step, vector<int> stepInfo) : Lattice(apothesis),
m_hasSteps(step),
m_stepInfo(stepInfo)
{
  ;
}

void FCC::setInitialHeight( int  height ) { m_iHeight = height; }

void FCC::build()
{
  
  if ( m_iSizeX%2 != 0 || m_iSizeX%2 != 0){
    cout << "The size of the lattice must be an even number in each direction." << endl;
    EXIT;
    }

  if ( m_iSizeX == 0 || m_iSizeY == 0) {
    m_errorHandler->error_simple_msg("The lattice size cannot be zero in either dimension.");
    EXIT;
    }

  if ( m_iHeight < 5) {
    m_errorHandler->warningSimple_msg("The lattice initial height is too small.Consider revising.");
    }

  // The sites of the lattice.
  m_vSites.resize( getSize() );
  for ( int i = 0; i < m_vSites.size(); i++)
    m_vSites[ i ] = new Site(this);

  //  m_pSites = new Site[ m_iSizeX*m_iSizeY];

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
      

  mf_neigh();
}

FCC::~FCC()
{
  for ( int i = 0; i<getSize(); i++)
    delete m_vSites[i];
}


void FCC::mf_neigh()
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

Site* FCC::getSite( int id){ return m_vSites[ id]; }

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
