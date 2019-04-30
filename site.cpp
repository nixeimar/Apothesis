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

#ifndef SITE_CPP
#define SITE_CPP

#include "site.h"

namespace SurfaceTiles{

Site::Site(){;}

Site::~Site(){;}

vector<Site* > Site::getNeighs() { return m_vNeigh; }

void Site::setNeigh( Site* s){
  m_vNeigh.push_back( s);
  }

void Site::setID(int id) { m_iID = id;}
int Site::getID() { return m_iID;}

void Site::setHeight(int h){ m_iHeight = h;}

int Site::getHeight() { return m_iHeight; }

void Site::setNeighboursNum( int num){ m_iNumNeighs = num;
                                       }
void Site::setNeighPosition( Site* s, NeighPoisition np){ m_mapNeigh[ np] = s; }

Site* Site::getNeighPosition( NeighPoisition np){ return m_mapNeigh[ np];}

void Site::storeActivationSite( Site* s, ActivationSite as){ m_mapAct[ as] = s;}

Site* Site::getActivationSite( ActivationSite as){ return m_mapAct[ as];}

void Site::addProcess(Process* process) { m_lProcs.push_back( process ); }

void Site::removeProcess(Process* process) { m_lProcs.remove( process ); }

}


#endif
