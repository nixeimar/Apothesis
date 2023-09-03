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

#include "diamond.h"

void Diamond::setInitialHeight( int height ) { m_iHeight = height; }

void Diamond::build(){
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

    // Here we set the label of the species
    for (int i = 0; i < m_iSizeY; i++)
    {
        for (int j = 0; j < m_iSizeX; j++)
            m_vSites[i * m_iSizeX + j]->setLabel(m_sLabel);
    }
}

/// Builds a  stepped surface
void Diamond::buildSteps(){}

void Diamond::writeXYZ(string){}

void Diamond::mf_neigh() {

    // The sites of the lattice.
    m_vSites.resize(getSize());

    int L = getSize();

    for (int site = 0; site < getSize(); site++) {
        int x = site / L;
        int y = site % L;

        m_vSites[ site ]->setNeigh( m_vSites[ ((x + 1) % L) * L + y ]); //Right
        m_vSites[ site ]->setNeigh( m_vSites[ x * L + ((y + 1) % L) ]); //Down
        m_vSites[ site ]->setNeigh( m_vSites[ ((x - 1 + L) % L) * L + y ]); //Left
        m_vSites[ site ]->setNeigh( m_vSites[ x * L + ((y - 1 + L) % L) ]); //Up
    }
}

unordered_map<string, double> Diamond::computeCoverages( vector<string> species){ }
