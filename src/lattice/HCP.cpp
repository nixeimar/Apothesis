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
    ;
}

void HCP::setInitialHeight(int height) { m_iHeight = height; }

void HCP::build()
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
    for (int i = 1; i < m_iSizeY - 1; i++)
    {
        for (int j = 1; j < m_iSizeX - 1; j++)
        {
            m_vSites[i * m_iSizeX + j]->setNeigh(m_vSites[(i - 1) * m_iSizeX + j]);
            m_vSites[i * m_iSizeX + j]->setNeighPosition(m_vSites[(i - 1) * m_iSizeX + j], Site::WEST_UP);

            m_vSites[i * m_iSizeX + j]->setNeigh(m_vSites[(i + 1) * m_iSizeX + j]);
            m_vSites[i * m_iSizeX + j]->setNeighPosition(m_vSites[(i + 1) * m_iSizeX + j], Site::WEST_DOWN);

            m_vSites[i * m_iSizeX + j]->setNeigh(m_vSites[i * m_iSizeX + j + 1]);
            m_vSites[i * m_iSizeX + j]->setNeighPosition(m_vSites[i * m_iSizeX + j + 1], Site::EAST);

            m_vSites[i * m_iSizeX + j]->setNeigh(m_vSites[i * m_iSizeX + j - 1]);
            m_vSites[i * m_iSizeX + j]->setNeighPosition(m_vSites[i * m_iSizeX + j - 1], Site::WEST);

            m_vSites[i * m_iSizeX + j]->setNeigh(m_vSites[(i - 1) * m_iSizeX + j + 1]);
            m_vSites[i * m_iSizeX + j]->setNeighPosition(m_vSites[(i - 1) * m_iSizeX + j + 1], Site::EAST_UP);

            m_vSites[i * m_iSizeX + j]->setNeigh(m_vSites[(i + 1) * m_iSizeX + j + 1]);
            m_vSites[i * m_iSizeX + j]->setNeighPosition(m_vSites[(i + 1) * m_iSizeX + j + 1], Site::EAST_DOWN);
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
            m_vSites[j]->setNeighPosition(m_vSites[j + m_iSizeX], Site::EAST_DOWN);

            m_vSites[j]->setNeigh(m_vSites[iThirdCorner + j]);
            m_vSites[j]->setNeighPosition(m_vSites[iThirdCorner + j], Site::EAST_UP);

            m_vSites[j]->setNeigh(m_vSites[iThirdCorner + j - 1]);
            m_vSites[j]->setNeighPosition(m_vSites[iThirdCorner + j - 1], Site::WEST_UP);

            m_vSites[j]->setNeigh(m_vSites[m_iSizeX + j - 1]);
            m_vSites[j]->setNeighPosition(m_vSites[m_iSizeX + j - 1], Site::WEST_DOWN);
        }
        else if (j == iFirstCorner)
        {
            m_vSites[j]->setNeigh(m_vSites[iSecondCorner]);
            m_vSites[j]->setNeighPosition(m_vSites[iSecondCorner], Site::WEST);

            m_vSites[j]->setNeigh(m_vSites[1]);
            m_vSites[j]->setNeighPosition(m_vSites[1], Site::EAST);

            m_vSites[j]->setNeigh(m_vSites[iSecondCorner + 1]);
            m_vSites[j]->setNeighPosition(m_vSites[iSecondCorner + 1], Site::EAST_DOWN);

            m_vSites[j]->setNeigh(m_vSites[iThirdCorner]);
            m_vSites[j]->setNeighPosition(m_vSites[iThirdCorner], Site::EAST_UP);

            m_vSites[j]->setNeigh(m_vSites[iForthCorner]);
            m_vSites[j]->setNeighPosition(m_vSites[iForthCorner], Site::WEST_UP);

            m_vSites[j]->setNeigh(m_vSites[2 * iSecondCorner + 1]);
            m_vSites[j]->setNeighPosition(m_vSites[2 * iSecondCorner + 1], Site::WEST_DOWN);
        }
        else if (j == iSecondCorner)
        {
            m_vSites[j]->setNeigh(m_vSites[j - 1]);
            m_vSites[j]->setNeighPosition(m_vSites[j - 1], Site::WEST);

            m_vSites[j]->setNeigh(m_vSites[0]);
            m_vSites[j]->setNeighPosition(m_vSites[0], Site::EAST);

            m_vSites[j]->setNeigh(m_vSites[2 * m_iSizeX - 1]);
            m_vSites[j]->setNeighPosition(m_vSites[2 * m_iSizeX - 1], Site::EAST_DOWN);

            m_vSites[j]->setNeigh(m_vSites[iForthCorner]);
            m_vSites[j]->setNeighPosition(m_vSites[iForthCorner], Site::EAST_UP);

            m_vSites[j]->setNeigh(m_vSites[iForthCorner]);
            m_vSites[j]->setNeighPosition(m_vSites[iForthCorner - 1], Site::WEST_UP);

            m_vSites[j]->setNeigh(m_vSites[2 * iSecondCorner]);
            m_vSites[j]->setNeighPosition(m_vSites[2 * iSecondCorner], Site::WEST_DOWN);
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
            m_vSites[j]->setNeighPosition(m_vSites[iFirstCorner + iPos], Site::WEST_DOWN);

            m_vSites[j]->setNeigh(m_vSites[j - m_iSizeX]);
            m_vSites[j]->setNeighPosition(m_vSites[j - m_iSizeX], Site::WEST_UP);

            m_vSites[j]->setNeigh(m_vSites[j - m_iSizeX + 1]);
            m_vSites[j]->setNeighPosition(m_vSites[j - m_iSizeX + 1], Site::EAST_UP);

            m_vSites[j]->setNeigh(m_vSites[iFirstCorner + iPos + 1]);
            m_vSites[j]->setNeighPosition(m_vSites[iFirstCorner + iPos + 1], Site::WEST_DOWN);

            iPos++;
        }
        else if (j == iThirdCorner)
        {
            m_vSites[j]->setNeigh(m_vSites[iForthCorner]);
            m_vSites[j]->setNeighPosition(m_vSites[iForthCorner], Site::WEST);

            m_vSites[j]->setNeigh(m_vSites[iThirdCorner + 1]);
            m_vSites[j]->setNeighPosition(m_vSites[iThirdCorner + 1], Site::EAST);

            m_vSites[j]->setNeigh(m_vSites[iFirstCorner]);
            m_vSites[j]->setNeighPosition(m_vSites[iFirstCorner], Site::WEST_DOWN);

            m_vSites[j]->setNeigh(m_vSites[iThirdCorner - m_iSizeX]);
            m_vSites[j]->setNeighPosition(m_vSites[iThirdCorner - m_iSizeX], Site::WEST_UP);

            m_vSites[j]->setNeigh(m_vSites[iThirdCorner - m_iSizeX + 1]);
            m_vSites[j]->setNeighPosition(m_vSites[iThirdCorner - m_iSizeX + 1], Site::EAST_UP);

            m_vSites[j]->setNeigh(m_vSites[iFirstCorner + 1]);
            m_vSites[j]->setNeighPosition(m_vSites[iFirstCorner + 1], Site::EAST_DOWN);
        }
        else if (j == iForthCorner)
        {
            m_vSites[j]->setNeigh(m_vSites[iForthCorner - 1]);
            m_vSites[j]->setNeighPosition(m_vSites[iForthCorner - 1], Site::WEST);

            m_vSites[j]->setNeigh(m_vSites[iThirdCorner]);
            m_vSites[j]->setNeighPosition(m_vSites[iThirdCorner], Site::EAST);

            m_vSites[j]->setNeigh(m_vSites[iSecondCorner]);
            m_vSites[j]->setNeighPosition(m_vSites[iSecondCorner], Site::WEST_DOWN);

            m_vSites[j]->setNeigh(m_vSites[iThirdCorner - 1]);
            m_vSites[j]->setNeighPosition(m_vSites[iThirdCorner - 1], Site::WEST_UP);

            m_vSites[j]->setNeigh(m_vSites[iThirdCorner - m_iSizeX]);
            m_vSites[j]->setNeighPosition(m_vSites[iThirdCorner - m_iSizeX], Site::EAST_UP);

            m_vSites[j]->setNeigh(m_vSites[iFirstCorner]);
            m_vSites[j]->setNeighPosition(m_vSites[iFirstCorner], Site::EAST_DOWN);
        }
    }

    /* First column*/
    for (int j = iFirstCorner + m_iSizeX; j < iThirdCorner; j += m_iSizeX)
    {
        if (j % 20 == 0)
        {
            m_vSites[j]->setNeigh(m_vSites[j + m_iSizeX - 1]);
            m_vSites[j]->setNeighPosition(m_vSites[j + m_iSizeX - 1], Site::WEST);

            m_vSites[j]->setNeigh(m_vSites[j + 1]);
            m_vSites[j]->setNeighPosition(m_vSites[j + 1], Site::EAST);

            m_vSites[j]->setNeigh(m_vSites[j + (2 * m_iSizeX) - 1]);
            m_vSites[j]->setNeighPosition(m_vSites[j + (2 * m_iSizeX) - 1], Site::WEST_DOWN);

            m_vSites[j]->setNeigh(m_vSites[j - 1]);
            m_vSites[j]->setNeighPosition(m_vSites[j - 1], Site::WEST_UP);

            m_vSites[j]->setNeigh(m_vSites[j + m_iSizeX]);
            m_vSites[j]->setNeighPosition(m_vSites[j + m_iSizeX], Site::EAST_DOWN);

            m_vSites[j]->setNeigh(m_vSites[j - m_iSizeX]);
            m_vSites[j]->setNeighPosition(m_vSites[j - m_iSizeX], Site::EAST_UP);
        }
        else
        {
            m_vSites[j]->setNeigh(m_vSites[j + m_iSizeX - 1]);
            m_vSites[j]->setNeighPosition(m_vSites[j + m_iSizeX - 1], Site::WEST);

            m_vSites[j]->setNeigh(m_vSites[j + 1]);
            m_vSites[j]->setNeighPosition(m_vSites[j + 1], Site::EAST);

            m_vSites[j]->setNeigh(m_vSites[j + m_iSizeX]);
            m_vSites[j]->setNeighPosition(m_vSites[j + m_iSizeX], Site::WEST_DOWN);

            m_vSites[j]->setNeigh(m_vSites[j - m_iSizeX]);
            m_vSites[j]->setNeighPosition(m_vSites[j - m_iSizeX], Site::WEST_UP);

            m_vSites[j]->setNeigh(m_vSites[j + m_iSizeX + 1]);
            m_vSites[j]->setNeighPosition(m_vSites[j + m_iSizeX + 1], Site::EAST_DOWN);

            m_vSites[j]->setNeigh(m_vSites[j - m_iSizeX + 1]);
            m_vSites[j]->setNeighPosition(m_vSites[j - m_iSizeX + 1], Site::EAST_UP);
        }
    }

    /* Last column*/
    for (int j = iSecondCorner + m_iSizeX; j < iForthCorner; j += m_iSizeX)
    {
        if ((j + 1) % 20 == 0)
        {
            m_vSites[j]->setNeigh(m_vSites[j - 1]);
            m_vSites[j]->setNeighPosition(m_vSites[j - 1], Site::WEST);

            m_vSites[j]->setNeigh(m_vSites[j - m_iSizeX + 1]);
            m_vSites[j]->setNeighPosition(m_vSites[j - m_iSizeX + 1], Site::EAST);

            m_vSites[j]->setNeigh(m_vSites[j + m_iSizeX]);
            m_vSites[j]->setNeighPosition(m_vSites[j + m_iSizeX], Site::WEST_DOWN);

            m_vSites[j]->setNeigh(m_vSites[j - m_iSizeX]);
            m_vSites[j]->setNeighPosition(m_vSites[j - m_iSizeX], Site::WEST_UP);

            m_vSites[j]->setNeigh(m_vSites[j - (2 * m_iSizeX) + 1]);
            m_vSites[j]->setNeighPosition(m_vSites[j - (2 * m_iSizeX) + 1], Site::EAST_UP);

            m_vSites[j]->setNeigh(m_vSites[j + 1]);
            m_vSites[j]->setNeighPosition(m_vSites[j + 1], Site::EAST_DOWN);
        }
        else
        {
            m_vSites[j]->setNeigh(m_vSites[j - 1]);
            m_vSites[j]->setNeighPosition(m_vSites[j - 1], Site::WEST);

            m_vSites[j]->setNeigh(m_vSites[j - m_iSizeX + 1]);
            m_vSites[j]->setNeighPosition(m_vSites[j - m_iSizeX + 1], Site::EAST);

            m_vSites[j]->setNeigh(m_vSites[j + m_iSizeX]);
            m_vSites[j]->setNeighPosition(m_vSites[j + m_iSizeX], Site::EAST_DOWN);

            m_vSites[j]->setNeigh(m_vSites[j - m_iSizeX]);
            m_vSites[j]->setNeighPosition(m_vSites[j - m_iSizeX], Site::EAST_UP);

            m_vSites[j]->setNeigh(m_vSites[j - (2 * m_iSizeX) + 1]);
            m_vSites[j]->setNeighPosition(m_vSites[j - (2 * m_iSizeX) + 1], Site::WEST_UP);

            m_vSites[j]->setNeigh(m_vSites[j + m_iSizeX - 1]);
            m_vSites[j]->setNeighPosition(m_vSites[j + m_iSizeX - 1], Site::WEST_DOWN);
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
