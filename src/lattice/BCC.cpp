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
#include "read.h"

BCC::BCC(Apothesis *apothesis) : Lattice(apothesis)
{
	;
}

// TODO: Should "hasSteps" be migrated to lattice base class?
BCC::BCC(Apothesis *apothesis, bool step, vector<int> stepInfo) : Lattice(apothesis),
m_hasSteps(step),
m_stepInfo(stepInfo)
{
	;
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
	m_vSites.resize(getSize());
	for (int i = 0; i < m_vSites.size(); i++)
		m_vSites[i] = new Site(this);

	//  m_pSites = new Site[ m_iSizeX*m_iSizeY];

	for (int i = 0; i < m_iSizeX; i++)
	{
		for (int j = i * m_iSizeY; j < (m_iSizeY + i * m_iSizeY); j++)
		{
			m_vSites[j]->setID(j);
			m_vSites[j]->setHeight(m_iHeight - 1);
			m_vSites[j]->setLatticeType(Site::LatticeType::BCC);
		}
	}

	for (int i = 0; i < m_iSizeX; i++)
	{
		for (int j = i * m_iSizeY; j < (m_iSizeY + i * m_iSizeY); j++)
			cout << m_vSites[j]->getHeight() << " ";

		cout << " " << endl;
	}

	for (int i = 0; i < m_iSizeX; i++)
	{
		for (int j = i * m_iSizeY; j < (m_iSizeY + i * m_iSizeY); j++)
			cout << m_vSites[j]->getID() << " ";

		cout << " " << endl;
	}

	if (m_hasSteps)
		mf_buildSteps();

	mf_neigh();
}

BCC::~BCC()
{
	for (int i = 0; i < getSize(); i++)
		delete m_vSites[i];
}

void BCC::setSteps(bool hasSteps)
{
	m_hasSteps = hasSteps;
}

void BCC::setStepInfo(int sizeX, int sizeY, int sizeZ)
{
	m_iStepX = sizeX;
	m_iStepY = sizeY;
	m_iStepZ = sizeZ;
}

void BCC::mf_buildSteps()
{
	// Pick dimension of step
	// TODO: Can we assume that the largest value is the dimension of stepping?
	// Find the initial height from arbitrary site
	int initialHeight = m_vSites[0]->getHeight();
	vector<int> currentDimensions{m_iSizeX, m_iSizeY, initialHeight};

	int iteration = 0;
	int stepDimension = 0, stepSoFar = 0, growthDimension = 0, growthSoFar = 0, latentDimension = 0;
	for (auto& dim : m_stepInfo)
	{
		// If the step information is the same value as lattice dim, this will not be the step/growth dimension
		if (dim != currentDimensions[iteration])
		{
			// The step dimension will be the largest value
		    if (dim > stepSoFar) 
			{
				stepDimension = iteration;
				stepSoFar = dim;
			}
			else
			{
				growthDimension = iteration;	
			}
		}
		else
		{
			latentDimension = iteration;
		}
		
		iteration++;
	}

	for (unsigned int firstDim = 0; firstDim < currentDimensions[latentDimension]; ++firstDim)
	{
		for (unsigned int secondDim = 0; secondDim < currentDimensions[stepDimension]; ++secondDim)
		{
			// Calculate how much we increase the step by.
			// Calculation is split up to ensure we have proper integer division in the first step.
			int growth = secondDim/m_stepInfo[stepDimension];
			growth *= m_stepInfo[growthDimension];
			int index = firstDim * currentDimensions[stepDimension] + secondDim;
			m_vSites[index]->increaseHeight(growth);
		}
	}

	/* if (m_iSizeX % m_iStepX != 0)
	{
		m_errorHandler->error_simple_msg("ERROR: The number of steps you provided doesn't conform with the lattice size ");
		exit(0);
	} */
	//if (m_iStepY != 0) // Be sure that we do have steps. If indi_y = 0 (1 0 0) then we have an initial flat surface
	//{
	//	unsigned int steps = m_iSizeX / m_iStepX;
	//	for (unsigned int step = 1; step < steps; step++)
	//		for (unsigned int i = step * m_iStepX; i < (step + 1) * m_iStepX; i++)
	//			for (unsigned int j = 0; j < m_iSizeY; j++)
	//				m_vSites[i * m_iStepY + j] += m_iStepY * step;
	//	//(*mesh)[i][j] += m_iStepY * step;///
	//	cout << "Number of steps:" << steps << endl;
	//}
}

void BCC::mf_neigh()
{
	/* All except the boundaries */
	for (int i = 1; i < m_iSizeY - 1; i++)
	{
		for (int j = 1; j < m_iSizeX - 1; j++)
		{
			m_vSites[i * m_iSizeX + j]->setNeigh(m_vSites[(i - 1) * m_iSizeX + j]);
			m_vSites[i * m_iSizeX + j]->setNeighPosition(m_vSites[(i - 1) * m_iSizeX + j], Site::SOUTH);

			m_vSites[i * m_iSizeX + j]->setNeigh(m_vSites[(i + 1) * m_iSizeX + j]);
			m_vSites[i * m_iSizeX + j]->setNeighPosition(m_vSites[(i + 1) * m_iSizeX + j], Site::NORTH);

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
			m_vSites[j]->setNeighPosition(m_vSites[j + m_iSizeX], Site::NORTH);

			m_vSites[j]->setNeigh(m_vSites[iThirdCorner + j]);
			m_vSites[j]->setNeighPosition(m_vSites[iThirdCorner + j], Site::SOUTH);
		}
		else if (j == iFirstCorner)
		{
			m_vSites[j]->setNeigh(m_vSites[iSecondCorner]);
			m_vSites[j]->setNeighPosition(m_vSites[iSecondCorner], Site::WEST);

			m_vSites[j]->setNeigh(m_vSites[1]);
			m_vSites[j]->setNeighPosition(m_vSites[1], Site::EAST);

			m_vSites[j]->setNeigh(m_vSites[iSecondCorner + 1]);
			m_vSites[j]->setNeighPosition(m_vSites[iSecondCorner + 1], Site::NORTH);

			m_vSites[j]->setNeigh(m_vSites[iThirdCorner]);
			m_vSites[j]->setNeighPosition(m_vSites[iThirdCorner], Site::SOUTH);
		}
		else if (j == iSecondCorner)
		{
			m_vSites[j]->setNeigh(m_vSites[j - 1]);
			m_vSites[j]->setNeighPosition(m_vSites[j - 1], Site::WEST);

			m_vSites[j]->setNeigh(m_vSites[0]);
			m_vSites[j]->setNeighPosition(m_vSites[0], Site::EAST);

			m_vSites[j]->setNeigh(m_vSites[2 * m_iSizeX - 1]);
			m_vSites[j]->setNeighPosition(m_vSites[2 * m_iSizeX - 1], Site::NORTH);

			m_vSites[j]->setNeigh(m_vSites[iForthCorner]);
			m_vSites[j]->setNeighPosition(m_vSites[iForthCorner], Site::SOUTH);
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
			m_vSites[j]->setNeighPosition(m_vSites[iFirstCorner + iPos], Site::NORTH);

			m_vSites[j]->setNeigh(m_vSites[j - m_iSizeX]);
			m_vSites[j]->setNeighPosition(m_vSites[j - m_iSizeX], Site::SOUTH);
			iPos++;
		}
		else if (j == iThirdCorner)
		{
			m_vSites[j]->setNeigh(m_vSites[iForthCorner]);
			m_vSites[j]->setNeighPosition(m_vSites[iForthCorner], Site::WEST);

			m_vSites[j]->setNeigh(m_vSites[iThirdCorner + 1]);
			m_vSites[j]->setNeighPosition(m_vSites[iThirdCorner + 1], Site::EAST);

			m_vSites[j]->setNeigh(m_vSites[iFirstCorner]);
			m_vSites[j]->setNeighPosition(m_vSites[iFirstCorner], Site::NORTH);

			m_vSites[j]->setNeigh(m_vSites[iThirdCorner - m_iSizeX]);
			m_vSites[j]->setNeighPosition(m_vSites[iThirdCorner - m_iSizeX], Site::SOUTH);
		}
		else if (j == iForthCorner)
		{
			m_vSites[j]->setNeigh(m_vSites[iForthCorner - 1]);
			m_vSites[j]->setNeighPosition(m_vSites[iForthCorner - 1], Site::WEST);

			m_vSites[j]->setNeigh(m_vSites[iThirdCorner]);
			m_vSites[j]->setNeighPosition(m_vSites[iThirdCorner], Site::EAST);

			m_vSites[j]->setNeigh(m_vSites[iSecondCorner]);
			m_vSites[j]->setNeighPosition(m_vSites[iSecondCorner], Site::NORTH);

			m_vSites[j]->setNeigh(m_vSites[iThirdCorner - 1]);
			m_vSites[j]->setNeighPosition(m_vSites[iThirdCorner - 1], Site::SOUTH);
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
		m_vSites[j]->setNeighPosition(m_vSites[j + m_iSizeX], Site::NORTH);

		m_vSites[j]->setNeigh(m_vSites[j - m_iSizeX]);
		m_vSites[j]->setNeighPosition(m_vSites[j - m_iSizeX], Site::SOUTH);
	}

	/* Last column */
	for (int j = iSecondCorner + m_iSizeX; j < iForthCorner; j += m_iSizeX)
	{
		m_vSites[j]->setNeigh(m_vSites[j - 1]);
		m_vSites[j]->setNeighPosition(m_vSites[j - 1], Site::WEST);

		m_vSites[j]->setNeigh(m_vSites[j - m_iSizeX + 1]);
		m_vSites[j]->setNeighPosition(m_vSites[j - m_iSizeX + 1], Site::EAST);

		m_vSites[j]->setNeigh(m_vSites[j + m_iSizeX]);
		m_vSites[j]->setNeighPosition(m_vSites[j + m_iSizeX], Site::NORTH);

		m_vSites[j]->setNeigh(m_vSites[j - m_iSizeX]);
		m_vSites[j]->setNeighPosition(m_vSites[j - m_iSizeX], Site::SOUTH);
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

Site *BCC::getSite(int id) { return m_vSites[id]; }

void BCC::check()
{
	int k = 0;

	cout << "Checking lattice..." << endl;

	int test = 2;
	cout << test << ": ";
	cout << "W:" << getSite(test)->getNeighPosition(Site::WEST)->getID() << " ";
	cout << "Wu:" << getSite(test)->getNeighPosition(Site::WEST_UP)->getID() << " ";
	cout << "WD:" << getSite(test)->getNeighPosition(Site::WEST_DOWN)->getID() << " ";
	cout << "E:" << getSite(test)->getNeighPosition(Site::EAST)->getID() << " ";
	cout << "EU:" << getSite(test)->getNeighPosition(Site::EAST_UP)->getID() << " ";
	cout << "ED:" << getSite(test)->getNeighPosition(Site::EAST_DOWN)->getID() << " ";
	cout << "N:" << getSite(test)->getNeighPosition(Site::NORTH)->getID() << " ";
	cout << "S:" << getSite(test)->getNeighPosition(Site::SOUTH)->getID() << endl;

	cout << "Activation: " << endl;

	cout << "N:" << getSite(test)->getActivationSite(Site::ACTV_NORTH)->getID() << " ";
	cout << "S:" << getSite(test)->getActivationSite(Site::ACTV_SOUTH)->getID() << " ";
	cout << "E:" << getSite(test)->getActivationSite(Site::ACTV_EAST)->getID() << " ";
	cout << "W:" << getSite(test)->getActivationSite(Site::ACTV_WEST)->getID() << endl;
}

void BCC::updateNeighbours(Site* site)
{
	int siteHeight = site->getHeight();

  	int totalNeigh = 0;
  	site->m_clearNeighbourList();

	// Check NESW sites, see if the heights are the same. If same, add to list of neighbours.
    bool isActiveEAST = false;
    isActiveEAST = siteHeight <= site->getNeighPosition(Site::EAST)->getHeight();
    if (isActiveEAST)
    {
      site->m_addSite(site->getNeighPosition(Site::EAST));
      totalNeigh++;
    }

    bool isActiveWEST = false;
    isActiveWEST = siteHeight <= site->getNeighPosition(Site::WEST)->getHeight();
    if (isActiveWEST)
    {
      site->m_addSite(site->getNeighPosition(Site::WEST));
      totalNeigh++;
    }

    bool isActiveNORTH = false;
    isActiveNORTH = siteHeight <= site->getNeighPosition(Site::NORTH)->getHeight();
    if (isActiveNORTH)
    {
      site->m_addSite(site->getNeighPosition(Site::NORTH));
      totalNeigh++;
    }

    bool isActiveSOUTH = false;
    isActiveSOUTH = siteHeight <= site->getNeighPosition(Site::SOUTH)->getHeight();
    if (isActiveSOUTH)
    {
      site->m_addSite(site->getNeighPosition(Site::SOUTH));
      totalNeigh++;
    }
}
