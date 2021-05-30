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

#include "species.h"

Species::Species
(
    string name,
    double mw,
    double stoicCoeff,
    int id
)
:
m_name(name),
m_mw(mw),
m_stoicCoeff(stoicCoeff),
m_id(id)
{  
    ;
}

Species::Species
(
    string name,
    double mw,
    int id
)
:
m_name(name),
m_mw(mw),
m_id(id)
{  
    ;
}

// Deconstructor for species class
Species::~Species(){}

// Returns name of the species
string Species::getName()
{
    return m_name;
}

// Returns molecular weight of species
double Species::getMW()
{
    return m_mw;
}

// Return 
double Species::getStoicCoeff()
{
    return m_stoicCoeff;
}

const int Species::getId()
{
    return m_id;
}
