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

#include "diffusion.h"
#include "register.cpp"

namespace MicroProcesses{

REGISTER_PROCESS_IMPL( Diffusion)

Diffusion::Diffusion(){;}

Diffusion::~Diffusion(){;}

void Diffusion::init(){;}

//void Diffusion::setName( string name ){;}

string Diffusion::getName(){;}

void Diffusion::activeSites( Lattice* ){;}

void Diffusion::selectSite(){;}

void Diffusion::perform(){;}

double Diffusion::getProbability(){;}

list<Site* > Diffusion::getActiveList(){;}

void Diffusion::setProcessMap(map<Process *, list<Site *> *> *){;}

void Diffusion::test(){;}

}
