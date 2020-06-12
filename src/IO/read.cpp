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

#include "read.h"

typedef rapidjson::Document Document;

Read::Read(Apothesis* apothesis):Pointers( apothesis),
                 m_sLatticeType("NONE"),
                 m_sProcess("process"),
                 m_sLattice("lattice"),
                 m_sTemperature("temperature"),
                 m_sPressure("pressure"),
                 m_sIterations("num_iterations"),
                 m_sCommentLine("#")
  {

  }

Read::~Read(){}

void Read::init( int argc, char* argv[] )
{

}


Document Read::readInputFile(string filename)
{

}

