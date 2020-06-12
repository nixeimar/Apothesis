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

// typedef the relevant rapidjson classes
typedef rapidjson::Document Document;
typedef rapidjson::FileReadStream FileReadStream;

Read::Read(Apothesis* apothesis):Pointers( apothesis),
                 m_sLatticeType("NONE"),
                 m_sProcess("process"),
                 m_sLattice("lattice"),
                 m_sTemperature("temperature"),
                 m_sPressure("pressure"),
                 m_sIterations("num_iterations"),
                 m_sCommentLine("#"),
                 m_input(readInputFile("input.kmc"))
  {
    //Initialize the map for the lattice
    m_mLatticeType[ "NONE" ] = Lattice::NONE;
    m_mLatticeType[ "BCC" ] = Lattice::BCC;
    m_mLatticeType[ "FCC" ] = Lattice::FCC;

    string lattice_type = m_input["lattice"]["type"].GetString();

  }

Read::~Read(){}

Document Read::readInputFile(string filename)
{
  FILE* fp = fopen(filename.c_str(), "r");

  // Throws an error if the file cannot be opened
  if (fp == NULL)
  {
    throw std::runtime_error(std::strerror(errno));
  }
  char readBuffer[65536];
  FileReadStream is(fp, readBuffer, sizeof(readBuffer));

  Document doc;
  doc.ParseStream(is);

  fclose(fp);

  return doc;  
}