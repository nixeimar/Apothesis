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
typedef rapidjson::Value Value;

Read::Read(Apothesis* apothesis):Pointers( apothesis),
                 m_sLatticeType("NONE"),
                 m_sProcess("process"),
                 m_sLattice("lattice"),
                 m_sTemperature("temperature"),
                 m_sPressure("pressure"),
                 m_sIterations("num_iterations"),
                 m_sCommentLine("#"),
                 m_input(readInputFile("input.kmc")),
                 m_dimensions(3)
  {
    //Initialize the map for the lattice
    m_LatticeType[ "NONE" ] = Lattice::NONE;
    m_LatticeType[ "BCC" ] = Lattice::BCC;
    m_LatticeType[ "FCC" ] = Lattice::FCC;

    
    string lattice_type = m_input["Lattice"]["Type"].GetString();
    std::cout<<"lattice_type "<< lattice_type << std::endl;
    
    m_lattice->setType(lattice_type);

    // Initialize the vector holding dimensions of the lattice
    vector<int> latticeDimensions;

    // Read the elements in the input file
    //auto dimensions = m_input["lattice"]["dims"].GetObject();

    for (int i = 0; i < m_dimensions; i++)
    {
      //latticeDimensions.push_back(latticeDimensions[i]);
      latticeDimensions.push_back(10);
    }

    std::cout<<"Setting lattice dimensions " << std::endl;
    m_lattice->setX(latticeDimensions[0]);
    m_lattice->setY(latticeDimensions[1]);
    m_lattice->setInitialHeight(latticeDimensions[2]);
    
    std::cout<<"Setting Iterations " << std::endl;
    // Set the iterations, temperature, and pressure
    m_parameters->setIterations(m_input["Iterations"].GetInt());

    std::cout<<"Setting Temperature " << std::endl;
    m_parameters->setTemperature(m_input["Temperature"].GetDouble());

    std::cout<<"Setting Pressure " << std::endl;
    m_parameters->setPressure(m_input["Pressure"].GetDouble());

    // Storing all processes into apothesis class
    Value& process = m_input["process"];
    for (Value::ConstMemberIterator itr = process.MemberBegin(); itr != process.MemberEnd(); ++itr)
    {
      apothesis->addProcess(itr->name.GetString());
    }

    // Storing all processes into apothesis class
    Value& speciesName = m_input["Species"];
    for (Value::ConstMemberIterator itr = process.MemberBegin(); itr != process.MemberEnd(); ++itr)
    {
      m_speciesName.push_back(itr->name.GetString());
    }

    for (vector<string>::const_iterator itr = m_speciesName.begin(); itr != m_speciesName.end(); ++itr)
    {
      const char* name = (*itr).c_str();
      double mw = m_input["Species"][name]["mw"].GetDouble();
      m_MWs.push_back(mw);
    }
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

Document& Read::getDoc()
{
  return m_input; 
}

vector<string> Read::getSpeciesNames()
{
  return m_speciesName;
}

vector<double> Read::getMWs()
{
  return m_MWs;
}