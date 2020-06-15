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
//typedef rapidjson::GenericValue<rapidjson::Encoding, rapidjson::Allocator> ValueType;
//typedef rapidjson::GenericObject<false, ValueType> Object;

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

    // How to show diffusion?
    
    vector <double> diffusionParameters;
    //auto process = m_input["Process"].GetObject();

    std::cout<<"Setting Diffusion parameter 1 " << std::endl;

    double param1 = m_input["Process"]["Diffusion"]["param_1"].GetDouble();
    double param2 = m_input["Process"]["Diffusion"]["param_2"].GetDouble();
    
    diffusionParameters.push_back(param1);
    diffusionParameters.push_back(param2);

    std::cout<<"Setting Diffusion parameter 2 " << std::endl;
    diffusionParameters.push_back(m_input["Process"]["Diffusion"]["param_2"].GetDouble());

    std::cout<<"Setting Diffusion Parameters " << std::endl;
    m_parameters->setProcess( "Diffusion", diffusionParameters );
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