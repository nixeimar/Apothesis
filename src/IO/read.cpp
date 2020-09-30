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

Read::Read(Apothesis *apothesis)
    : Pointers(apothesis),
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
  m_LatticeType["NONE"] = Lattice::NONE;
  m_LatticeType["BCC"] = Lattice::BCC;
  m_LatticeType["FCC"] = Lattice::FCC;

  string lattice_type = m_input["Lattice"]["Type"].GetString();
  std::cout << "lattice_type " << lattice_type << std::endl;

  Value &pLattice = m_input["Lattice"];

  bool steps = false;
  vector<int> vStepInfo = {0, 0, 0};
  if (pLattice.HasMember("Step"))
  {
    steps = true;
    const Value &stepInfo = m_input["Lattice"]["Step"];
    assert(stepInfo.IsArray());
    vStepInfo.clear();
    vStepInfo = {stepInfo[0].GetInt(), stepInfo[1].GetInt(), stepInfo[2].GetInt()};
  }

  switch (m_LatticeType[lattice_type])
  {
  case Lattice::FCC:
  {
    // TODO I think STEPS parameter needs to be incorporated along with the lattice
    FCC *lattice = new FCC(apothesis);
    apothesis->pLattice = lattice;
    m_lattice->setType("FCC");
    break;
  }
  case Lattice::BCC:
  {
    if (steps)
    {
      BCC *lattice = new BCC(apothesis, true, vStepInfo);
      apothesis->pLattice = lattice;
    }
    else
    {
      BCC *lattice = new BCC(apothesis);
      apothesis->pLattice = lattice;
    }
    m_lattice->setType("BCC");

    break;
  }
  //case Lattice::FCC: apothesis->pLattice = new BCC(apothesis);
  //TODO: In case of default, create error
  default:
  {
    ;
  }
  }

  m_lattice = apothesis->pLattice;
  m_lattice->setType(lattice_type);

  // Initialize the vector holding dimensions of the lattice
  vector<int> latticeDimensions;

  // Read the elements in the input file
  Value &dimensions = m_input["Lattice"]["dims"];
  apothesis->logSuccessfulRead(dimensions.IsArray(), "Lattice dimensions");
  for (int i = 0; i < m_dimensions; i++)
  {
    if (!dimensions[i].IsInt())
      apothesis->pErrorHandler->error_simple_msg("Lattice dimension is not a number");
    latticeDimensions.push_back((int)dimensions[i].GetDouble());
  }

  //    std::cout<<"Setting lattice dimensions " << std::endl;
  m_lattice->setX(latticeDimensions[0]);
  m_lattice->setY(latticeDimensions[1]);
  m_lattice->setInitialHeight(latticeDimensions[2]);

  // Set the iterations, temperature, and pressure
  apothesis->logSuccessfulRead(m_input.HasMember("Iterations"), "Iterations");
  m_parameters->setIterations(m_input["Iterations"].GetInt());

  apothesis->logSuccessfulRead(m_input.HasMember("Temperature"), "Temperature");
  m_parameters->setTemperature(m_input["Temperature"].GetDouble());

  apothesis->logSuccessfulRead(m_input.HasMember("Pressure"), "Pressure");
  m_parameters->setPressure(m_input["Pressure"].GetDouble());

  // Storing all processes names into apothesis class
  apothesis->logSuccessfulRead(m_input.HasMember("Process"), "Process");
  Value &process = m_input["Process"];

  // Add processes into a vector, stored under apothesis
  for (Value::ConstMemberIterator itr = process.MemberBegin(); itr != process.MemberEnd(); ++itr)
  {
    apothesis->addProcess(itr->name.GetString());
  }

  // Storing all species into apothesis class
  apothesis->logSuccessfulRead(m_input.HasMember("Species"), "Species");
  Value &speciesName = m_input["Species"];

  // Push names of existing species
  for (Value::ConstMemberIterator itr = speciesName.MemberBegin(); itr != speciesName.MemberEnd(); ++itr)
  {
    if (itr->name.IsString())
      m_speciesName.push_back(itr->name.GetString());
  }

  // Read species and corresponding molecular weights into member vars within apothesis
  for (vector<string>::const_iterator itr = m_speciesName.begin(); itr != m_speciesName.end(); ++itr)
  {
    const char *name = (*itr).c_str();
    apothesis->logSuccessfulRead(speciesName.HasMember(name), name);

    // TODO: Include a map with default values in case nothing is specified
    double mw = m_input["Species"][name]["mw"].GetDouble();
    m_MWs.push_back(mw);
  }

  // Read special configuration settings
  // Storing all processes names into apothesis class

  try
  {
    Value &configPtr = m_input["config"];
    apothesis->logSuccessfulRead(configPtr.HasMember("debug"), "debug");
    string debugMode = configPtr["debug"].GetString();
    if (!debugMode.compare("On") || !debugMode.compare("on") || !debugMode.compare("True") || !debugMode.compare("true") || !debugMode.compare("debug"))
    {
      apothesis->setDebugMode(true);
    }
  }
  catch (const std::exception &e)
  {
    // do nothing
  }

  // Add processes into a vector, stored under apothesis
  for (Value::ConstMemberIterator itr = process.MemberBegin(); itr != process.MemberEnd(); ++itr)
  {
    apothesis->addProcess(itr->name.GetString());
  }
}

Read::~Read() {}

Document Read::readInputFile(string filename)
{
  FILE *fp = fopen(filename.c_str(), "r");

  // Throws an error if the file cannot be opened
  if (fp == NULL)
  {
    throw std::runtime_error(std::strerror(errno));
  }
  else
  {
    std::cout << filename << " successfully opened." << std::endl;
  }

  char readBuffer[65536];
  FileReadStream is(fp, readBuffer, sizeof(readBuffer));

  Document doc;
  doc.ParseStream(is);

  doc.IsObject() ? std::cout << "Successfully parsed input file" << std::endl
                 : std::cout << "Error in format of input file" << std::endl;

  fclose(fp);

  return doc;
}

Document &Read::getDoc()
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