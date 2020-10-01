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

#include "pointers.h"
#include "apothesis.h"
#include "lattice.h"
#include "FCC.h"
#include "io.h"
#include "read.h"
#include "errorhandler.h"
#include "parameters.h"
#include "process.h"
#include "species.h"
#include "string.h"
#include "adsorption.h"
#include "desorption.h"
#include "diffusion.h"
#include "SurfaceReaction.h"
#include <numeric>

using namespace MicroProcesses;

typedef rapidjson::Document Document;
typedef rapidjson::Value Value;
typedef rapidjson::SizeType SizeType;

//using namespace Utils;

Apothesis::Apothesis(int argc, char *argv[])
    : pLattice(0),
      pRead(0),
      m_debugMode(false)
{
  m_iArgc = argc;
  m_vcArgv = argv;

  pParameters = new Utils::Parameters(this);

  /* This must be constructed before the input */

  // Create input instance
  pIO = new IO(this);

  pRead = new Read(this);

  vector<string> pName = pRead->getSpeciesNames();

  // Build the lattice. This should always follow the read input

  std::cout << "Building the lattice" << std::endl;
  pLattice->build();
  std::cout << "Finished building the lattice" << std::endl;
}

Apothesis::~Apothesis()
{
  delete pIO;
  delete pRead;
  delete pLattice;

  // Delete the processes created by the factory method
  for (vector<Process *>::iterator it = m_vProcesses.begin();
       it != m_vProcesses.end(); it++)
  {
    delete *it;
  }
}

void Apothesis::init()
{
  cout << "Opening output file" << endl;
  /// The output file name will come from the user and will have the extenstion .log
  /// This would come as a parameter from the user from the args (also the input).
  /// Now both are hard copied.

  if (!pIO->outputOpen())
    pIO->openOutputFile("Output");

  // Processes in this case
  vector<string> pProc = m_processes;

  cout << pProc[0] << endl;

  Document &doc = pRead->getDoc();

  pIO->writeLogOutput("Initializing instances of species");
  if (std::find(pProc.begin(), pProc.end(), "Reaction") == pProc.end())
  {
    vector<double> mws = pRead->getMWs();
    vector<string> names = pRead->getSpeciesNames();

    for (int i = 0; i < mws.size(); ++i)
    {
      Species *s = new Species(names[i], mws[i]);
      m_species[names[i]] = s;
    }
  }

  // Initializing interactions between species
  pIO->writeLogOutput("Reading interactions between species");
  Value &speciesName = doc["Species"];
  for (Value::ConstMemberIterator itr = speciesName.MemberBegin(); itr != speciesName.MemberEnd(); ++itr)
  {
    const char *name = itr->name.GetString();
    Value &singleSpecies = speciesName[name];
    if (singleSpecies.HasMember("Interactions"))
    {
      Value &interactions = singleSpecies["Interactions"];

      //TODO better logging function
      logSuccessfulRead(interactions.IsArray(), "Interactions");

      int numInters = interactions.Size();
      for (int i = 0; i < numInters; i++)
      {
        m_interactions.push_back(make_tuple(name, interactions[i].GetString()));
      }
    }
  }

  pIO->writeLogOutput("Initializing processes");

  if (std::find(pProc.begin(), pProc.end(), "Adsorption") != pProc.end())
  {
    pIO->writeLogOutput("Initializing Adsorption");

    // Read parameters for Adsorption
    Value &specie = doc["Process"]["Adsorption"]["Species"];
    Value &stick = doc["Process"]["Adsorption"]["Sticking"];
    Value &mFrac = doc["Process"]["Adsorption"]["MassFraction"];

    // Verify presence of each parameter in input file
    logSuccessfulRead(specie.IsArray(), "Adsorption species");
    logSuccessfulRead(stick.IsArray(), "Adsorption sticking coefficients");
    logSuccessfulRead(mFrac.IsArray(), "Adsorption mass fraction");

    // Initialize vectors
    vector<string> species;
    vector<double> sticking;
    vector<double> massFraction;

    // Sum of mass fraction. Later used to normalize.
    double sum = 0;

    for (SizeType i = 0; i < specie.Size(); i++)
    {
      // Output possible errors
      if (!specie[i].IsString())
        pErrorHandler->error_simple_msg("Species format is not a string");
      if (!stick[i].IsNumber())
        pErrorHandler->error_simple_msg("Sticking coefficient format is not a double");
      if (!mFrac[i].IsNumber())
        pErrorHandler->error_simple_msg("Mass fraction format is not a double");

      // Push values to corresponding vectors
      species.push_back(specie[i].GetString());
      sticking.push_back(stick[i].GetDouble());
      massFraction.push_back(mFrac[i].GetDouble());
      sum += massFraction[i];
    }

    // Normalize the values of the mass fraction
    for (vector<double>::iterator itr = massFraction.begin(); itr != massFraction.end(); ++itr)
    {
      *itr = *itr / sum;
    }

    for (int i = 0; i < species.size(); ++i)
    {
      Adsorption *a = new Adsorption(this, species[i], m_species[species[i]], sticking[i], massFraction[i]);

      // Keep two separate vectors: one for all processes, one for adsorption processes only
      m_vAdsorption.push_back(a);
      m_vProcesses.push_back(a);
    }

    // Resolving Conflicts

    pIO->writeLogOutput("...Done initializing Adsorption process.");
  }
  if (std::find(pProc.begin(), pProc.end(), "Desorption") != pProc.end())
  {
    pIO->writeLogOutput("Initializing Desorption");

    // Read parameters for Desorption
    Value &vSpecie = doc["Process"]["Desorption"]["Species"];
    Value &vEnergy = doc["Process"]["Desorption"]["Energy"];
    Value &vFreq = doc["Process"]["Desorption"]["Frequency"];

    // Verify presence of each parameter in input file
    logSuccessfulRead(vSpecie.IsArray(), "Desorption species");
    logSuccessfulRead(vEnergy.IsArray(), "Desorption energy");
    logSuccessfulRead(vFreq.IsArray(), "Desorption frequency");

    // Initialize vectors
    vector<string> species;
    vector<double> energy;
    vector<double> frequency;

    for (SizeType i = 0; i < vSpecie.Size(); i++)
    {
      // Output possible errors
      if (!vSpecie[i].IsString())
        pErrorHandler->error_simple_msg("Species format is not a string");
      if (!vEnergy[i].IsNumber())
        pErrorHandler->error_simple_msg("Desorption energy format is not a number");
      if (!vFreq[i].IsNumber())
        pErrorHandler->error_simple_msg("Desorption frequency format is not a number");

      // Push values to corresponding vectors
      species.push_back(vSpecie[i].GetString());
      energy.push_back(vEnergy[i].GetDouble());
      frequency.push_back(vFreq[i].GetDouble());
    }

    // Add process to m_vProcesses
    for (int i = 0; i < species.size(); ++i)
    {
      // Call function to find adsorption class
      Adsorption *pAdsorption = findAdsorption(species[i]);

      // Create new instance of desorption class
      Desorption *d = new Desorption(this, species[i], m_species[species[i]], energy[i], frequency[i]);

      // Check if associated adsorption class is a valid (non-null) pointer. Associate desorption class with appropriate adsorption and vice versa
      if (pAdsorption)
      {
        d->setAdsorptionPointer(pAdsorption);
        pAdsorption->setDesorptionPointer(d);
        pAdsorption->setDesorption(true);
      }

      // if not, create associations in adsorption and desorption

      // else, output warning

      m_vDesorption.push_back(d);
      m_vProcesses.push_back(d);
    }
    pIO->writeLogOutput("...Done initializing desorption process.");
  }
  if (std::find(pProc.begin(), pProc.end(), "Diffusion") != pProc.end())
  {
    pIO->writeLogOutput("Initializing Diffusion");

    // Read parameters for Diffusion
    Value &vSpecie = doc["Process"]["Diffusion"]["Species"];
    Value &vEnergy = doc["Process"]["Diffusion"]["Energy"];
    Value &vFreq = doc["Process"]["Diffusion"]["Frequency"];

    // Verify presence of each parameter in input file
    logSuccessfulRead(vSpecie.IsArray(), "Diffusion species");
    logSuccessfulRead(vEnergy.IsArray(), "Diffusion energy");
    logSuccessfulRead(vFreq.IsArray(), "Diffusion frequency");

    // Initialize vectors
    vector<string> species;
    vector<double> energy;
    vector<double> frequency;

    for (SizeType i = 0; i < vSpecie.Size(); i++)
    {
      // Output possible errors
      if (!vSpecie[i].IsString())
        pErrorHandler->error_simple_msg("Species format is not a string");
      if (!vEnergy[i].IsNumber())
        pErrorHandler->error_simple_msg("Diffusion energy format is not a number");
      if (!vFreq[i].IsNumber())
        pErrorHandler->error_simple_msg("Diffusion frequency format is not a number");

      // Push values to corresponding vectors
      species.push_back(vSpecie[i].GetString());
      energy.push_back(vEnergy[i].GetDouble());
      frequency.push_back(vFreq[i].GetDouble());
    }

    // Add process to m_vProcesses
    for (int i = 0; i < species.size(); ++i)
    {
      // Call function to find adsorption class
      Adsorption *pAdsorption = findAdsorption(species[i]);
      Desorption *pDesorption = findDesorption(species[i]);

      Diffusion *diff = new Diffusion(this, species[i], energy[i], frequency[i]);
      diff->setAdsorptionPointer(pAdsorption);
      diff->setDesorptionPointer(pDesorption);

      pAdsorption->setDiffusion(true);
      pAdsorption->setDiffusionPointer(diff);

      pDesorption->setDiffusion(true);
      pDesorption->setDiffusionPointer(diff);

      m_vProcesses.push_back(diff);
    }
    pIO->writeLogOutput("...Done initializing diffusion process.");
  }
  if (std::find(pProc.begin(), pProc.end(), "Reaction") != pProc.end())
  {
    pIO->writeLogOutput("Initializing Reaction");

    // Read parameters for Reaction
    Value &pRxn = doc["Process"]["Reaction"];
    Value &vSpecie = doc["Process"]["Reaction"]["Species"];
    Value &vStoich = doc["Process"]["Reaction"]["Stoichiometry"];
    Value &vEnergy = doc["Process"]["Reaction"]["Energy"];
    Value &vPreExp = doc["Process"]["Reaction"]["PreExp"];

    // Verify presence of each parameter in input file
    logSuccessfulRead(vSpecie.IsArray(), "Reaction species");
    logSuccessfulRead(vStoich.IsArray(), "Reaction Stoichiometry");
    logSuccessfulRead(vEnergy.IsNumber(), "Enthalpy of reaction");
    logSuccessfulRead(vPreExp.IsNumber(), "Pre-exponential factor for the reaction");

    // Initialize vectors
    vector<Species *> species;
    vector<double> stoichiometry;
    double energy;
    double preexp;

    auto mws = pRead->getMWs();

    // Loop through species
    for (SizeType i = 0; i < vSpecie.Size(); i++)
    {
      // Output possible errors
      if (!vSpecie[i].IsString())
        pErrorHandler->error_simple_msg("Species format is not a string");
      if (!vStoich[i].IsNumber())
        pErrorHandler->error_simple_msg("Diffusion energy format is not a number");

      // Push values to corresponding vectors
      Species *s = getSpecies(vSpecie[i].GetString());
      species.push_back(s);
      stoichiometry.push_back(vStoich[i].GetDouble());

      string name = vSpecie[i].GetString();
      // If the species is not already defined, push new member onto maps
      if (m_species[name] == NULL)
      {
        Species *s = new Species(name, mws[i], vStoich[i].GetDouble());
        m_species[name] = s;
      }
    }

    // Store value for energy and pre-exponential factor
    energy = vEnergy.GetDouble();
    preexp = vPreExp.GetDouble();

    // Check mass balance on reaction
    double cumulativemass = 0;
    for (map<string, Species *>::iterator itr = m_species.begin(); itr != m_species.end(); ++itr)
    {
      Species *s = itr->second;
      cumulativemass += s->getMW() * s->getStoicCoeff();
    }

    if (abs(cumulativemass) > 1e-10)
    {
      cout << "Warning! Mass balance of Reaction is not balanced" << endl;
    }

    // Read immobilization variable
    bool immobilized = true;
    if (pRxn.HasMember("Immobilize"))
    {
      immobilized = pRxn["Immobilize"].GetBool();
    }

    SurfaceReaction *s = new SurfaceReaction(this, species, stoichiometry, energy, preexp, immobilized);
    m_vProcesses.push_back(s);
    m_vSurfaceReaction.push_back(s);
    pIO->writeLogOutput("...Done initializing reaction.");
  }

  // Initialize interactions between adsorption species and classes
  vector<tuple<string, string>>::iterator itr = m_interactions.begin();
  // For each pair of interactions, add the possible interaction (2-way) within the appropriate pointers
  for (itr; itr != m_interactions.end(); ++itr)
  {
    tuple<string, string> temp = *itr;
    string spec1 = get<0>(temp);
    string spec2 = get<1>(temp);

    Species *s1 = m_species[spec1];
    Species *s2 = m_species[spec2];

    //TODO: Add error catch statements and more comments here
    Adsorption *pAdsorption1 = findAdsorption(spec1);
    pAdsorption1->addInteraction(s2);

    Adsorption *pAdsorption2 = findAdsorption(spec2);
    pAdsorption2->addInteraction(s1);
  }
  vector<Adsorption*> adsorptionpointers = getAdsorptionPointers();
  for (vector<Adsorption *>::iterator iter = adsorptionpointers.begin(); iter != adsorptionpointers.end(); ++iter)
  {
    Adsorption *pAds = *iter;
    vector<Species *> possibleInteractions = pAds->getInteractions();
    bool found = false;
    for (int i = 0; i < possibleInteractions.size(); ++i)
    {
      cout<< possibleInteractions[i]->getName() << ", ";
    }
  }
  cout<<endl;
  cout<<endl;
  /// First the processes that participate in the simulation

  /// that were read from the file input and the I/O functionality
  //m_vProcesses[0]->setInstance( this );
  for (vector<Process *>::iterator itr = m_vProcesses.begin(); itr != m_vProcesses.end(); ++itr)
  {
    Process *p = *itr;
    p->activeSites(pLattice);
  }
}

void Apothesis::exec()
{
  ///Perform the number of KMC steps read from the input.
  int iterations = pParameters->getIterations();

  if (iterations == 0)
  {
    pErrorHandler->error_simple_msg("Zero iterations found.");
    EXIT;
  }
  else
  {
    pIO->writeLogOutput("Running " + to_string(iterations) + " iterations");
  }

  /// Get list of possible processes
  for (int i = 0; i < iterations; ++i)
  {
    /// Print to output
    pIO->writeLogOutput("Time step: " + to_string(i));

    //TODO: If (debug) print out the id of the lattice site
    vector<Process *> processes = m_vProcesses;
    /// Find probability of each process
    vector<double> probabilities = calculateProbabilities(m_vProcesses);
    /// Pick random number with 3 digits
    double random = (double)rand() / RAND_MAX;

    /// Pick Process
    Process *p = pickProcess(probabilities, random, processes);

    // Site should be picked here
    p->selectSite();

    /// Perform process on that site
    p->perform();

    // The frequency that the various information are written in the file
    // must befined by the user. Fix it ...
    // The user should also check if the messages are written on the terminal or not.
    pIO->writeLogOutput(p->getName() + " ");
    pIO->writeLatticeHeights();

    if (i == 0)
    {
      cout << m_vProcesses[0]->getName() << " process is being performed..." << endl;
    }
    // TODO: check if all processes are 0
  }
}

void Apothesis::addProcess(string process)
{
  m_processes.push_back(process);
}

void Apothesis::logSuccessfulRead(bool read, string parameter)
{
  if (!pIO->outputOpen())
  {
    pIO->openOutputFile("Output");
  }

  read ? pIO->writeLogOutput("Reading " + parameter)
       : pErrorHandler->error_simple_msg("No " + parameter + " found in input file");
}

map<string, Species *> Apothesis::getAllSpecies()
{
  return m_species;
}

Species *Apothesis::getSpecies(string species)
{
  return m_species[species];
}

vector<double> Apothesis::calculateProbabilities(vector<Process *> pProcesses)
{
  /// Calculate probabilities for each process
  int numProcesses = pProcesses.size();

  /// Calculate probability of each
  vector<double> probability;
  double total = 0;
  vector<Process *>::iterator itr = pProcesses.begin();

  // Find the probability for each process. Push onto prob.
  for (; itr != pProcesses.end(); ++itr)
  {
    Process *process = *itr;
    double prob = process->getProbability();
    probability.push_back(prob + total);
    total += prob;
  }

  // Normalize all values
  for (vector<double>::iterator itr = probability.begin(); itr != probability.end(); ++itr)
  {
    *itr = *itr / total;
  }

  return probability;
}

// May be possible to delete (if pProcesses is a vector, can simply access by index)
Process *Apothesis::getProcessAt(int index, vector<Process *> pProcesses)
{
  vector<Process *>::iterator it = pProcesses.begin();
  std::advance(it, index);
  return *it;
}

Process *Apothesis::pickProcess(vector<double> probabilities, double random, vector<Process *> pProcesses)
{
  for (int index = 0; index < probabilities.size(); ++index)
  {
    if (random < probabilities[index])
    {
      return pProcesses[index];
    }
  }
  return pProcesses[probabilities.size() - 1];
}

Adsorption *Apothesis::findAdsorption(string species)
{
  for (vector<Adsorption *>::iterator itr = m_vAdsorption.begin(); itr != m_vAdsorption.end(); ++itr)
  {
    Adsorption *a = *itr;
    if (!species.compare(a->getSpeciesName()))
    {
      return a;
    }
    // find name of adsorption species
  }
  cout << "Warning! Could not find instance of Adsorption class for desorbed species " << species << endl;
}

Desorption *Apothesis::findDesorption(string species)
{
  for (vector<Desorption *>::iterator itr = m_vDesorption.begin(); itr != m_vDesorption.end(); ++itr)
  {
    Desorption *d = *itr;
    if (!species.compare(d->getSpeciesName()))
    {
      return d;
    }
    // find name of adsorption species
  }
  cout << "Warning! Could not find instance of Adsorption class for desorbed species " << species << endl;
}

IO *Apothesis::getIOPointer()
{
  return pIO;
}

void Apothesis::setDebugMode(bool ifDebug)
{
  m_debugMode = ifDebug;
}

bool Apothesis::getDebugMode()
{
  return m_debugMode;
}

void Apothesis::setLatticePointer(Lattice *lattice)
{
  pLattice = lattice;
}

vector<Adsorption *> Apothesis::getAdsorptionPointers()
{
  return m_vAdsorption;
}

vector<SurfaceReaction *> Apothesis::getReactionPointers()
{
  return m_vSurfaceReaction;
}
