#include "reader.h"

void Reader::parseFile(){

    openInputFile(m_inputPath);
    vector<string> fLines=inputFileLines();

    for(int i=0; i < fLines.size(); ++i) {
        string currentLine;
        currentLine=fLines[i];

        // split the line using space separator and store them to a vector of tokens
        vector<string> vsTokens;
        vsTokens = split(currentLine, string(" "));

        if (vsTokens[0].compare(m_sBuildKey) == 0)
        {
            m_fsetLattice(vsTokens);
        }
        if (vsTokens[0].compare(m_sReadKey) == 0)
        {
            //Cml reader constructor
            //xyz reader constructor
        }

        if(vsTokens[0].compare(m_sNSpeciesKey)==0){
            vector<string> specLines;
            int k=i+1;
            for(int j=k; j < k+toInt(vsTokens[1]); j++){
                specLines.push_back(fLines[j]);
            }
            m_fsetSpecies(specLines);
        }

        if(vsTokens[0].compare(m_sNProcKey)==0){
            vector<string> procLines;
            int k=i+1;
            for(int j=k; j < k+toInt(vsTokens[1]); j++){
                procLines.push_back(fLines[j]);
            }
            m_fsetProcesses(procLines);
        }

        if (vsTokens[0].compare(m_sPressureKey) == 0)
        {
            m_fsetPressure(vsTokens[1]);
        }

        if (vsTokens[0].compare(m_sStepKey) == 0)
        {
            m_fsetSteps(vsTokens);
        }


        if (vsTokens[0].compare(m_sTemperatureKey) == 0)
        {
            m_fsetTemperature(vsTokens[1]);
        }

        if (vsTokens[0].compare(m_sTimeKey) == 0)
        {
            m_fsetTime(vsTokens[1]);
        }

        if (vsTokens[0].compare(m_sDebugKey) == 0)
        {
            m_fsetDebugMode(vsTokens[1]);
        }

    }

    initializeLattice();
}

void Reader::openInputFile(string path){
    m_inputFile.open(path, ios::in);

    if (!m_inputFile.is_open())
    {
      m_errorHandler->error_simple_msg("Cannot open file input.kmc.");
      EXIT;
    }
}

void Reader::initializeLattice(){

    switch (m_LatticeType[m_sLatticeType])
    {
    case Lattice::FCC:
    {
      if (m_bSteps)
      {
          ;
        //FCC *lattice = new FCC(m_apothesis, true, m_vSteps);
        //m_apothesis->pLattice = lattice;
      }
      else
      {
        FCC *lattice = new FCC(m_apothesis);
        m_apothesis->pLattice = lattice;
        m_lattice->setType("FCC");
        break;
      }

    }
    case Lattice::BCC:
    {
      if (m_bSteps)
      {
        BCC *lattice = new BCC(m_apothesis, true, m_vSteps);
        m_apothesis->pLattice = lattice;

      }
      else
      {
        BCC *lattice = new BCC(m_apothesis);
        m_apothesis->pLattice = lattice;
      }
      m_lattice->setType("BCC");

      break;
    }
    default:
    {
      cout<<"Unresolvable lattice type found. Exiting..."<<endl;
      EXIT;
    }
    }

    m_lattice = m_apothesis->pLattice;
    m_lattice->setType(m_sLatticeType);
    m_lattice->setX(m_vLatticeDims[0]);
    m_lattice->setY(m_vLatticeDims[1]);
    m_lattice->setInitialHeight(m_vLatticeDims[2]);


}

vector<string> Reader::inputFileLines()
{

  string line;
  vector<string> lines;
  while (getline(m_inputFile, line))
  {

    // Remove any tabs, weird spaces etc.
    line = simplified(line);

    // split the line using space separator and store them to a vector of tokens
    vector<string> vsTokens;
    vsTokens = split(line, string(" "));

    //We do not care about empty lines
    if (vsTokens.size() == 0)
      continue;

    //We do not care about comments.
    if (startsWith(vsTokens[0], m_sCommentLine))
      continue;

    //Fill the vector
    lines.push_back(line);

  }
  return lines;
}

void Reader::m_fsetLattice(vector<string> vsTokens){

    m_sLatticeType=vsTokens[1];
    std::cout << "lattice type : "<< m_sLatticeType << std::endl;

    if (isNumber(vsTokens[2]))
    {
     m_vLatticeDims.push_back(toInt(vsTokens[2]));
     std::cout << "lattice x : "<<toInt(vsTokens[2]) << std::endl;

    }
    else
    {
      m_errorHandler->error_simple_msg("The x dimension of lattice is not a number.");
      EXIT;
    }

    if (isNumber(vsTokens[3]))
    {
      m_vLatticeDims.push_back(toInt(vsTokens[3]));
      std::cout << "lattice y : "<<toInt(vsTokens[3]) << std::endl;
    }
    else
    {
      m_errorHandler->error_simple_msg("The y dimension of lattice is not a number.");
      EXIT;
    }

    if (isNumber(vsTokens[4]))
    {
       m_vLatticeDims.push_back(toInt(vsTokens[4]));
       std::cout << "lattice initial height : "<<toInt(vsTokens[4]) << std::endl;
    }
    else
    {
      m_errorHandler->error_simple_msg("The height must be a number.");
      EXIT;
    }
}

void Reader::m_fsetSteps(vector<string> vsTokens){
    m_bSteps=false;
    if (isNumber(vsTokens[1]) && isNumber(vsTokens[2]) && isNumber(vsTokens[3]))
    {
        m_vSteps={toInt(vsTokens[1]), toInt(vsTokens[2]), toInt(vsTokens[3])};
        if(toInt(vsTokens[1])==0 && toInt(vsTokens[2])==0 && toInt(vsTokens[3])==0)
            m_bSteps=false;
        else
            m_bSteps=true;

    }
}

void Reader::m_fsetSpecies(vector<string> lines){

    for(string& line:lines){
        vector<string> vsTokens;
        vsTokens = split(line, string(" "));
        if(vsTokens.size() < 2){
            m_errorHandler->error_simple_msg("Missing input fields for species "+ vsTokens[0]);
        }
        else {
            if(isNumber(vsTokens[1])){
                m_mSpecies.insert({vsTokens[0],toDouble(vsTokens[1])});
            }
            else{
                m_errorHandler->error_simple_msg("Missing mw for species "+ vsTokens[0] + "is not a number");
            }
        }

    }

}

void Reader::m_fsetProcesses(vector<string> lines){
    int procId=0;
    for(string & line: lines){
        cout << line << endl;
        m_fidentifyProcess(line,procId);
        procId++;
    }
}

void Reader::m_fidentifyProcess(string processKey, int id){
    vector<string> process = split(processKey, string(","));
    vector<string> reactants, products;

    tie(reactants,products)=m_fsplitReactionKey(process[0]);
    vector<string> species=m_fprocSpecies(reactants,products);
    vector<double> energetics= m_fprocEnergetics(process[1]);
    vector<double> stoichiometry;

    string procName;
    if(m_bisAdsorption(reactants)){
        procName="Adsorption"+to_string(id);
    }else if(m_bisDesorption(products)){
        procName="Desorption"+to_string(id);
    }else if(m_bisDiffusion(reactants,products)){
        procName="Diffusion"+to_string(id);
    }else{
        procName="Reaction"+to_string(id);
        stoichiometry=m_fprocStoichiometry(reactants,products);

    }
    m_fsetProcInfo(procName,species,energetics,stoichiometry);
}

void Reader::m_fsetProcInfo(string procName, vector<string> species, vector<double> energetics, vector<double> stoichiometry){

    m_mProcSpecies.insert({procName,species});
    m_mProcEnergetics.insert({procName,energetics});
    if(!stoichiometry.empty())
        m_mProcStoichiometry.insert({procName,stoichiometry});
}


vector<string> Reader::m_fprocSpecies(vector<string> reactants, vector<string> products){
    vector<string> species;
    for(const auto& [key,value]: m_mSpecies){
         //std::cout << key << '\n';
         for(string& reactant: reactants){
            if(contains(reactant,key)){
                if (find(species.begin(), species.end(), key) == species.end()) {
                  species.push_back(key);
                }
            }
         }
         for(string& product: products){
            if(contains(product,key)){
                if (find(species.begin(), species.end(), key) == species.end()) {
                  species.push_back(key);
                }
            }
         }
    }

    return species;
}


vector<string> Reader::mapKeys(map<string,double> inMap){

    std::vector<string> vsKeys;

    for(const auto& [key,value]: inMap){
        vsKeys.push_back(key);
    }
    return vsKeys;
}

vector<double> Reader::m_fprocStoichiometry(vector<string> reactants, vector<string> products){
    vector<double> stoichiometry;
    vector<string> speciesKeys=mapKeys(m_mSpecies);

    //std::cout << key << '\n';
    for(string& reactant : reactants){
       string longestMatch=lcMatch(reactant,speciesKeys);
       eraseSubstring(reactant, longestMatch);
       if (isNumber(simplified(reactant)))
            stoichiometry.push_back(toDouble(simplified(reactant)));

    }
    for(string& product : products){
        string longestMatch=lcMatch(product,speciesKeys);
        eraseSubstring(product, longestMatch);
        if (isNumber(simplified(product)))
            stoichiometry.push_back(toDouble(simplified(product)));
    }

    return stoichiometry;
}


void Reader::eraseSubstring(string & mainStr, const string & toErase)
{
    size_t pos = string::npos;
    // Search for the substring in string in a loop untill nothing is found
    while ((pos  = mainStr.find(toErase) )!= string::npos)
    {
        // If found then erase it from string
        mainStr.erase(pos, toErase.length());
    }
}

vector<double> Reader::m_fprocEnergetics(string energeticsKey){
    vector<double> vdEnergetics;
    vector<string> vsEnergetics = split(energeticsKey, string(" "));
    for(string& str: vsEnergetics) {
        if (isNumber(str)){
            vdEnergetics.push_back(toDouble(str));
        }
    }
    return vdEnergetics;
}

pair<vector<string>,vector<string>> Reader::m_fsplitReactionKey(string reactionKey){
    vector<string> reaction = split(reactionKey, string("->"));
    vector<string> reactants=split(reaction[0],"+");
    vector<string> products=split(reaction[1],"+");
    return make_pair(reactants,products);
}

bool Reader::m_bisAdsorption(vector<string> reactants){
    return (reactants.size()>1) ? (simplified(reactants[0])==m_ssiteKey || simplified(reactants[1])==m_ssiteKey): false ;
}

bool Reader::m_bisDesorption(vector<string> products){
    return (products.size()>1) ? (simplified(products[0])==m_ssiteKey || simplified(products[1])==m_ssiteKey): false ;
}

bool Reader::m_bisDiffusion(vector<string> reactants, vector<string> products){
    return (reactants.size()==1 && products.size()==1 && contains(reactants[0],m_ssiteKey) && contains(products[0],m_ssiteKey)) ;
}

string Reader::lcMatch(string X, vector<string> vsY)
{
    map<int,string> mMatches;

    for(string& str: vsY){
        if(contains(X,str)){
            mMatches.insert({str.size(),str});
        }
    }
    int maxKey=0;
    for(const auto& [key,value]: mMatches){
        if(key > maxKey)
            maxKey=key;
    }

    return mMatches[maxKey];

}


void Reader::m_fsetTime(string time){
    if (isNumber(time))
    {
      m_dTime=toDouble(time);
      m_parameters->setEndTime(m_dTime);
      std::cout << "Time "<<toDouble(time)<<std::endl;
    }
    else
    {
      m_errorHandler->error_simple_msg("Could not read number of KMC simulation time from input file. Is it a number?");
      EXIT;
    }
}

void Reader::m_fsetPressure(string pressure){
    if (isNumber(pressure))
    {
      m_dPressure=toDouble(pressure);
      m_parameters->setPressure(m_dPressure);
      std::cout << "Pressure : "<<toDouble(pressure) << std::endl;
    }
    else
    {
      m_errorHandler->error_simple_msg("Could not read pressure from input file. Is it a number?");
      EXIT;
    }
}

void Reader::m_fsetTemperature(string temperature){
    if (isNumber(temperature))
    {
      m_dTemperature=toDouble(temperature);
      m_parameters->setTemperature(m_dTemperature);
      std::cout << "Temperature : "<<toDouble(temperature) << std::endl;
    }
    else
    {
      m_errorHandler->error_simple_msg("Could not read temperature from input file. Is it a number?");
      EXIT;
    }
}

void Reader::m_fsetDebugMode(string debugMode){
    cout << debugMode << endl;
    if (contains(debugMode,"on", Insensitive) || contains(debugMode,"off", Insensitive))
    {
      m_sDebugMode=debugMode;
      cout << "Debug mode: "<< debugMode<< endl;
    }
    else
    {
      m_errorHandler->error_simple_msg("Could not read debug mode.");
      EXIT;
    }
}

string Reader::simplified(string str)
{
  string s;
  bool is_white = false;
  bool was_white = false;
  bool append = false;
  size_t n = str.size();
  size_t nm = n - 1;

  for (size_t i = 0; i < n; i++)
  {
    char c;
    c = str[i];
    switch (c)
    {
    case ' ':
    case '\n':
    case '\t':
    case '\v':
    case '\f':
    case '\r':
      is_white = true;
      c = ' ';
      append = (was_white || i == 0 || i == nm) ? false : true;
      break;
    default:
      is_white = false;
      append = true;
      break;
    };
    if (append)
      s.push_back(c);
    was_white = is_white;
  }
  return s;
}

bool Reader::isNumber(string str)
{

  size_t n = str.size();

  bool bPlusFound = false;
  bool bPlusFound2 = false;
  bool bMinusFound = false;
  bool bMinusFound2 = false;
  bool be = false;
  bool bE = false;
  bool bd = false;
  bool bD = false;
  bool bf = false;
  bool bF = false;
  bool bDotFound = false;

  for (size_t i = 0; i < n; i++)
  {
    char c;
    c = str[i];
    switch (c)
    {
    case '+':
      if (i == 0)
      {
        bPlusFound = true;
        break;
      }

      if (i != 0 && (bPlusFound || bMinusFound) && !bD && !bd && !be && !bE && !bf && !bF)
        return false;

      if ((bD || bd || be || bE || bf || bF) && (!bPlusFound2))
        bPlusFound2 = true;
      else if ((bD || bd || be || bE || bf || bF) && bPlusFound2)
        return false;
      else if ((bD || bd || be || bE || bf || bF) && bMinusFound2)
        return false;

      break;
    case '-':
      if (i == 0)
      {
        bMinusFound = true;
        break;
      }

      if (i != 0 && (bPlusFound || bMinusFound) && !bD && !bd && !be && !bE && !bf && !bF)
        return false;

      if ((bD || bd || be || bE || bf || bF) && (!bPlusFound2))
        bPlusFound2 = true;
      else if ((bD || bd || be || bE || bf || bF) && bPlusFound2)
        return false;
      else if ((bD || bd || be || bE || bf || bF) && bMinusFound2)
        return false;

      break;
    case 'e':
      if (!bD && !bd && !be && !bE && !bf && !bF)
        be = true;
      else
        return false;
      break;
    case 'E':
      if (!bD && !bd && !be && !bE && !bf && !bF)
        bE = true;
      else
        return false;
      break;
    case 'D':
      if (!bD && !bd && !be && !bE && !bf && !bF)
        bd = true;
      else
        return false;
      break;
    case 'd':
      if (!bD && !bd && !be && !bE && !bf && !bF)
        bD = true;
      else
        return false;
      break;
    case 'F':
      if (!bD && !bd && !be && !bE && !bf && !bF)
        bf = true;
      else
        return false;
      break;
    case 'f':
      if (!bD && !bd && !be && !bE && !bf && !bF)
        bF = true;
      else
        return false;
      break;
    case '.':
      if (!bDotFound && !bD && !bd && !be && !bE && !bf && !bF)
      {
        bDotFound = true;
      }
      else if (!bDotFound && (bD || bd || be || bE || bf || bF))
        return false;
      else if (bDotFound)
        return false;
      break;
    default:
      if (!isdigit(c))
        return false;
    };
  }
  return true;
}

double Reader::toDouble(string str)
{
  if (isNumber(str))
    return stod(str);
  else
    return 0;
}

vector<string> Reader::split(string str, string delim)
{
  vector<string> r;
  size_t prev = 0, pos = 0;
  do
  {
    pos = str.find(delim, prev);
    if (pos == string::npos)
      pos = str.length();
    string token = str.substr(prev, pos - prev);
    if (!token.empty())
      r.push_back(token);
    prev = pos + delim.length();
  } while (pos < str.length() && prev < str.length());

  return r;
}

int Reader::toInt(string str)
{
  if (isNumber(str))
    return stof(str);
  else
    return 0;
}

bool Reader::contains(string str1, string str2, CASE cas)
{
  switch (cas)
  {
  case Insensitive:
  {
    char c;
    char csmall;
    string smstr, sstr;
    size_t n = str1.size();

    for (size_t i = 0; i < n; i++)
    {
      c = str1[i];
      csmall = (char)tolower(c);
      smstr.push_back(csmall);
    }

    for (size_t i = 0; i < str2.size(); i++)
    {
      c = str2[i];
      csmall = (char)tolower(c);
      sstr.push_back(csmall);
    }

    if (smstr.find(sstr) != string::npos)
      return true;
    return false;
  }
  default:
    if (str1.find(str2) != string::npos)
      return true;
    return false;
  }
}


double Reader::getTemperature(){
    return m_dTemperature;
}

double Reader::getPressure(){
    return m_dPressure;
}

double Reader::getTime(){
    return m_dTime;
}

string Reader::getDebugMode(){
    return m_sDebugMode;
}

map<string,double> Reader::getSpecies(){
    return m_mSpecies;
}


map<string,vector<string>> Reader::getProcSpecies(){
    return m_mProcSpecies;
}

map<string,vector<double>> Reader::getProcEnergetics(){
    return m_mProcEnergetics;
}


map<string,vector<double>> Reader::getProcStoichiometry(){
    return m_mProcStoichiometry;
}
