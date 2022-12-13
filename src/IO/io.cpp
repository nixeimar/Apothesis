//============================================================================
//    Apothesis: A kinetic Monte Calro (KMC) code for deposition processes.
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

#include "io.h"

IO::IO(Apothesis* apothesis):Pointers(apothesis),
    m_sLatticeType("NONE"),
    m_sProcess("process"),
    m_sLattice("lattice"),
    m_sTemperature("temperature"),
    m_sPressure("pressure"),
    m_sTime("time"),
    m_sSteps("steps"),
    m_sRandom("random"),
    m_sSpecies("species"),
    m_sWrite("write"),
    m_sAdsorption("Adsorption"),
    m_sDesorption("Desorption"),
    m_sDiffusion("Diffusion"),
    m_sCommentLine("#")
{
    //Initialize the map for the lattice
    m_mLatticeType[ "NONE" ] = Lattice::NONE;
    m_mLatticeType[ "SimpleCubic" ] = Lattice::SimpleCubic;
    m_mLatticeType[ "FCC" ] = Lattice::FCC;
}

IO::~IO(){}

void IO::init(int argc, char* argv[])
{
    openInputFile("input.kmc");
}


string IO::getInputPath() const {;}


void IO::readInputFile()
{
    list< string > lKeywords{ m_sLattice, m_sPressure, m_sTemperature, m_sTime, m_sSteps, m_sRandom, m_sSpecies, m_sWrite,
                              m_sAdsorption, m_sDesorption, m_sDiffusion};

    string sLine;
    while ( getline( m_InputFile, sLine ) ) {

        // Remove any tabs, weird spaces etc.
        sLine = simplified( sLine );
        //cout << sLine << endl;

        //We do not care about comments.
        if ( startsWith( sLine, m_sCommentLine ) ) continue;

        // split the line using : separator and store them to a vector of tokens
        vector<string> vsTokensBasic;
        vsTokensBasic = split( sLine, string( ":" ) );

        //We do not care about empty lines
        if ( vsTokensBasic.size() == 0 ) continue;

        bool bComment = false;
        for ( int i = 0; i< vsTokensBasic.size(); i++){
            if ( !bComment && startsWith( vsTokensBasic[ i ], m_sCommentLine ) )
                bComment = true;

            // Remove the comments from the tokens so not to consider them
            if ( bComment )
                vsTokensBasic[ i ].clear();
        }

        // Remove any empty parts of the vector
        vector<string>::iterator it = remove_if( vsTokensBasic.begin(), vsTokensBasic.end(), mem_fun_ref(&string::empty) );
        vsTokensBasic.erase( it, vsTokensBasic.end() );

        // Check if a token is not a keyword
        if ( find( lKeywords.begin(), lKeywords.end(), vsTokensBasic[ 0 ] ) == lKeywords.end() ){
            string msg = "Unknown keyword ( " + vsTokensBasic[ 0 ] + " )";
            m_errorHandler->error_simple_msg( msg );
            EXIT;
        }

        if ( vsTokensBasic[ 0].compare(  m_sLattice ) == 0 ){

            vector<string> vsTokens;
            vsTokens = split( vsTokensBasic[ 1 ], string( " " ) );

            bool bComment = false;
            for ( int i = 0; i< vsTokens.size(); i++){
                if ( !bComment && startsWith( vsTokens[ i ], m_sCommentLine ) )
                    bComment = true;

                // Remove the comments from the tokens so not to consider them
                if ( bComment )
                    vsTokens[ i ].clear();
            }

            int count = 0;

            m_lattice->setType( vsTokens[ 0 ] );

            if ( isNumber( vsTokens[ 1 ] ) ){
                m_lattice->setX( toInt(  trim( vsTokens[ 1 ] ) ) );
            }
            else {
                m_errorHandler->error_simple_msg("The x dimension of lattice is not a number.");
                EXIT;
            }

            if ( isNumber( vsTokens[ 2 ] ) ){
                m_lattice->setY( toInt( trim( vsTokens[ 2 ] ) ) );
            }
            else {
                m_errorHandler->error_simple_msg("The y dimension of lattice is not a number.");
                EXIT;
            }

            if ( isNumber( vsTokens[ 3 ] ) ){
                m_lattice->setInitialHeight( toInt( trim( vsTokens[ 3 ] ) ) );
            }
            else {
                m_errorHandler->error_simple_msg("The height must be a  number.");
                EXIT;
            }
        }

        if ( vsTokensBasic[ 0].compare(  m_sSteps ) == 0 ){
            m_lattice->setSteps( true );

            vector<string> vsTokens;
            vsTokens = split( vsTokensBasic[ 1 ], string( " " ) );

            bool bComment = false;
            for ( int i = 0; i< vsTokens.size(); i++){
                if ( !bComment && startsWith( vsTokens[ i ], m_sCommentLine ) )
                    bComment = true;

                // Remove the comments from the tokens so not to consider them
                if ( bComment )
                    vsTokens[ i ].clear();
            }

            if ( isNumber( vsTokens[ 0 ] ) ){
                m_lattice->setNumSteps( toInt( vsTokens[ 0 ] ) );
            }
            else {
                m_errorHandler->error_simple_msg("The x dimension of step is not a number.");
                EXIT;
            }

            if ( isNumber( vsTokens[ 1 ] ) ){
                m_lattice->setStepHeight( toInt( vsTokens[ 1 ] ) );
            }
            else {
                m_errorHandler->error_simple_msg("The y dimension of step is not a number.");
                EXIT;
            }
        }

        if ( vsTokensBasic[ 0].compare(  m_sTemperature ) == 0 ){

            bool bComment = false;
            for ( int i = 0; i< vsTokensBasic.size(); i++){
                if ( !bComment && startsWith( vsTokensBasic[ i ], m_sCommentLine ) )
                    bComment = true;

                // Remove the comments from the tokens so not to consider them
                if ( bComment )
                    vsTokensBasic[ i ].clear();
            }

            if ( isNumber( trim(vsTokensBasic[ 1 ] ) )){
                m_parameters->setTemperature( toDouble( trim(vsTokensBasic[ 1] ) ) );
            }
            else {
                m_errorHandler->error_simple_msg("Could not read temperature from input file. Is it a number?");
                EXIT;
            }
        }

        if ( vsTokensBasic[ 0].compare(  m_sPressure ) == 0 ){

            bool bComment = false;
            for ( int i = 0; i< vsTokensBasic.size(); i++){
                if ( !bComment && startsWith( vsTokensBasic[ i ], m_sCommentLine ) )
                    bComment = true;

                // Remove the comments from the tokens so not to consider them
                if ( bComment )
                    vsTokensBasic[ i ].clear();
            }

            if ( isNumber( trim(vsTokensBasic[ 1 ] ) ) ){
                // Is this an error?
                m_parameters->setPressure( toDouble( trim(vsTokensBasic[ 1] ) ));
            }
            else {
                m_errorHandler->error_simple_msg("Could not read pressure from input file. Is it a number?");
                EXIT;
            }
        }

        if ( vsTokensBasic[ 0].compare( m_sTime ) == 0 ){

            bool bComment = false;
            for ( int i = 0; i< vsTokensBasic.size(); i++){
                if ( !bComment && startsWith( vsTokensBasic[ i ], m_sCommentLine ) )
                    bComment = true;

                // Remove the comments from the tokens so not to consider them
                if ( bComment )
                    vsTokensBasic[ i ].clear();
            }

            if ( isNumber(  trim( vsTokensBasic[ 1 ] ) ) ){
                m_parameters->setEndTime( toDouble( trim( vsTokensBasic[ 1] ) ) );
            }
            else {
                m_errorHandler->error_simple_msg("Could not read number of KMC simulation time from input file. Is it a number?");
                EXIT;
            }
        }

        if ( vsTokensBasic[ 0].compare( m_sRandom ) == 0){
            m_parameters->setRandGenInit( toDouble( trim(vsTokensBasic[ 1] ) ) );
        }

        if ( vsTokensBasic[ 0].compare( m_sWrite ) == 0){

            vector<string> vsTokens;
            vsTokens = split( vsTokensBasic[ 1 ], string( " " ) );

            bool bComment = false;
            for ( int i = 0; i< vsTokens.size(); i++){
                if ( !bComment && startsWith( vsTokens[ i ], m_sCommentLine ) )
                    bComment = true;

                // Remove the comments from the tokens so not to consider them
                if ( bComment )
                    vsTokens[ i ].clear();
            }

            if ( vsTokens[ 0].compare( "log") == 0 ) {
                if ( isNumber( trim(vsTokens[ 1 ] ) ) ){
                    m_parameters->setWriteLogTimeStep( toDouble( trim(vsTokens[ 1 ] )) );
                }
                else {
                    m_errorHandler->error_simple_msg("Could not read number for writing to log. Is it a number?");
                    EXIT;
                }
            }
            else if ( vsTokens[ 0 ].compare( "lattice") == 0 ) {
                if ( isNumber( trim(vsTokens[ 1 ] ) ) ){
                    m_parameters->setWriteLatticeTimeStep( toDouble( trim(vsTokens[ 1 ] ) ) );
                }
                else {
                    m_errorHandler->error_simple_msg("Could not read number for writing the lattice. Is it a number?");
                    EXIT;
                }
            }
            else {
                m_errorHandler->error_simple_msg("Not correct keyword for writer. Available selections are: \"log\" and \"lattice\"");
                EXIT;
            }
        }


        if ( vsTokensBasic[ 0].compare( "Adsorption" ) == 0  ||
             vsTokensBasic[ 0].compare( "Desorption" ) == 0  ||
             vsTokensBasic[ 0].compare( "Diffusion" ) == 0 ) {
            //Set the processes to be created along with their parameters
            vector< string > tempVec;
            // We want the parameters
            // First is the keyword process, then the name of the process as this is defined in the REGISTER_PROCESS( e.g. Adsorption)
            // and then after 2 the parameters follow

            vector<string> vsTokens;
            vsTokens = split( vsTokensBasic[ 1 ], string( " " ) );

            bool bComment = false;
            for ( int i = 0; i< vsTokens.size(); i++){
                if ( !bComment && startsWith( vsTokens[ i ], m_sCommentLine ) )
                    bComment = true;

                // Remove the comments from the tokens so not to consider them
                if ( bComment )
                    vsTokens[ i ].clear();
            }

            for ( int i = 0; i < vsTokens.size(); i++ ){
                tempVec.push_back(  vsTokens[ i ] );
            }
            m_parameters->setProcess( vsTokensBasic[ 0 ], tempVec );
        }

/*        if ( vsTokens[ 0].compare( m_sProcess ) == 0){
            //Set the processes to be created along with their parameters
            vector< string > tempVec;
            // We want the parameters
            // First is the keyword process, then the name of the process as this is defined in the REGISTER_PROCESS( e.g. Adsorption)
            // and then after 2 the parameters follow
            for ( int i = 2; i < vsTokens.size(); i++ ){
                tempVec.push_back(  vsTokens[ i ] );
            }
            m_parameters->setProcess( vsTokens[ 1 ], tempVec );
        }*/
    }//Reading the lines
}

void IO::openInputFile( string file )
{
    m_InputFile.open(file, ios::in );

    if ( !m_InputFile.is_open() ) {
        m_errorHandler->error_simple_msg( "Cannot open file input.kmc" ) ;
        EXIT;
    }
}

void IO::openRoughnessFile( string file )
{
    m_RoughnessFile.open(file, ios::out );

    if ( !m_RoughnessFile.is_open() ) {
        m_errorHandler->error_simple_msg( "Cannot open file " + file ) ;
        EXIT;
    }
}

string IO::simplified( string str )
{
    string s;
    char c;
    bool is_white = false;
    bool was_white = false;
    bool append = false;
    size_t i = 0;
    size_t n = str.size();
    size_t nm = n - 1;
    string::iterator it = str.begin();

    for (size_t i = 0; i < n; i++) {
        c = str[i];
        switch (c) {
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
        default: is_white = false;
            append = true;
            break;
        };
        if (append) s.push_back(c);
        was_white = is_white;
    }
    return s;
}

bool IO::isNumber( string str){
    char c;
    size_t i = 0;
    size_t n = str.size();
    string::iterator it = str.begin();

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

    for (size_t i = 0; i < n; i++) {
        c = str[i];
        switch (c) {
        case '+':
            if ( i == 0 ){
                bPlusFound = true;
                break;
            }

            if ( i != 0 && (bPlusFound || bMinusFound) && !bD && !bd && !be && !bE && !bf && !bF)
                return false;

            if ( (bD || bd || be || bE || bf || bF) && (!bPlusFound2) )
                bPlusFound2 = true;
            else if ( (bD || bd || be || bE || bf || bF) && bPlusFound2 )
                return false;
            else if ( (bD || bd || be || bE || bf || bF) && bMinusFound2 )
                return false;

            break;
        case '-':
            if ( i == 0 ){
                bMinusFound = true;
                break;
            }

            if ( i != 0 && (bPlusFound || bMinusFound) && !bD && !bd && !be && !bE && !bf && !bF)
                return false;

            if ( (bD || bd || be || bE || bf || bF) && (!bPlusFound2) )
                bPlusFound2 = true;
            else if ( (bD || bd || be || bE || bf || bF) && bPlusFound2 )
                return false;
            else if ( (bD || bd || be || bE || bf || bF) && bMinusFound2 )
                return false;

            break;
        case 'e':
            if ( !bD && !bd && !be && !bE && !bf && !bF)
                be = true;
            else
                return false;
            break;
        case 'E':
            if ( !bD && !bd && !be && !bE && !bf && !bF)
                bE = true;
            else
                return false;
            break;
        case 'D':
            if ( !bD && !bd && !be && !bE && !bf && !bF)
                bd = true;
            else
                return false;
            break;
        case 'd':
            if ( !bD && !bd && !be && !bE && !bf && !bF)
                bD = true;
            else
                return false;
            break;
        case 'F':
            if ( !bD && !bd && !be && !bE && !bf && !bF)
                bf = true;
            else
                return false;
            break;
        case 'f':
            if ( !bD && !bd && !be && !bE && !bf && !bF)
                bF = true;
            else
                return false;
            break;
        case '.':
            if ( !bDotFound && !bD && !bd && !be && !bE && !bf && !bF){
                bDotFound = true;
            }
            else if ( !bDotFound && (bD || bd || be || bE || bf || bF))
                return false;
            else if ( bDotFound )
                return false;
            break;
        default:
            if ( !isdigit( c))
                return false;
        };
    }
    return true;
}

/// Splits a string to a vector of strings
vector<string> IO::split( string str, string delim ){
    vector<string> r;
    size_t prev = 0, pos = 0;
    do {
        pos = str.find(delim, prev);
        if (pos == string::npos) pos = str.length();
        string token = str.substr(prev, pos-prev);
        if (!token.empty()) r.push_back(token);
        prev = pos + delim.length();
    }
    while (pos < str.length() && prev < str.length());

    return r;
}


double IO::toDouble( string str )
{
    if ( isNumber( str ) )
        return stod( str);
    else
        return 0;
}

int IO::toInt( string str )
{
    if ( isNumber( str ) )
        return stof(str);
    else
        return 0;
}

bool IO::contains(string str1, string str2, CASE cas )
{
    switch (cas) {
    case Insensitive:{
        char c;
        char csmall;
        string smstr, sstr;
        size_t n = str1.size();

        for (size_t i = 0; i < n; i++) {
            c = str1[i];
            csmall = (char)tolower( c );
            smstr.push_back( csmall);
        }

        for (size_t i = 0; i < str2.size(); i++) {
            c = str2[i];
            csmall = (char)tolower( c );
            sstr.push_back( csmall);
        }

        if ( smstr.find(  sstr ) != string::npos )
            return true;
        return false;
    }
    default:
        if ( str1.find( str2 ) != string::npos )
            return true;
        return false;
    }
}

/// Opens the output file
bool IO::openOutputFile( string name )
{
    m_OutFile.open( name + ".log" , ios::out );
    if ( m_OutFile.is_open() )
        return true;

    m_errorHandler->error_simple_msg( "Cannot open file log for writting." ) ;
    exit(-1);
    //    return false;
}


void IO::writeInOutput( string toWrite )
{
    m_OutFile << toWrite << endl;
    m_OutFile.flush();
}

void IO::closeOutputFile()
{
    if ( m_OutFile.is_open( ) )
        m_OutFile.close();
}

void IO::writeLogOutput( string str )
{
    m_OutFile << str << endl;
}

void IO::writeRoughness( double t, double r)
{
    m_RoughnessFile <<  t  << "\t" << r << endl;
}

void IO::writeLatticeInfo()
{
    m_OutFile << "Lattice size: " << m_lattice->getX() << "x" << m_lattice->getY();

    if ( m_lattice->getType() == Lattice::FCC )
        m_OutFile << "Lattice type: " << "FCC";
}

void IO::writeLatticeHeights( double time  )
{
    std::string name="Lattice_" + std::to_string( time ) + ".data";
    std::ofstream file(name);

    file << "Time (s): " << time << endl;

    for (int i = 0; i < m_lattice->getY(); i++){
        for (int j = 0; j < m_lattice->getX(); j++)
            file << m_lattice->getSite( i*m_lattice->getX() + j )->getHeight() << " " ;

        file << endl;
    }
}

string IO::GetCurrentWorkingDir()
{
    char buff[FILENAME_MAX];
    string current_working_dir( buff );
    return current_working_dir;
}

Lattice::Type IO::getLatticeType()
{
    map<string, Lattice::Type>::iterator itMap = m_mLatticeType.find( m_sLatticeType);
    if ( itMap == m_mLatticeType.end() )
        return m_mLatticeType[ "NONE"];

    return m_mLatticeType[ m_sLatticeType ];
}

bool IO::outputOpen()
{
    return m_OutFile.is_open();
}

void IO::closeRoughnessFile()
{
    if ( m_RoughnessFile.is_open( ) )
        m_RoughnessFile.close();
}


