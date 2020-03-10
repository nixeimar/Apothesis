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

#ifndef IO_H
#define IO_H

#include <fstream>


#include <string>
#include <vector>
#include <map>
#include <algorithm>


#include "pointers.h"
#include "apothesis.h"
#include "lattice.h"

#include "lattice.h"
#include "errorhandler.h"
#include "parameters.h"

#if defined( _WIN32) || defined( _WIN64)
#include <direct.h>
#include <io.h>
#define GetCurrentDir _getcwd
#define BAC '\\'
#elif ( __linux__)
#include <unistd.h>
#include <dirent.h>
#define GetCurrentDir getcwd
#define BAC '/'
#endif

using namespace std;
using namespace Utils;

/** Tha class for hnadlin input/output operations */

class IO: public Pointers
  {
  public:
    enum CASE{ Sensitive, Insensitive };

    IO();
    IO( Apothesis *apothesis);

    virtual ~IO();

    /// Initialization of the files: Assigns paths and file names
    /// If they are not given by the user, the default values are used
    /// Dafault: The path of the executable
    /// Input file name: input.kmc
    /// Output file name: output.log
    void init( int argc, char* argv[] );

    /// Returns the path of the input file.
    string getInputPath() const;

    /// Returns the path of the output file.
    const string& getOutputPath() const;

    /// Returns the name of the input file if it is user defined.
    const string& getInputFilename() const;

    /// Returns the name of the output file if it is user defined.
    const string& getOutputFilename() const;

    /// Opens the input file.
    void openInputFile();

    /// Opens the output file with the name name.
    bool openOutputFile( string name );

    /// Closes the output file.
    void closeOutputFile();

    /// Reads the input file " .kmc".
    void readInputFile();

    /// Converts a string to double.
    double toDouble( string );

    /// Converts a string to int.
    int toInt( string );

    /// Check if a string contains another string. TODO: This should be transferred to a generic string class).
    bool contains( string, string, CASE cas = Insensitive );

    /// Splits a string to a vector of strings. TODO: This should be transferred to a generic string class).
    vector<string> split( string , string );

    /// Given a string returns a string with all the delimeters replaced. TODO: This should be transferred to a generic string class).
    string simplified( string );

    /// Returns if the given string in number. TODO: This should be transferred to a generic string class).
    bool isNumber( string );

    /// Checks if a file exists.
    bool exists(const string& s);

    /// Converts a double to string. TODO: This should be transferred to a generic string class).
    std::string convertToString( double x);

    /// Converts an integer to string. TODO: This should be tranwriteIterationInfosferred to a generic string class).
    std::string convertToString( int x);

    /// Return the working directory
    string GetCurrentWorkingDir();

    /// Specific input for the lattice.
    Lattice::Type getLatticeType();

    ///  Check of a string starts with a certain string. TODO: This should be transferred to a generic string class).
    inline bool startsWith( string str, string substr )
      {
      if ( str.find( substr ) != string::npos )
        return true;
      else
        return false;
      }

    /// Write in the output file.
    /// For this an appropriate format must be defined.
    /// Now it is only writes a simple string.
    void writeLogOutput( string );

    /// Write lattice info.
    void writeLatticeInfo();

    /// Write the height of each site
    void writeLatticeHeights();

    /// Export the lattice in xyz format. Not implemented yet
    void exportLatticeXYZ();

    /// Export the lattice in cml format. Not implemented yet.
    void exportLatticeCML();

  protected:
    /// The type of lattice
    string m_sLatticeType;

    /// Supported lattice types
    map< string, Lattice::Type> m_mLatticeType;

    /// The processes to be constructed with their parameters (per process)
    map< string, list< string > >  m_mProc;

    /// The input file
    ifstream m_InputFile;

    /// The output file
    ofstream m_OutFile;

    /// Keywords:
    /// Process keyword
    string m_sProcess;

    /// Lattice keyword
    string m_sLattice;

    /// Temperature keyword
    string m_sTemperature;

    /// Pressure keyword
    string m_sPressure;

    ///  Num of iterations keyword.
    string m_sIterations;

    /// Comment
    string m_sCommentLine;

  private:
    /// The number of flags given by the user
    int m_iArgc;

    /// The flags given by the user
    char** m_vcArgv;

  };

#endif // IO_H
