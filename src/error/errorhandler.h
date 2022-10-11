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

#ifndef ERRORHANDLER_H
#define ERRORHANDLER_H

#include "pointers.h"
#include "apothesis.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>


namespace Utils {

/** A class for handling the development and user errors.
  * Different ways to handle an error depending on its nature.
    Many of the functions are not implememented yet.*/

class ErrorHandler: public Pointers
  {
  public:
    /// Constructor
    ErrorHandler( Apothesis *apothesis );

    /// Destructor
    virtual ~ErrorHandler();

    /// Different ways to handle an error depending on its nature.

    /// Error message simple text
    void error_simple_msg( string msg );

    /// Error message
    void errorMsg( const string &file, const string &filetype, int line, const string &msg );

    /// Error message
    void errorMsg( const string &file, const string &filetype, int line, int column, const string &msg );

    /// Error message
    void errorMsg( const string &file, const string &filetype, const string &msg );

    /// Warning message
    void warningMsg( const string &file, const string &filetype, const string &msg );

    ///  Warning message
    void warningMsg( const string &movetype, const string &msg );

    ///  Warning message simple text
    void warningSimple_msg( const string &msg );
  };
}

#endif // ERRORHANDLER_H
