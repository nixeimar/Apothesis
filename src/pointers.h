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

#ifndef POINTERS_H
#define POINTERS_H

#include "apothesis.h"

/** The pointers class contains pointers to the classes defined in kmc.cpp.
 * Itcreates a common space of all the class instances that are defined in the kmc class.
 * As long as a kmc instance exists these classes are "united" in this common space. So,
 * if someone wants to use the lattice instance from the io instance all that has to do is
 * to call the insance as this defined in the kmc.cpp.
 * For example in kcc.cpp an instance of the class lattice is pLattice and an object of the
 * io class is pIO. If someones wants to use the pIO from pLattice all (s)he has to do is to
 * call pIO-> and the relevant public functionalities of the class. This is very handy to
 * bind together different instances. (Taken from LAMMPS ;) )*/

class Pointers
{

public:
    /// Contructor
    Pointers( Apothesis* apothesis):m_apothesis(apothesis),
                        m_lattice( apothesis->pLattice),
                        m_io(apothesis->pIO),
                        m_read(apothesis->pRead),
                        m_errorHandler(apothesis->pErrorHandler),
                        m_parameters(apothesis->pParameters),
                        m_randomGen( apothesis->pRandomGen )
    {}

protected:
    /// Pointers to the classes of apothesis.cpp
    Apothesis*& m_apothesis;

    /// Pointers to the classes of apothesis.cpp
    Lattice*& m_lattice;

    /// Pointers to the classes of apothesis.cpp
    IO*& m_io;

    /// Pointers to the classes of apothesis.cpp
    Read*& m_read;

    /// Pointers to the classes of apothesis.cpp
    Utils::ErrorHandler*& m_errorHandler;

    /// Pointers to the classes of kmc.cpp
    Utils::Parameters*& m_parameters;

    /// Random generator
    RandomGen::RandomGenerator *&m_randomGen;

};

#endif // POINTERS_H
