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

#ifndef REGISTER_LATTICE_H
#define REGISTER_LATTICE_H

#include "abstract_lattice.h"
#include "lattice.h"

/** A template class used by the factory method for creating the processes.*/

template<class T>
class RegisterLattice: public AbstractLattice// Pointers class contains ptrs to master copy of
{
public:
    /// Contructor
//    Register<T>( const std::string& name);
    Register<T>(const std::string& name):AbstractLattice( name ){}

    /// Destructor
    virtual ~Register<T>(){}

    /// Creator of the lattice.
    virtual Lattice* create(){
        return new T;
    }

};


#endif // REGISTER_LATTICE_H
