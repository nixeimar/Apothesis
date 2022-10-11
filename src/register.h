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

#ifndef REGISTER_H
#define REGISTER_H

#include "abstract_process.h"
#include "process.h"

/** A template class used by the factory method for creating the processes.*/

template<class T>
class Register: public AbstractProcess// Pointers class contains ptrs to master copy of
{
public:
    /// Contructor
//    Register<T>( const std::string& name);
    Register<T>(const std::string& name):AbstractProcess( name ){}

    /// Destructor
    virtual ~Register<T>(){}

    /// Creator of the process.
    //virtual MicroProcesses::Process* create();
    virtual MicroProcesses::Process* create(){
        return new T;
    }

};

#endif // REGISTER_H
