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

#ifndef FACTORY_LATTICE_H
#define FACTORY_LATTICE_H

#include <string>
#include <map>

#include "register_lattice.h"

class Lattice;
class AbstractLattice;
// class Apothesis;

using namespace std;

class FactoryLattice
{
    public:
    /// Constructor
    FactoryLattice();

    /// Destructor
    virtual ~FactoryLattice();

    /// Register the new lattice (Store it in the map table since a map table cannot be initialized)
    static void registerThis( const string& , AbstractLattice* );

    /// Create the lattice
    static Lattice* createLattice( const string& );

private:
    /// The factory map
    static map< string, AbstractLattice* >& getTable();

};

#define REGISTER_LATTICE( __NAME__ ) \
    private: static const RegisterLattice< __NAME__> creator;

#define REGISTER_LATTICE_IMPL( __NAME__) \
    const RegisterLattice< __NAME__ > __NAME__::creator(#__NAME__);

#endif // FACTORY_LATTICE_H