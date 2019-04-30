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

#ifndef ABSTRACT_PROCESS_H
#define ABSTRACT_PROCESS_H

#include <string>

namespace MicroProcesses { class Process; }

/** The abstract class which is used for the process factory **/
class AbstractProcess
  {
  public:
    /// Constructor
    AbstractProcess(){};

    /// Destructor
    virtual ~AbstractProcess(){}

    /// Pure virtual method for creating a process.
    virtual MicroProcesses::Process* create()= 0;
  };

#endif // ABSTRACT_PROCESS_H
