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


#include "cvd.h"

CVD::CVD() { ; }
CVD::~CVD()
{

}

void CVD::init()
{
    if ( io->getParameters()->getGrowthSpecies().empty()  ) {
        apothesis->pErrorHandler->error_simple_msg("At least one growing species must be defined for a CVD process."
                                                   "To define a growth species insert \"growth: X\" in the input file"
                                                   " where X is the species of the growing surface.");
        EXIT
    }

    if ( !io->getParameters()->getEtchedSpecies().size()  ) {
        apothesis->pErrorHandler->error_simple_msg("No etching species must be defined in a CVD process");
        EXIT
    }
}
