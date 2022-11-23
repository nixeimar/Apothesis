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

#include "reaction.h"

reaction::reaction(){}
reaction::~reaction(){}


void reaction::perform(Site *)
{
    ;

}

bool reaction::rules(Site *s)
{
    // Check if the site
    if ( m_mStoichiometry.find( s->getLabel() ) == m_mStoichiometry.end() )
        return false;

    for (auto const& [key, val] : m_mStoichiometry) {
        int iCount = 0;

        if (s->getLabel().compare( key ) == 0)
            iCount++;

        for (Site* site:s->getNeighs()){
            if ( site->getLabel().compare( key ) == 0)
                iCount++;
        }

        if ( iCount < val)
            return false;
    }

    return true;
}

double reaction::getProbability(){
    return 0.;
}


