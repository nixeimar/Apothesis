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


#include "ald.h"

ALD::ALD() { ; }
ALD::~ALD() { ; }

void ALD::init()
{

    cout << "Init CVD process ...";

    Utils::Parameters* basic = io->getParameters();

    for (auto const& [key, p ] : io->getCycles())
        io->mergeParameters( p, basic );

    m_mCycles = io->getCycles();

    apothesis->setIO( io );
}

bool ALD::createWorkingDir(const std::string& dirName) {
    try {
        // Create the directory if it doesn't exist
        if (!fs::exists(dirName)) {
            if (!fs::create_directory(dirName)) {
                std::cerr << "Failed to create directory: " << dirName << std::endl;
                return false;
            }
        }

        // Set the current working directory
        fs::current_path(dirName);
        std::cout << "Current working directory changed to: " << fs::current_path() << std::endl;

        return true;
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
        return false;
    }
}

void ALD::perform()
{
    int iCycle = 0;
    for (auto const& [key, p ] : m_mCycles ) {

        std::string originalDir = fs::current_path().string();
        createWorkingDir(  "Test_" + std::to_string(iCycle) );

        Lattice* l = apothesis->getLattice();
        apothesis->update( p, l );
        apothesis->exec();

        //Update the parameters
        p->setStartTime( apothesis->getEndTime() );
        p->setEndTime( apothesis->getStartTime() + p->getEndTime() );

        // Restore the original working directory
        fs::current_path(originalDir);

        iCycle++;
    }

}
