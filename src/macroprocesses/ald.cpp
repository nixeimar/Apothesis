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

    cout << "Init ALD process ...";

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
    int iNumCycles = io->getParameters()->getNumCycles();
    double startTime = 0.0;

    std::string originalDir = fs::current_path().string();

    apothesis->initRandom( io->getParameters() );

    for (int i = 0; i < iNumCycles; i++ ) {
        for (auto const& [key, p ] : m_mCycles ) {

            createWorkingDir(  "Cycle_" + std::to_string( startTime ) );

            //Open the output file
            Lattice* l = apothesis->getLattice();

            apothesis->update( p, l,  startTime );
            apothesis->exec();

            startTime = 0.0;
            startTime = apothesis->getEndTime();

            iCycle++;

            // Restore the original working directory
            fs::current_path(originalDir);
        }
    }


    std::cout << "Apothesis finished succesfully ..." << std::endl;

}
