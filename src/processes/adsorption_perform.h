#ifndef ADSORPTION_PERFORM_H
#define ADSORPTION_PERFORM_H

#include "adsorption.h"
#include "site.h"

namespace MicroProcesses
{

/// The process is PVD
void signleSpeciesSimpleAdsorption(Adsorption*, Site* );

/// The process is PVD for multiple sites
void signleSpeciesAdsorption(Adsorption*, Site*);

/// The process is CVD or ALD
void multiSpeciesSimpleAdsorption(Adsorption*, Site*);

/// The process is CVD or ALD for multiple sites
void multiSpeciesAdsorption(Adsorption*, Site*);

}

#endif // ADSORPTION_PERFORM_H
