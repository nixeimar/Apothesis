#ifndef DESORPTION_PERFORM_H
#define DESORPTION_PERFORM_H

#include "desorption.h"
#include "site.h"

namespace MicroProcesses
{

/// The process is PVD
void singleSpeciesSimpleDesorption(Desorption*, Site*);

/// The process is CVD or ALD
void multiSpeciesSimpleDesorption(Desorption*, Site*);

}

#endif // DESORPTION_PERFORM_H
