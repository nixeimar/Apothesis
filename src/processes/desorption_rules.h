#ifndef DESORPTION_RULES_H
#define DESORPTION_RULES_H

#include "desorption.h"
#include "site.h"

namespace MicroProcesses
{

/// If the keyword 'all' is used then the rule is based on the neighbours
bool allRule(Desorption*, Site* s);

/// Returns always true - this is actually as having uncoditional acceptance
bool basicRule(Desorption*, Site* s);

/// If the keyword myltilayer is included in diffusion
bool diffusionMultilayerUp(Desorption*, Site* s);

/// If the keyword myltilayer is included in diffusion
bool diffusionMultilayerDown(Desorption*, Site* s);

/// For desorbing different species the site must be occupied
bool difSpeciesRule(Desorption*, Site* s);

}

#endif // DESORPTION_RULES_H
