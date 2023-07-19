#ifndef DIFFUSION_RULES_H
#define DIFFUSION_RULES_H

#include "diffusion.h"

namespace MicroProcesses
{

/**  This is the basic rule: For any atom X which does not belong to the growing film
 *   check if there is a vacant site that can be diffused to.
**/
bool diffusionBasicRule(Diffusion*, Site* s);

/**  This is the rule when the user has used the "all" keyword in the input file.
 *   It is applied only to the atoms that belong to the growing film. (PVD only)
**/
bool diffusionAllRule(Diffusion*, Site* s);

}

#endif // DIFFUSION_RULES_H
