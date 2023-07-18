#ifndef DIFFUSION_PERFORM_H
#define DIFFUSION_PERFORM_H

#include "diffusion.h"
#include "site.h"

/** This is the simplest of diffusion.
 *  It takes particle X and moves it in a vacant site from its first neighbors
**/
void simpleDiffusion(Diffusion*, Site*);

/// The process is PVD as in Lam and Vlachos (2000)
void performPVD(Diffusion*, Site*);

/// ToDo: Add dimer diffusion

#endif // DIFFUSION_PERFORM_H
