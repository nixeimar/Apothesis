#ifndef DIFFUSION_PERFORM_H
#define DIFFUSION_PERFORM_H

#include "diffusion.h"

namespace MicroProcesses
{

/** This is the simplest of diffusion.
 *  It takes particle X and moves it in a vacant site from its first neighbors
**/
void simpleDiffusion( Diffusion*, Site*);

/// Diffuse one layer down
void pefrormMultilayerUp( Diffusion*, Site*);

/// Diffuse one layer up
void pefrormMultilayerDown( Diffusion*, Site*);

/// The process is PVD as in Lam and Vlachos (2000)
void performPVD( Diffusion*, Site*);

/// ToDo: Add dimer diffusion

}

#endif // DIFFUSION_PERFORM_H
