#ifndef DIFFUSION_TYPES_H
#define DIFFUSION_TYPES_H

#include "diffusion.h"

namespace MicroProcesses
{

/// Constant value for the diffusion process rate i.e.
double constantType( Diffusion* );

/// Arrhenius type
double arrheniusType( Diffusion* );

}

#endif // DIFFUSION_TYPES_H
