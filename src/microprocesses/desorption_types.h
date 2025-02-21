#ifndef DESORPTION_TYPES_H
#define DESORPTION_TYPES_H

#include "desorption.h"

namespace MicroProcesses
{

/// Arrhenius type rate
double arrheniusType( Desorption* );

/// Constant value for the adsorption process rate i.e. constant 1.0 [ML/s]
double constantType( Desorption* );


}
#endif // DESORPTION_TYPES_H
