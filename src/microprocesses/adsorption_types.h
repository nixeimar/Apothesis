#ifndef ADSORPTION_TYPES_H
#define ADSORPTION_TYPES_H

#include "adsorption.h"

namespace MicroProcesses
{

/// The simple type for the adsorption process rate i.e.
/// simple s0*f*P/(2*pi*MW*Ctot*kb*T) -> Sticking coefficient [-], f [-], C_tot [sites/m2], MW [kg/mol]
double simpleType( Adsorption* );

/// The arrhenius type for the adsorption process rate i.e.
/// arrhenius v0 A exp(-nE/kT), A = exp((E-Em)/kT) -> frequency v0 [-],  E (Joules), Em [Joules]
double arrheniusType( Adsorption* );

/// Constant value for the adsorption process rate i.e.
/// constant 1.0 [ML/s]
double constantType( Adsorption* );

}

#endif // ADSORPTION_TYPES_H
