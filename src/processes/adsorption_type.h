#ifndef ADSORPTION_TYPE_H
#define ADSORPTION_TYPE_H

#include "adsorption.h"

/// The simple type for the adsorption process rate i.e.
/// simple s0*f*P/(2*pi*MW*Ctot*kb*T) -> Sticking coefficient [-], f [-], C_tot [sites/m2], MW [kg/mol]
void simpleType(Adsorption*);

/// The arrhenius type for the adsorption process rate i.e.
/// arrhenius v0 A exp(-nE/kT), A = exp((E-Em)/kT) -> frequency v0 [-],  E (Joules), Em [Joules]
void arrheniusType(Adsorption* );

/// Constant value for the adsorption process rate i.e.
/// constant 1.0 [ML/s]
void constantType(Adsorption*);

#endif // ADSORPTION_TYPE_H
