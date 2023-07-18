#ifndef ADSORPTION_RULES_H
#define ADSORPTION_RULES_H

#include "site.h"
#include "adsorption.h"

/// The uncoditional rule. The process is accepted without checked.
inline bool uncoRule(Adsorption*, Site*) { return true; }

/// The basic rule for accepting this process.
/// Check if the site is empty (i.e. the label is the same as the lattice species)
/// then returns true (the processes can be performed).
bool basicRule( Adsorption*, Site*);

/// For adsorbing different species in a single site must not be occupied (and TODO: the height must be the same)
bool multiSpeciesSimpleRule(Adsorption*, Site*);

/// For adsorbing different species the sites must not be occupied (and TODO: the height must be the same)
bool multiSpeciesRule(Adsorption*, Site*);


#endif // ADSORPTION_RULES_H
