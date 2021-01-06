#include "adsorption_new.h"

adsorption_new::adsorption_new(){}

adsorption_new::~adsorption_new(){}

void adsorption_new::perform()
{
    m_pLattice->adsorp( m_Site->getID(), m_Species );
}

