#include "adsorption_new.h"

Adsorption_new::Adsorption_new(){}

Adsorption_new::~Adsorption_new(){}

void Adsorption_new::perform()
{
    m_pLattice->adsorp( m_Site->getID(), m_Species );
}

