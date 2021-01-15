#include "desorption_new.h"

Desorption_new::Desorption_new(){}

Desorption_new::~Desorption_new(){}

void Desorption_new::perform()
{
    m_pLattice->desorp( m_Site->getID(), m_Species );
}
