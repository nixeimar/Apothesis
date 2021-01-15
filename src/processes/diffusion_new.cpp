#include "diffusion_new.h"

Diffusion_new::Diffusion_new(){}
Diffusion_new::~Diffusion_new(){}

void Diffusion_new::perform()
{
   m_pLattice->desorp( m_originSite->getID(), m_Species );
   m_pLattice->adsorp( m_targetSite->getID(), m_Species );
}
