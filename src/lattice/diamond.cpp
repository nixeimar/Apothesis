diamond::diamond(Apothesis* apothesis):Lattice(apothesis){}

void diamond::setInitialHeight( int  height ) { m_iHeight = height; }

void diamond::build(){
 // The sites of the lattice.
    m_vSites.resize(getSize());
    for (int i = 0; i < m_vSites.size(); i++)
        m_vSites[i] = new Site();

    // This is OK
    for (int i = 0; i < m_iSizeX; i++)
    {
        for (int j = i * m_iSizeY; j < (m_iSizeY + i * m_iSizeY); j++)
        {
            m_vSites[j]->setID(j);
            m_vSites[j]->setHeight(m_iHeight);
        }
    }

    mf_neigh();

    // Here we set the label of the species
    for (int i = 0; i < m_iSizeY; i++)
    {
        for (int j = 0; j < m_iSizeX; j++)
            m_vSites[i * m_iSizeX + j]->setLabel(m_sLabel);
    }
}

void diamond:mf_neigh() {
	
   // The sites of the lattice.
    m_vSites.resize(getSize());

        for (int site = 0; site < getSize(); site++) {
	    int x = site / L;
            int y = site % L;
         
            m_vSites[ site ]->setNeigh( m_vSites[ ((x + 1) % L) * L + y ]); //Right
            m_vSites[ site ]->setNeigh( m_vSites[ x * L + ((y + 1) % L) ]); //Down
            m_vSites[ site ]->setNeigh( m_vSites[ ((x - 1 + L) % L) * L + y ]); //Left
            m_vSites[ site ]->setNeigh( m_vSites[ x * L + ((y - 1 + L) % L) ]) //Up
    }
}
