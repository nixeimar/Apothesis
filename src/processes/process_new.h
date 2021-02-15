#ifndef PROCESS_NEW_H
#define PROCESS_NEW_H

#include <iostream>
#include <string>
#include <map>
#include <any>
#include "lattice.h"
#include "site.h"
#include "species_new.h"

using namespace std;

class Process_new
{
public:
    Process_new();
    ~Process_new();

    inline void setName( string procName ){ m_sProcName = procName; }
    inline string getName(){ return  m_sProcName; }

    inline void setID( int id ){ m_iID = id; }
    inline int getID(){ return m_iID; }

    inline void setType( int id ){ m_iID = id; }
    inline string getType( int id ){ return m_sType; }

    inline void setTargetSite( int id );

    ///Get probability
    virtual double getProbability() = 0;

    /// Perform this process in the site
    virtual void perform( int siteID ) = 0;

    inline void setLattice( Lattice* lattice ){ m_pLattice = lattice; }

    /// Counts how many times this process happens
    inline void eventHappened(){ m_iHappened++; }

    /// Returns how many times this process happens
    int getNumEventHappened(){ return m_iHappened; }

    /// The rules for this type of process e.g. diffusion from the site 1 to site 2 etc.
    virtual void rules( set< string, std::any > ) = 0;

protected:
    /** Pointer to the lattice of the process */
    Lattice* m_pLattice;

private:
    /// The name of this prcess
    string m_sProcName;

    /// The id of the process
    int m_iID;

    ///The type of the process
    string m_sType;

    /// Counts the times that this processes happened
    int m_iHappened;
};

#endif // PROCESS_NEW_H
