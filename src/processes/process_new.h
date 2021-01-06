#ifndef PROCESS_NEW_H
#define PROCESS_NEW_H

#include <iostream>
#include <string>
#include <map>
#include <site.h>
#include <species.h>

using namespace std;

class process_new
{
public:
    process_new();
    ~process_new();

    inline void setName( string procName ){ m_sProcName = procName; }
    inline string getName(){ return  m_sProcName; }

    inline void setID( int id ){ m_iID = id; }
    inline int getID(){ return m_iID; }

    inline void setType( int id ){ m_iID = id; }
    inline string getType( int id ){ return m_sType; }

    /// Perform this process in the site
    virtual void perform() = 0;

    inline void setLattice( Lattice* lattice ){ m_pLattice = lattice; }

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

};

#endif // PROCESS_NEW_H
