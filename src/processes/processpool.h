#ifndef PRORCESSPOOL_H
#define PRORCESSPOOL_H

#include <string>
#include <map>

#include "process_new.h"

using namespace std;

namespace newDesign
{

class ProrcessPool
{
public:
    ProrcessPool();

    inline Process_new* getProcessByName(string name) { return m_mapProcs[ name ]; }
    inline Process_new* getProcessByID( int id ){ return m_mapProcsIDs[ id ]; }

    void addProcess( string name, Process_new* proc);
    void addProcess( int id, Process_new* proc);
    inline int getProcessNum() { return m_mapProcs.size(); }

private:
    map<string, Process_new*> m_mapProcs;
    map<int, Process_new*> m_mapProcsIDs;
};

}

#endif // PRORCESSPOOL_H
