#ifndef ABSTRACT_LATTICE_H
#define ABSTRACT_LATTICE_H

#include <string>

class AbstractLattice {
public:
    
    AbstractLattice(const std::string&);
    
    virtual ~AbstractLattice(){} 

    virtual Lattice* create() = 0;
};

#endif // ABSTRACT_LATTICE_H