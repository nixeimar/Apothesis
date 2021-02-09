#ifndef RANDOM_GENERATOR_H
#define RANDOM_GENERATOR_H

#include "apothesis.h"
#include "pointers.h"

#include <iostream>
#include "time.h"

#include "aux/randomc.h"

class CRandomMersenne;

namespace RandomGen {

class RandomGenerator : public Pointers
  {
  public:
    /// Constructor
    RandomGenerator( Apothesis* apothesis );

    void init( const int& seed );

    /// Destructor
    virtual ~RandomGenerator();

    /// Returns a random floating point number from a normal distribution between 0 and 1
    double getDoubleRandom();

    /// Returns a random integer number from the interval [Min,Max]
    int getIntRandom( int Min, int Max );

  private:
    /// The random generator used in the computations
    CRandomMersenne* m_mersenne;
  };

}

#endif
