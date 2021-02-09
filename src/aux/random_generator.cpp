#include "random_generator.h"

namespace RandomGen {

RandomGenerator::RandomGenerator( Apothesis *apothesis ):Pointers( apothesis )
{
    m_mersenne = new CRandomMersenne( time( 0 ) );
}

RandomGenerator::~RandomGenerator() { delete m_mersenne; }

void RandomGenerator::init( const int& seed )
{
    if ( seed != 0 )
        m_mersenne->RandomInit( seed );
}

double RandomGenerator::getDoubleRandom() { return m_mersenne->Random(); }

int RandomGenerator::getIntRandom( int Min, int Max ) { return m_mersenne->IRandom( Min, Max ); }

}
