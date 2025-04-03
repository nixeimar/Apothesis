#ifndef MACROPROCESS_H
#define MACROPROCESS_H

#include "apothesis.h"
#include "parameters.h"
#include "io.h"

class Macroprocess
{
public:
    Macroprocess();
    virtual ~Macroprocess();

    virtual void init() = 0;
    virtual void perform();

    void setIO( IO* io);

protected:

    Apothesis* apothesis;
    Parameters* parameters;
    IO* io;

};

#endif // MACROPROCESS_H
