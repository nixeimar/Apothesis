#include "macroprocess.h"

Macroprocess::Macroprocess()
{
    apothesis = new Apothesis();
}

Macroprocess::~Macroprocess() {
    if (apothesis)
        delete apothesis;

    if (parameters)
        delete parameters;
}

void Macroprocess::setIO(IO *newIO)
{
    io = newIO;
}

void Macroprocess::perform(){

    cout << "Apothesis runnning ..." << endl;
    apothesis->exec();
    cout << "Apothesis finished succesfully." << endl;
}
