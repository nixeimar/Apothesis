#include "adsorption_types.h"

namespace MicroProcesses
{

double simpleType(Adsorption* proc)
{
    double pi = proc->getParameters()->dPi;
    double Na = proc->getParameters()->dAvogadroNum; // Avogadro's number [1/mol]
    double mass = proc->getMolecularWeight()/Na; //[kg/mol]
    double T = proc->getParameters()->getTemperature(); //[K]
    double P = proc->getParameters()->getPressure(); //[Pa]
    double f = proc->getMolarFraction();
    double s = proc->getStickingCoef();

    double ctot = proc->getSitesConc();
    double kb = proc->getParameters()->dkBoltz;

    double numerator = s*f*P;
    double denumerator = ctot*sqrt(2.0e0*pi*mass*kb*T);

    return numerator/denumerator;
}

//ToDo: To be implemented and checked
double arrheniusType(Adsorption* proc) {
    return 0.0;
}

double constantType(Adsorption* proc) {
    return proc->getAdsorptionRate()*proc->getNumVacantSites();
}

}
