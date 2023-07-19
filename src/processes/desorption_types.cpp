#include "desorption_types.h"

namespace MicroProcesses
{

double constantType(Desorption* proc){
    return  proc->getDesorptionRate()*proc->getNumVacantSites();
}

double arrheniusType(Desorption* proc)
{
    double T = proc->getParameters()->getTemperature();
    double k = proc->getParameters()->dkBoltz;
    double Ed = proc->getActivationEnergy()/proc->getParameters()->dAvogadroNum;

    return proc->getVibrationalFrequency()*exp(-(double)(   proc->getNumNeighs() + 1)*Ed/(k*T));
}

}
