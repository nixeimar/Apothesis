#include "adsorption_type.h"

void simpleType(Adsorption* proc)
{
    double pi = proc->getParameters()->dPi;
    double Na = proc->getParameters()->dAvogadroNum; // Avogadro's number [1/mol]
    double mass = proc->getMolecularWeight()/Na; //[kg/mol]
    double T = proc->getParameters()->getTemperature(); //[K]
    double P = proc->getParameters()->getPressure(); //[Pa]

    double deno = proc->getStickingCoef()*proc->getMolarFraction()*proc->getParameters()->getPressure();
    double ddeno = proc->getSitesConc()*sqrt(2.0e0*pi*mass*proc->getParameters()->dkBoltz*T);

    proc->setRateCoefficient( deno/ddeno );

//    m_dProb = m_dStick*m_dF*P/(m_dCtot*sqrt(2.0e0*pi*mass*m_pUtilParams->dkBoltz*T) );
}

//ToDo: To be implemented and checked
void arrheniusType(Adsorption* proc) {
    ;
}

void constantType(Adsorption* proc){

    proc->setRateCoefficient( proc->getAdsorptionRate()*proc->getNumVacantSites() );
//    m_dProb = m_dAdsorptionRate*m_iNumVacant;
}

