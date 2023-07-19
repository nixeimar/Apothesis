#include "adsorption.h"

namespace MicroProcesses
{

void Adsorption::simpleType()
{
    double pi = m_pUtilParams->dPi;
    double Na = m_pUtilParams->dAvogadroNum; // Avogadro's number [1/mol]
    double mass = m_dMW/Na; //[kg/mol]
    double T = m_pUtilParams->getTemperature(); //[K]
    double P = m_pUtilParams->getPressure(); //[Pa]

    m_dRateConstant = m_dStick*m_dF*P/(m_dCtot*sqrt(2.0e0*pi*mass*m_pUtilParams->dkBoltz*T) );
}

//ToDo: To be implemented and checked
void Adsorption::arrheniusType() {
    ;
}

void Adsorption::constantType(){
    m_dRateConstant = m_dAdsorptionRate*m_iNumVacant;
}

}
