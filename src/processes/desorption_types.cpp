#include "desorption.h"

namespace MicroProcesses
{

void Desorption::constantType(){
    m_dRateConstant = m_dDesorptionRate*m_iNumVacant;
}

void Desorption::arrheniusType()
{
    double T = m_pUtilParams->getTemperature();
    double k = m_pUtilParams->dkBoltz;
    double Ed = m_dEd/m_pUtilParams->dAvogadroNum;

    m_dRateConstant = m_dv0*exp(-(double)(m_iNumNeighs + 1)*Ed/(k*T));
}

}
