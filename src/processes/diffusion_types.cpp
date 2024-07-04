#include "diffusion_types.h"

namespace MicroProcesses
{

double constantType(Diffusion* proc){
    return proc->getDiffusionRate()*proc->getNumVacantSites();
}

double arrheniusType(Diffusion* proc)
{
    /*--- Taken from  Lam and Vlachos (2000)PHYSICAL REVIEW B, VOLUME 64, 035401 - DOI: 10.1103/PhysRevB.64.035401 ---*/
    /*    double Na = 6.0221417930e+23;				// Avogadro's number [1/mol]
    double P = 101325;					// [Pa]
    double T = any_cast<double>(m_vParams[0]); //500;						// [K]
    double k = any_cast<double>(m_vParams[1]); // 1.3806503e-23;			// Boltzmann's constant [j/K]
    double s0 = 0.1;
    double C_tot = 1.0e+19;				// [sites/m^2] Vlachos code says [moles sites/m^2]
    double E_d = any_cast<double>(m_vParams[2]); //(7.14e+4)/Na;			// [j]
    double E = 71128/Na;   //(7.14e+4)/Na;			// [j] -> 17 kcal
    double m = 32e-3/Na;				// [kg]
    double E_m = any_cast<double>(m_vParams[3]); //(4.28e+4)/Na;			// [j]
    double k_d = 1.0e+13;				// [s^-1]
    double y = 2.0e-3;					// Mole fraction of the precursor on the wafer
    /*--------------------------------------------------*/

    //   double v0 = k_d; //*exp(-E/(k*T));
    //   double A = exp( (E_d-E_m)/(k*T) );

    //--------------------- Transitions probability ----------------------------------------//
    //  return 0;// A*v0*exp( -(double)any_cast<int>(m_mParams["neighs"])*E/(k*T) );
    //----------------------------------------------------------------------------------------//

    double v0 = proc->getVibrationalFrequency();
    double E = proc->getActivationEnergy();
    double Em = proc->getDifActivationEnergy();
    double T = proc->getParameters()->getTemperature();
    int n = proc->getNumNeighs()+1;

    double k = proc->getParameters()->dkBoltz;
    E = E/proc->getParameters()->dAvogadroNum;
    Em = Em/proc->getParameters()->dAvogadroNum;
    double A = exp(E-Em)/(k*T);

    return v0*A*exp(-(double)n*E/(k*T));
}

}
