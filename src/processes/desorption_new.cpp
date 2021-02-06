#include "desorption_new.h"

Desorption_new::Desorption_new():m_iNeigh(0){}

Desorption_new::~Desorption_new(){}

void Desorption_new::perform( int siteID  )
{
    m_pLattice->desorp( siteID, m_Species );
}

double Desorption_new::getProbability(){

    //These must trenafered in the global definitions
    /*--- Taken from  Lam and Vlachos (2000)PHYSICAL REVIEW B, VOLUME 64, 035401 - DOI: 10.1103/PhysRevB.64.035401 ---*/
    double Na = 6.0221417930e+23;				// Avogadro's number [1/mol]
    double P = 101325;					// [Pa]
    double T = 500;						// [K]
    double k = 1.3806503e-23;			// Boltzmann's constant [j/K]
    double s0 = 0.1;
    double C_tot = 1.0e+19;				// [sites/m^2] Vlachos code says [moles sites/m^2]
    double E_d = (7.14e+4)/Na;			// [j]
    double E = 71128/Na;   //(7.14e+4)/Na;			// [j] -> 17 kcal
    double m = 32e-3/Na;				// [kg]
    double E_m = (4.28e+4)/Na;			// [j]
    double k_d = 1.0e+13;				// [s^-1]
    double y = 2.0e-3;					// Mole fraction of the precursor on the wafera
    /*--------------------------------------------------*/

    double v0 = k_d; //*exp(-E/(k*T));
    double A = 0.0e0; //exp((E_d-E_m)/(k*T));

    return v0*exp(-m_iNeigh*E/(k*T));			//Desorption 1 neigh
}
