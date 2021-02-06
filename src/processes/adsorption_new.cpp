#include "adsorption_new.h"

Adsorption_new::Adsorption_new():m_Species(0){}

Adsorption_new::~Adsorption_new(){}

void Adsorption_new::perform( int siteID )
{
    m_pLattice->adsorp( siteID, m_Species );
}

/*--- Taken from  Lam and Vlachos (2000)PHYSICAL REVIEW B, VOLUME 64, 035401 - DOI: 10.1103/PhysRevB.64.035401

        double Na = 6.0221417930e+23;				// Avogadro's number [1/mol]
        double P = 101325;					// [Pa]
        double T = 500;						// [K]
        double k = 1.3806503e-23;			// Boltzmann's constant [j/K]
        double s0 = 0.1;
        double C_tot = 1.0e+19;				// [sites/m^2] Vlachos code says [moles sites/m^2]
        double E_d = (7.14e+4)/Na;			// [j]
        double E = 71128/Na;   //(7.14e+4)/Na;			// [j]
        double m = 32e-3/Na;				// [kg]
        double E_m = (4.28e+4)/Na;			// [j]
        double k_d = 1.0e+13;				// [s^-1]
        double y = 2.0e-3;					// Mole fraction of the precursor on the wafer


//pa = s0*P*y/(C_tot*sqrt(2.0e0*3.14159265*m*k*T));

/*-------------------------------------------------

        double v0 = k_d; //*exp(-E/(k*T));
        double A = 0.0e0; //exp((E_d-E_m)/(k*T));

        pa = s0*P*y/(C_tot*sqrt(2.0e0*3.14159265*m*k*T));

//--------------------- Transitions probabilities ----------------------------------------//
        (*prob)[0] = pa*Nx*Ny;									//Adsorption
        (*prob)[1] = group[0].size()*v0*exp(-1.0e0*E/(k*T));			//Desorption 1 neigh
        (*prob)[3] = group[1].size()*v0*exp(-2.0e0*E/(k*T));        //Desorption 2 neigh
        (*prob)[5] = group[2].size()*v0*exp(-3.0e0*E/(k*T));		//Desorption 3 neigh
        (*prob)[7] = group[3].size()*v0*exp(-4.0e0*E/(k*T));		//Desorption 4 neigh
        (*prob)[9] = group[4].size()*v0*exp(-5.0e0*E/(k*T));		//Desorption 5 neigh
        (*prob)[2] = A*(*prob)[1];	  							//Diffusion  1 neisgh
        (*prob)[4] = A*(*prob)[3];								//Diffusion  2 neisgh
        (*prob)[6] = A*(*prob)[5];								//Diffusion  3 neisgh
        (*prob)[8] = A*(*prob)[7];								//Diffusion  4 neisgh
        (*prob)[10] = A*(*prob)[9];								//Diffusion  5 neisgh
//----------------------------------------------------------------------------------------//

        (*p_tot) = 0;
        for (unsigned int i=0; i<nof_trans_prob; i++)
                 *p_tot += (*prob)[i];

/*        for (unsigned int i=0; i<nof_trans_prob; i++)
                cout<< "NO"<< "\t" <<(*prob)[i] << endl;

        cout<<"***************"<< endl;*/

//---------Canonical-form--------------//
//       for (unsigned int i=0; i<nof_trans_prob; i++)
//               (*prob)[i] /= (*p_tot);
//-------------------------------------//

/*	for (unsigned int i=0; i<nof_trans_prob; i++)
                cout<<"CAN"<< "\t" <<(*prob)[i] << endl;
        cout<<"***************"<< endl;
        cout<<"PROBS"<<endl;
        system("pause");
--------------------------------------------------*/

double Adsorption_new::getProbability(){

    //These must trenafered in the global definitions
    double Na = 6.0221417930e+23;				// Avogadro's number [1/mol]
    double P = 101325;					// [Pa]
    double T = 500;						// [K]
    double k = 1.3806503e-23;			// Boltzmann's constant [j/K]
    double s0 = 0.1;
    double C_tot = 1.0e+19;				// [sites/m^2] Vlachos code says [moles sites/m^2]
    double E_d = (7.14e+4)/Na;			// [j]
    double E = 71128/Na;   //(7.14e+4)/Na;			// [j]
    double m = 32e-3/Na;				// [kg]
    double E_m = (4.28e+4)/Na;			// [j]
    double k_d = 1.0e+13;				// [s^-1]
    double y = 2.0e-3;					// Mole fraction of the precursor on the wafer

    return s0*y*P/(2.*C_tot*sqrt(2.0e0*3.14159265*m*k*T) );
}
