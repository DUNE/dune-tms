#ifndef _BETHEBLOCH_H_SEEN_
#define _BETHEBLOCH_H_SEEN_

#include <cmath>
#include <string>
#include <iostream>

#include "Material.h"

// Almost entirely from the PDG
// and pdg.lbl.gov/AtomicNuclearProperties
namespace BetheBloch_Utils {

  // Fine structure
  const double a = 1/137.035999139;
  // Fine structure squared
  const double a_2 = a*a;
  // Electron mass in MeV
  const double Me = 0.5109989461;
  // Squared electron mass in MeV2
  const double Me_2 = Me*Me;
  // Muon mass in MeV
  const double Mm = 105.658389;
  // Muon mass in MeV2
  const double Mm_2 = Mm*Mm;
  // Compton wavelength of electron (1E11 cm);
  const double Le = 3.8616;
  // Compton wavelength of electron2 (1E22 cm2);
  const double Le_2 = Le*Le;
  // good old pi
  const double pi = M_PI;
  // Avogadros number in 10E-23
  const double Na = 6.0222140857;

  // Sternheimer factor
  const double Stern2ln10 = 2*log(10);

  // coefficient for dE/dx in /mol cm2
  const double K = 0.307075;

  inline double RelativisticBeta(double m, double E) {
    if (E < m) return 0;

    double E2 = E*E;
    double m2 = m*m;

    return sqrt(E2-m2)/E;
  }

  inline double RelativisticGamma(double m, double E) {
    if (E < m) return 0;

    return E/m;
  }

  // Convert energy to momentum
  inline double EnergyToMomentum(double m, double E) {
    if (E < m) return 0;
    return sqrt( E*E - m*m );
  }

  // Maximum energy transfer to electron from incoming muon (causing ionisation), in MeV
  inline double MaximumEnergyTransfer(double E) {
    double me = Me;
    double me_2 = Me_2;
    double mm = Mm;
    double mm_2 = Mm_2;
    double p = EnergyToMomentum(mm, E);
    double p_2 = p*p;

    return 2 * me * p_2 / ( me_2 + mm_2 + 2*me*E );
  }
};

class BetheBloch_Calculator {

  public:

    BetheBloch_Calculator() = delete;

    BetheBloch_Calculator(Material::MaterialType type) :
      fMaterial(type) {
      };

    // Density correction factor a la PDG (Sternheimer) MeV
    // eq 34.7 pdg
    inline double DensityCorrectionFactor(double E) {
      // Update from muons in Iron from pdg
      const double x0 = fMaterial.x0;
      const double x1 = fMaterial.x1;
      const double a = fMaterial.a;
      const double m = fMaterial.m;
      const double Cbar = fMaterial.Cbar;
      const double d0 = fMaterial.d0;

      const double beta = BetheBloch_Utils::RelativisticBeta(BetheBloch_Utils::Mm, E);
      const double gamma = BetheBloch_Utils::RelativisticGamma(BetheBloch_Utils::Mm, E);
      const double X = log10( beta*gamma );

      if (x0 < X && X < x1) {
        return ( BetheBloch_Utils::Stern2ln10 * X - Cbar + a * pow(x1-X,m) );
      }
      if (X > x1) {
        return ( BetheBloch_Utils::Stern2ln10 * X - Cbar );
      }

      // Conductor
      if (X < x0) {
        return d0*pow(10,2*(X-x0));
      }

      return 0;
    }

    // Calculate the ionsiation for muons in Bethe Bloch in MeV/(gr/cm2)
    // Eq 34.5 of PDG 2020
    inline double Calc_dEdx(double E) {
      // z over A
      double Z_A = fMaterial.Z_A;
      double beta = BetheBloch_Utils::RelativisticBeta(BetheBloch_Utils::Mm, E);
      double beta_2 = beta*beta;
      double gamma = BetheBloch_Utils::RelativisticGamma(BetheBloch_Utils::Mm, E);
      double gamma_2 = gamma*gamma;
      // Ionisation from PDG
      double I = fMaterial.I * 1E-6; //286E-6; // 286 eV -> 286E-6 in MeV
      double I_2 = I*I;
      double Em = BetheBloch_Utils::MaximumEnergyTransfer(E);
      double d = DensityCorrectionFactor(E);

      double eps = (BetheBloch_Utils::K)*Z_A/beta_2; // convenient prefactor

      // K = 4pi NA re2 me c2
      double de_dx = eps *
        (0.5*log( 2*BetheBloch_Utils::Me*beta_2*gamma_2*Em/I_2 ) - beta_2 - 0.5*d );

      // Correction term in Groom not present in PDG
      //double de_dx2 = 10* a_2 * 2*pi * Na * Le_2 * Z_A * (Me / beta_2) *
      //( log( 2*Me*beta_2*gamma_2*Em/I_2 ) - 2*beta_2 + 0.25*(Em_2/E_2) - d );

      return de_dx; // in MeV/(gr/cm^2)
    }

    // Calculate the ionsiation for muons in Bethe Bloch in MeV/(gr/cm2)
    // Gaussian CLT
    // Equation 2.96 of William R Leo - Technques for Nuclear and Particle Physics
    // and ICRU Report 49
    inline double Calc_dEdx_Straggling(double E) {
      // z over A
      double Z_A = fMaterial.Z_A;
      double beta = BetheBloch_Utils::RelativisticBeta(BetheBloch_Utils::Mm, E);
      double beta_2 = beta*beta;
      double Em = BetheBloch_Utils::MaximumEnergyTransfer(E);

      double Xi = ((BetheBloch_Utils::K/2)*Z_A/beta_2);
      double kappa = Xi/Em;
      double sigma = Xi*sqrt( (1-beta_2)/(2*kappa) ); // eq 2.96 of Leo

      return sigma;
    }


    // Calculate the ionsiation for muons in Bethe Bloch MeV/(gr/cm2)
    // eq 34.12 pdg
    // Also eq 2.95 in Leo Techniques for Nuclear and Particle Physics
    inline double Calc_dEdx_mostprob(double E) {
      // z over A
      double Z_A = fMaterial.Z_A;
      double beta = BetheBloch_Utils::RelativisticBeta(BetheBloch_Utils::Mm, E);
      double beta_2 = beta*beta;
      double gamma = BetheBloch_Utils::RelativisticGamma(BetheBloch_Utils::Mm, E);
      double gamma_2 = gamma*gamma;
      // Ionisation from PDG
      double I = fMaterial.I*1E-6; //286E-6; // 286 eV -> 286E-6 in MeV
      double d = DensityCorrectionFactor(E);

      // Set some thickness
      //const double thick = 1./7.85; // thickness in g/cm2 -> try 1cm with 7.85 g/cm3 density
      //const double thick = 7.85;
      //const double thick = fMaterial.rho * 1; // 1 cm of material
      const double j = 0.200; // from pdg

      //double eps = (BetheBloch_Utils::K/2)*Z_A*thick/beta_2; // convenient
      double eps = (BetheBloch_Utils::K/2)*Z_A/beta_2; // convenient

      double de_dx = eps *
        (log( 2*BetheBloch_Utils::Me*beta_2*gamma_2/I ) + log(eps/I) + j - beta_2 - d );

      //double de_dx2 = 10* a_2 * 2*pi * Na * Le_2 * Z_A * (Me / beta_2) *
      //( log( 2*Me*beta_2*gamma_2*Em/I_2 ) - 2*beta_2 + 0.25*(Em_2/E_2) - d );

      return de_dx; // in MeV/(gr/cm^2)
    }

    Material fMaterial;
};

#endif
