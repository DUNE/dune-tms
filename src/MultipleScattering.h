#ifndef _TMS_MSC_H_SEEN_
#define _TMS_MSC_H_SEEN_

#include "Material.h"
// Needed for some utils
#include "BetheBloch.h"

class MultipleScatter_Calculator {

  public:

    MultipleScatter_Calculator() = delete;
    MultipleScatter_Calculator(Material::MaterialType type) : fMaterial(type) {
      Reset();
    }

    // Implement highland scatter also?
    inline double Calc_MS_Highland() {return 0;}

    // Lynch and Dahl multiple scattererer... erer
    // E in MeV, steplength in g/cm2
    // Calculates contribution of this steplength and material to the total multiple scattering
    // Get result by calling Calc_MS_Sigma, which sums up the individual contributions
    inline void Calc_MS(double E, double steplength) {

      double p = sqrt(E*E-BetheBloch_Utils::Mm*BetheBloch_Utils::Mm);
      double p_2 = p*p;

      // Relativistic beta
      double beta = BetheBloch_Utils::RelativisticBeta(BetheBloch_Utils::Mm, E);
      double beta_2 = beta*beta;

      // Fine structure constant
      double alpha = BetheBloch_Utils::a;

      // Step length (amount of distance in this material)
      double X = steplength;

      // Material properties
      double Z = fMaterial.Z;
      double A = fMaterial.A;

      // Don't include z (charge) because it always comes in as z2 and unnecessary

      // eq 1, 2 of Lynch and Dahl
      // NB these are not chi2 in the statistical sense!
      chisq_c += 0.157*(Z*(Z+1)*X/A)/(p_2*beta_2);
      double chisq_a = 2.007E-5*pow(Z, 2./3.)*(1+3.34*(Z*Z*alpha*alpha/beta_2))/(p_2);

      // Sum over the scatterers when we have multiple materials between the scintillators
      // Increment this in eq 7, 8, 9 by following eq 11
      double weight = X*Z*(Z+1)/A;
      numerator += weight*log(sqrt(chisq_a)); // top of equation 11 (the numerator)
      denom += weight; // bottom of equation 11 (the denominator)
    }

    // Finally calculate equation 7,8,9
    inline double Calc_MS_Sigma() {
      // Once the calculation is finished we can calculate omega
      const double F = 0.98; // Can probably play with this to see effect: comes from page 7 in section 3 (just before section 4 begins)

      double chi_a_eff = exp(numerator/denom); // equation 11 of L&D (solve for chi_alpha
      double chisq_a_eff = chi_a_eff*chi_a_eff; // square it 
      double omega = chisq_c/chisq_a_eff; // Use the effective
      double v = 0.5*omega/(1-F);
      double sigma = chisq_c/(1+F*F)*(((1+v)/v)*log1p(v)-1); // log1p(x) is just a more accurate version of log(1+x), i.e. log1p(v) is doing log(v+1)

      // Reset at the end?
      Reset();

      return sigma;
    }

    inline void Reset() {
      chisq_c = 0.;
      numerator = 0.;
      denom = 0.;
    }

    Material fMaterial;

  private:
    double chisq_c;
    double numerator;
    double denom;

};

#endif
