#ifndef _MATERIAL_H_SEEN_
#define _MATERIAL_H_SEEN_

#include <string>
#include <iostream>

// Generic material class
// Only support iron and polystyrene for now
// Taken right from LBL and PDG
class Material {

  public:

    enum MaterialType { kPolyStyrene, kIron, kGraphite, kGArgon, kLArgon, kWater, kAir, 
      kUnknown };

    std::string MaterialName() {

      switch (fMaterialType) {
        case kPolyStyrene:
          return "Polystyrene";
          break;

        case kIron:
          return "Iron";
          break;

        case kGraphite:
          return "Graphite";
          break;

        case kGArgon:
          return "GAr";
          break;

        case kLArgon:
          return "LAr";
          break;

        case kWater:
          return "H_{2}O";
          break;

        case kAir:
          return "Air";
          break;

        default:
          return "unknown";
          break;
      }

      return "unknown";
    }

    // Need a constructor for the TMS materials
    // Could take strings but unnecessary string comparisons
    // Compare on density for now
    Material(double density) {
      fMaterialType = kUnknown;

      // Do an ugly check of the density to roughly ID the material to match to PDG properties
      // 1.05 density of the scintillator in TMS
      if (fabs(density-1.05) < 1E-3) {
        fMaterialType = kPolyStyrene;
        // 7.85 density of the steel TMS (maybe want other particle properties)
      } else if (fabs(density-7.85) < 1E-3) {
        fMaterialType = kIron;
        // Air in the TMS
      } else if (fabs(density-0.001225) < 1E-3) {
        fMaterialType = kAir;
      } else {
        fMaterialType = kUnknown;
      }

      SetProperties();
    }

    // The constructor with a type
    Material(MaterialType type) {
      fMaterialType = type;
      SetProperties();
    }

    // Set the properties of the material after a materialtype is chosen
    // When Z_A and Z and A may disagree, it's because PDG gives Z_A (and not Z and A separately), and Z and A are coming from the geometry
    inline void SetProperties() {

      switch (fMaterialType) {
        // Polystyrene: https://pdg.lbl.gov/2020/AtomicNuclearProperties/MUE/muE_polystyrene.pdf
        case kPolyStyrene:
          Z = 5.58;
          A = 11.0865;
          Z_A = 0.53768;
          rho = 1.060;
          I = 68.7;
          a = 0.16454;
          m = 3.2224;
          x0 = 0.1647;
          x1 = 2.5031;
          Cbar = 3.2999;
          d0 = 0.00;
          break;

          // Iron: https://pdg.lbl.gov/2020/AtomicNuclearProperties/MUE/muE_iron_Fe.pdf
        case kIron:
          // From SteelRMMS
          Z = 25.827;
          A = 55.461;
          Z_A = Z/A;
          rho = 7.874;
          I = 286.0;
          a = 0.14680;
          m = 2.9632;
          x0 = -0.0012;
          x1 = 3.1531;
          Cbar = 4.2911;
          d0 = 0.12;
          break;

          // Graphite
          // https://pdg.lbl.gov/2020/AtomicNuclearProperties/MUE/muE_carbon_graphite_C.pdf
        case kGraphite:
          Z = 6.;
          A = 12.0107;
          Z_A = Z/A;
          rho = 2.210;
          I = 78.0;
          a = 0.20762;
          m = 2.9532;
          x0 = -0.0090;
          x1 = 2.4817;
          Cbar = 2.8926;
          d0 = 0.14;
          break;

          // Gaseous Argon
          // https://pdg.lbl.gov/2020/AtomicNuclearProperties/MUE/muE_argon_gas_Ar.pdf
        case kGArgon:
          Z = 18.0;
          A = 39.948;
          Z_A = Z/A;
          rho = 1.662E-3;
          I = 188;
          a = 0.19714;
          m = 2.9618;
          x0 = 1.7635;
          x1 = 4.4855;
          Cbar = 11.9480;
          d0 = 0.00;
          break;

          // Liquid Argon
          // https://pdg.lbl.gov/2020/AtomicNuclearProperties/MUE/muE_liquid_argon.pdf
        case kLArgon:
          Z = 18.0;
          A = 39.948;
          Z_A = Z/A;
          rho = 1.396;
          I = 188;
          a = 0.19559;
          m = 3.0000;
          x0 = 0.2000;
          x1 = 3.0000;
          Cbar = 5.2146;
          d0 = 0.00;
          break;

          // Liquid water
          // https://pdg.lbl.gov/2020/AtomicNuclearProperties/MUE/muE_water_liquid.pdf
        case kWater:
          Z = 10.;
          A = 18.;
          Z_A = 0.55509;
          rho = 1.000;
          I = 79.7;
          a = 0.09116;
          m = 3.4773;
          x0 = 0.2400;
          x1 = 2.8004;
          Cbar = 3.5017;
          d0 = 0.00;
          break;

          // Air dry at 1atm
          // https://pdg.lbl.gov/2020/AtomicNuclearProperties/MUE/muE_air_dry_1_atm.pdf
        default: // TODO is this bad to do? for final production this should be reverted to below
        case kAir:
          Z = 7.31201;
          A = 14.7131;
          Z_A = 0.49919;
          rho = 1.205E-3;
          I = 85.7;
          a = 0.10914;
          m = 3.3994;
          x0 = 1.7418;
          x1 = 4.2759;
          Cbar = 10.5961;
          d0 = 0.00;
          break;

        //default:
        //  std::cerr << "Material not supported" << std::endl;
        //  std::cout << fMaterialType << std::endl;
        //  throw;
      }
    }

    double Z;
    double A;
    // Z/A
    double Z_A;
    // Ionisation
    double I; // eV
    // density
    double rho; // g/cm3

    // Numbers for density correction factors to Bethe
    double x0;
    double x1;
    double a;
    double m;
    double Cbar;
    double d0;

    MaterialType fMaterialType;
};

#endif
