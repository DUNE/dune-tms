#include "TMS_Kalman.h"

TMS_Kalman::TMS_Kalman() : 
  Bethe(Material::kPolyStyrene),
  MSC(Material::kPolyStyrene),
  ForwardFitting(false) {
}

// Take a collection of hits
TMS_Kalman::TMS_Kalman(std::vector<TMS_Hit> &Candidates) : 
  Bethe(Material::kPolyStyrene), 
  MSC(Material::kPolyStyrene),
  //ForwardFitting(false),
  ForwardFitting(true),
  Talk(false)
  //Talk(true)
{

  // Empty the KalmanStates
  KalmanNodes.clear();

  // Save the number of initial candidates
  int nCand = Candidates.size();
  // And muon mass
  mass = BetheBloch_Utils::Mm;

  std::cout << "Front: " << Candidates.front().GetZ() << std::endl;
  if (ForwardFitting) std::sort(Candidates.begin(), Candidates.end(), TMS_Hit::SortByZInc);
  else                std::sort(Candidates.begin(), Candidates.end(), TMS_Hit::SortByZ);
  std::cout << "After sort front: " << Candidates.front().GetZ() << std::endl;

  // Make a new Kalman state for each hit
  KalmanNodes.reserve(nCand);
  for (int i = 0; i < nCand; ++i) {

    TMS_Hit hit = Candidates[i];
    //if (hit.GetBar().GetBarType() != TMS_Bar::kYBar) continue;
    double x = hit.GetNotZ();
    double y = hit.GetTrueHit().GetY();
    double z = hit.GetZ();
    double future_z = (i+1 == nCand ) ? z : Candidates[i+1].GetZ();

    double DeltaZ = future_z-z;
    std::cout << "hit: " <<  i << " z: " << z << " deltaz: " << DeltaZ << std::endl;

    // This also initialises the state vectors in each of the nodes
    TMS_KalmanNode Node(x, y, z, DeltaZ);
    KalmanNodes.emplace_back(std::move(Node));
  }

  // Set the momentum seed for the first hit from its length
  if (ForwardFitting) {
    double startx = Candidates.front().GetNotZ();
    double endx = Candidates.back().GetNotZ();
    double startz = Candidates.front().GetZ();
    double endz = Candidates.back().GetZ();
    double KEest = GetKEEstimateFromLength(startx, endx, startz, endz);
    double momest = sqrt((KEest+mass)*(KEest+mass)-mass*mass);
    std::cout << "momentum estimate from length: " << momest << std::endl;
    KalmanNodes.front().PreviousState.qp = 1./momest;
    KalmanNodes.front().CurrentState.qp = 1./momest;
  }

  RunKalman();
}

// Used for seeding the starting KE for the Kalman filter from the start and end point of a track
double TMS_Kalman::GetKEEstimateFromLength(double startx, double endx, double startz, double endz) {
  // if in thin and thick target there's a different relationship
  double distx = endx-startx;
  double distz = endz-startz;
  double dist = sqrt(distx*distx+distz*distz);

  /* pol0 fit to xz distance of start and death point using truth, for track ending in THICK (death > 739 cm in z)
     p0                        =      101.506   +/-   4.05823     
     p1                        =     0.132685   +/-   0.00354222  
     */

  /* pol0 fit to xz distance of start and death point using truth, for track ending in THIN (death < 739 cm in z)
     p0                        =     -1.13305   +/-   0.129217    
     p1                        =     0.234404   +/-   0.00218364  
     */

  double KEest = 0;
  // If in thick region
  if (endz > TMS_Const::TMS_Thick_Start) KEest = 101.5+0.133*dist;
  else KEest = -1.13+0.234*dist;

  std::cout << "dist: " << dist << " KE: " << KEest << std::endl;

  return KEest;
}

// Update to the next step
void TMS_Kalman::RunKalman() {

  int nCand = KalmanNodes.size();
  for (int i = 1; i < nCand; ++i) {
    // Perform the update from the (i-1)th node's predicted to the ith node's previous
    Update(KalmanNodes[i-1], KalmanNodes[i]);
    Predict(KalmanNodes[i]);
  }
  std::cout << "Momentum of last hit: " << 1./KalmanNodes.back().CurrentState.qp << std::endl;
}

void TMS_Kalman::Update(TMS_KalmanNode &PreviousNode, TMS_KalmanNode &CurrentNode) {
  CurrentNode.PreviousState = PreviousNode.CurrentState;
}

// Predict the next step
// Here we use the previous state to predict the current state
// We also calculate the updated noise matrix from multiple scattering and energy loss
// Account for energy loss, multiple scattering, and bending due to the magnetic field
void TMS_Kalman::Predict(TMS_KalmanNode &Node) {

  // Right out of MINOS
  // Can probably do this better; e.g. calculate noise matrix at the same time as updating others (avoid duplicate calculations)
  /*
     BuildNoiseMatrix();
     ExtrapolateCovariance();
     KalmanGain();
     UpdateStateVector();
     UpdateCovariance();
     MoveArrays();
     */

  TMS_KalmanState &PreviousState = Node.PreviousState;
  TMS_KalmanState &CurrentState = Node.CurrentState;

  // Propagate the current state
  TMatrixD &Transfer = Node.TransferMatrix;
  std::cout << "Transfer matrix: " << std::endl;
  Transfer.Print();
  TVectorD PreviousVec(5);
  PreviousVec[0] = PreviousState.x;
  PreviousVec[1] = PreviousState.y;
  PreviousVec[2] = PreviousState.dxdz;
  PreviousVec[3] = PreviousState.dydz;
  PreviousVec[4] = PreviousState.qp;
  // Just the propagator matrix influence
  TVectorD UpdateVec = Transfer*(PreviousVec);
  // Now construct the current state (z of CurrentState is already set to be z+dz)
  CurrentState.x = UpdateVec[0];
  CurrentState.y = UpdateVec[1];
  CurrentState.dxdz = UpdateVec[2];
  CurrentState.dydz = UpdateVec[3];
  // Don't update q/p until later (when we've done the energy loss calculation)


  std::cout << "Previous vector: " << std::endl;
  PreviousState.Print();
  std::cout << "Current vector (before energy loss): " << std::endl;
  CurrentState.Print();


  // Update the energy
  double mom = 1./PreviousState.qp;
  // Initial energy before energy loss
  double en_initial = sqrt(mom*mom+mass*mass);
  // The energy we'll be changing
  double en = en_initial;

  if (Talk) std::cout << "mom: " << mom << std::endl;

  // Read the position between current point and extrapolated into next bar
  double xval = PreviousState.x;
  double yval = PreviousState.y;
  double zval = PreviousState.z;

  double xval2 = CurrentState.x;
  double yval2 = CurrentState.y;
  double zval2 = CurrentState.z;

  // Make TVector3s of the two points
  TVector3 start(xval,yval,zval); // Start
  TVector3 stop(xval2,yval2,zval2); // Stop

  if (Talk) {
    std::cout << "Going from " << start.X() << " " << start.Y() << " " << start.Z() << std::endl;
    std::cout << "To " << stop.X() << " " << stop.Y() << " " << stop.Z() << std::endl;
  }

  // Get the materials between the two points
  std::vector<std::pair<TGeoMaterial*, double> > Materials = TMS_Geom::GetInstance().GetMaterials(start, stop);

  if (Talk) std::cout << "Looping over " << Materials.size() << " materials" << std::endl;
  double TotalPathLength = 0;
  double TotalLength = 0;

  // Loop over the materials between the two projection points
  int counter = 0;
  double total_en_var = 0;
  if (Talk) std::cout << "Energy before: " << en_initial << std::endl;
  for (auto material : Materials) {

    // Read these directly from a TGeoManager
    // If the geometry is in mm (CLHEP units), then want to scale density to g/cm3 and thickness to cm
    // Otherwise, assume it's like that and then fix it with geometry scaling functions
    double density = material.first->GetDensity()/(CLHEP::g/CLHEP::cm3); 
    double thickness = material.second/10.; 
    // Potentially need to scale from g/cm3 to g/mm3, so find the scale factor and scale by 1/that^3.
    double scale_factor = TMS_Geom::GetInstance().Scale(1.0);
    density /= std::pow(scale_factor, 3);
    thickness = TMS_Geom::GetInstance().Scale(thickness);
    material.first->Print();

    if (Talk) {
      std::cout << "Material " << counter << " = " << material.first->GetName() << std::endl;
      std::cout << "  density: " << density << std::endl;
      std::cout << "  thickness: " << thickness << std::endl;
      std::cout << "  thickness*density = " << density*thickness << std::endl;
    }

    // Skip if density or thickness is small
    if (density*thickness < 0.1) {
      if (Talk) std::cout << "  Skipping material, to little path length to bother" << std::endl;
      continue;
    }

    TotalPathLength += density*thickness;
    TotalLength += thickness;

    // Update the Bethe Bloch calculator to use this material
    Material matter(density);

    Bethe.fMaterial = matter;

    if (Talk) {
      std::cout << "Mean: " << Bethe.Calc_dEdx(en)*density*thickness << std::endl;
      std::cout << "Most prob: " << Bethe.Calc_dEdx_mostprob(en)*density*thickness << std::endl;
      std::cout << "energy before: " << en << std::endl;
    }

    // Subtract off the energy loss for this material
    if (ForwardFitting) en -= Bethe.Calc_dEdx(en)*density*thickness;
    else                en += Bethe.Calc_dEdx(en)*density*thickness;
    if (Talk) std::cout << "energy after: " << en << std::endl;


    // Variance assuming Gaussian straggling
    double en_var = Bethe.Calc_dEdx_Straggling(en)*density*thickness;
    total_en_var += en_var*en_var;
    if (Talk) std::cout << "energy var: " << en_var*en_var << std::endl;

    // Set the material for the multiple scattering
    MSC.fMaterial = matter;
    // Calculate this before or after the energy subtraction/addition?
    MSC.Calc_MS(en, thickness*density);

    counter++;
  }
  if (Talk) {
    std::cout << "Total path length (g/cm2): " << TotalPathLength << std::endl;
    std::cout << "Total length (cm): " << TotalLength << std::endl;
  }

  // Updated momentum^2
  double p_2_up = en*en-BetheBloch_Utils::Mm*BetheBloch_Utils::Mm;
  double p_up;
  if (p_2_up > 0) p_up = sqrt(p_2_up);
  else {
    std::cerr << "negative momentum squared, setting momentum to 1 MeV" << std::endl;
    p_up = 1;
  }
  std::cout << "momentum: " << p_up << std::endl;

  // Update the state's q/p
  CurrentState.qp = 1./p_up;

  std::cout << "Current vector (after energy loss): " << std::endl;
  CurrentState.Print();

  if (Talk) {
    std::cout << "total energy variation: " << total_en_var << std::endl;
    std::cout << "energy after: " << en << std::endl;
    std::cout << "Delta(Energy): " << en-en_initial << std::endl;
  }

  double p_var = (2*en/p_up)*(2*en/p_up) * total_en_var;
  double qp_var = 1./(p_up*p_up*p_up*p_up) * p_var;

  if (Talk) {
    std::cout << "momentum after: " << p_up << std::endl;
    std::cout << "total momentum variation: " << p_var << std::endl;
  }
  // Momentum variance

  // Set the noise matrix
  TMatrixD &NoiseMatrix = Node.NoiseMatrix;

  NoiseMatrix(4,4) = qp_var;

  // Finally when we've iterated through the materials we can get the combined effect from the multiple-scattering
  double ms = MSC.Calc_MS_Sigma();

  if (Talk) std::cout << "Multiple scattering: " << ms << std::endl;

  // Now proceed with Wolin and Ho (Nucl Inst A329 1993 493-500)
  // covariance for multiple scattering
  // Also see MINOS note on Kalman filter (John Marshall, Nov 15 2005)
  double ax2 = CurrentState.dxdz; // get the unit vectors
  double ay2 = CurrentState.dydz; // get the unit vectors

  double norm = 1+ax2+ay2; // 1+P3^2+P4^2 in eq 16, 17, 18 in Wolin and Ho
  double covAxAx = norm*ms*(1+ax2);// eq 16 Wolin and Ho
  double covAyAy = norm*ms*(1+ay2);// eq 17 Wolin and Ho
  double covAxAy = norm*ms*ax2*ay2;// eq 18 Wolin and Ho

  // Check the z0 variable!
  if (!ForwardFitting) TotalPathLength = -1*TotalPathLength;

  // Build the covariance matrix in Wolin and Ho after eq 18
  // Equation 15 in Robert Harr Calculation of Track and Vertex Errors for Detector Design Studies
  // Also equation 48; kappa in Harr is "norm" here
  // See also Rainer Mankel, 1998 "ranger a pattern recognition algorithm for the HERA-B main tracking system, part iv the object-oriented track fit" Appendix B, HERA-B note 98-079
  double TotalPathLengthSq = TotalPathLength*TotalPathLength;
  NoiseMatrix(0,0) = covAxAx * TotalPathLengthSq / 3.;
  NoiseMatrix(1,1) = covAyAy * TotalPathLengthSq / 3.;
  NoiseMatrix(2,2) = covAxAx;
  NoiseMatrix(3,3) = covAyAy;

  // Negative signs depend on if we're doing backward or forward fitting (- sign for decreasing z, + sign for increasing z)
  int Sign = +1;
  if (!ForwardFitting) Sign *= -1;

  NoiseMatrix(1,0) = NoiseMatrix(0,1) =    covAxAy * TotalPathLengthSq/3.;
  NoiseMatrix(2,0) = NoiseMatrix(0,2) = (Sign)*covAxAx * TotalPathLength/2.;
  NoiseMatrix(3,0) = NoiseMatrix(0,3) = (Sign)*covAxAy * TotalPathLength/2.;

  NoiseMatrix(2,1) = NoiseMatrix(1,2) = (Sign)*covAxAy * TotalPathLength/2.;
  NoiseMatrix(3,1) = NoiseMatrix(3,1) = (Sign)*covAyAy * TotalPathLength/2.;

  NoiseMatrix(3,2) = NoiseMatrix(2,3) =    covAxAy;

  // I think that's it
  // Other than the B-field...!

}


