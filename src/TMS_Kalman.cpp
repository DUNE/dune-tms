#include "TMS_Kalman.h"
#include "TMS_Geom.h"

TMS_Kalman::TMS_Kalman() : 
  Bethe(Material::kPolyStyrene),
  MSC(Material::kPolyStyrene),
  ForwardFitting(false) {
}

// Take a collection of hits, return collection of Kalman nodes
TMS_Kalman::TMS_Kalman(std::vector<TMS_Hit> &Candidates) : 
  Bethe(Material::kPolyStyrene), 
  MSC(Material::kPolyStyrene),
  ForwardFitting(false),
  Talk(false)
{
  TRandom3* RNG = new TRandom3(1337); // TODO: Seed properly sometime

  // Empty the KalmanStates
  KalmanNodes.clear();

  // Save the number of initial candidates
  int nCand = Candidates.size();
  // And muon mass
  mass = BetheBloch_Utils::Mm;

  std::sort(Candidates.begin(), Candidates.end(), TMS_Hit::SortByZ);

  // Make a new Kalman state for each hit
  //KalmanNodes.reserve(nCand);
  for (int i = 0; i < nCand; ++i) {

    TMS_Hit hit = Candidates[i];
    double x_true = hit.GetTrueHit().GetX();
    double y_true = hit.GetTrueHit().GetY();
    double x = hit.GetRecoX();
    double y = hit.GetRecoY();
    double z = hit.GetZ();

    TVector3 vecc = TVector3(x,y,z);

    //if (! (TMS_Geom::GetInstance().IsInsideBox(vecc, TMS_Const::TMS_Start_Exact, TMS_Const::TMS_End_Exact))) // TODO should use exact, uncomment when geometry is fixed-ish
    if (! (TMS_Geom::GetInstance().IsInsideBox(vecc, TMS_Const::TMS_Start, TMS_Const::TMS_End)))
    {
      std::cerr << "[TMS_Kalman.cpp] Hit " << i << "/" << nCand << " position not within TMS before Kalman filter, (x,y,z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
      std::cerr << "[TMS_Kalman.cpp] Were reco x and y values set before running Kalman?" << std::endl;
      //throw; // yeet it
      ErrorDetVol = true;
    }

    int j;
    for (j=i; j<nCand; j++)
      if (abs(Candidates[j].GetZ() - z) > 1E-3)
      {
        break;
      }

    double future_z = (i+1 == nCand ) ? z : Candidates[i+1].GetZ();
    double DeltaZ = future_z-z;

    // This also initialises the state vectors in each of the nodes
    if (abs(DeltaZ) > 1E-3) // TODO: Only add one hit per z for now, noise breaks
    {
      // TODO: Combine multiple hits into a single 'node' <-> 'measurement'
      TMS_KalmanNode Node(x, y, z, DeltaZ);
      Node.SetTrueXY(x_true, y_true); // Add truth to enable reco to truth comparison
      Node.LayerOrientation = hit.GetBar().GetBarType();
      Node.LayerBarWidth    = hit.GetBar().GetBarWidth();
      Node.LayerBarLength   = hit.GetBar().GetBarLength();

      KalmanNodes.emplace_back(std::move(Node));
//    } else { // TODO: Handle layers with more than one hit, waiting on Asa to confirm potential structures
//      //std::cout << "more than one hit per layer? Kalman unhappy " << i << "\t " << j-i << std::endl;
//      for (int k = i; k<j; k++)
//      {
//      }
    }
  }

  int N_LAYER_BACK = 10;
  // Can't look back further than the first element
  if (Candidates.size() < (unsigned)N_LAYER_BACK)
    N_LAYER_BACK = Candidates.size();

  AverageXSlope = (Candidates[Candidates.size() - N_LAYER_BACK].GetRecoX() - Candidates.back().GetRecoX()) / (Candidates[Candidates.size() - N_LAYER_BACK].GetZ() - Candidates.back().GetZ());
  AverageYSlope = (Candidates.front().GetRecoY() - Candidates.back().GetRecoY()) / (Candidates.front().GetZ() - Candidates.back().GetZ());

  // Set the momentum seed for the first hit from its length
  if (ForwardFitting) {
    double startx = Candidates.front().GetNotZ();
    double endx = Candidates.back().GetNotZ();
    double startz = Candidates.front().GetZ();
    double endz = Candidates.back().GetZ();
    double KEest = GetKEEstimateFromLength(startx, endx, startz, endz);
    double momest = sqrt((KEest+mass)*(KEest+mass)-mass*mass);
    if (Talk) std::cout << "momentum estimate from length: " << momest << std::endl;
    KalmanNodes.front().PreviousState.qp = 1./momest;
    KalmanNodes.front().CurrentState.qp = 1./momest;
    KalmanNodes.back().CurrentState.dxdz = 0.0;//AverageXSlope;
    KalmanNodes.back().CurrentState.dydz = 0.0;//AverageYSlope;
    KalmanNodes.back().PreviousState.dxdz = 0.0;//AverageXSlope;
    KalmanNodes.back().PreviousState.dydz = 0.0;//AverageYSlope;
  } else { // TODO check if 0 is sane
    KalmanNodes.back().CurrentState.dxdz = 0.0;//AverageXSlope;
    KalmanNodes.back().CurrentState.dydz = 0.0;//AverageYSlope;
    KalmanNodes.back().PreviousState.dxdz = 0.0;//AverageXSlope;
    KalmanNodes.back().PreviousState.dydz = 0.0;//AverageYSlope;
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

  SetMomentum(1./KalmanNodes.back().CurrentState.qp);

  // Set start pos/dir
  SetStartPosition(KalmanNodes.back().CurrentState.x, KalmanNodes.back().CurrentState.y, KalmanNodes.back().CurrentState.z);
  SetStartDirection(KalmanNodes.back().CurrentState.dxdz, KalmanNodes.back().CurrentState.dydz);
  // Set end pos/dir
  if (KalmanNodes.size() > 1) {
    SetEndPosition(KalmanNodes.at(1).CurrentState.x, KalmanNodes.at(1).CurrentState.y, KalmanNodes.at(1).CurrentState.z);
    SetEndDirection(KalmanNodes.at(1).CurrentState.dxdz, KalmanNodes.at(1).CurrentState.dydz);
  } else { // Kalman output is rubbish in this case, but we don't crash :)
    SetEndPosition(KalmanNodes.front().CurrentState.x, KalmanNodes.front().CurrentState.y, KalmanNodes.front().CurrentState.z);
    SetEndDirection(KalmanNodes.front().CurrentState.dxdz, KalmanNodes.front().CurrentState.dydz);
  }

  if (std::isnan(momentum) || std::isinf(momentum)){
    std::cerr << "[TMS_Kalmann.cpp] Weirdness -- Momentum from fitter isn't a sane number: " << momentum << std::endl;
  }
}

void TMS_Kalman::Update(TMS_KalmanNode &PreviousNode, TMS_KalmanNode &CurrentNode) {
  CurrentNode.PreviousState = PreviousNode.CurrentState;
  CurrentNode.CovarianceMatrix = PreviousNode.CovarianceMatrix;
}

// Predict the next step
// Here we use the previous state to predict the current state
// We also calculate the updated noise matrix from multiple scattering and energy loss
//jAccount for energy loss, multiple scattering, and bending due to the magnetic field
void TMS_Kalman::Predict(TMS_KalmanNode &Node) {

  TMS_KalmanState &PreviousState = Node.PreviousState;
  TMS_KalmanState &CurrentState = Node.CurrentState;

  // Matrices to propagate the current state
  TMatrixD &Transfer = Node.TransferMatrix;
  TMatrixD &TransferT = Node.TransferMatrixT;

  // Initialise to something sane(-ish)
  if (PreviousState.dxdz ==  -999.9) // Only on initialisation?
  {
    //PreviousState.dxdz = TMS_Kalman::AverageXSlope;
    if ( abs(TMS_Kalman::AverageXSlope) > 2.0 )
    {
      std::cerr << "[TMS_Kalman.cpp] Excessive average X slope = " << TMS_Kalman::AverageXSlope << " of track (first to last hit), setting to 0" << std::endl;
      PreviousState.dxdz = 0.0;
    } else {
      PreviousState.dxdz = TMS_Kalman::AverageXSlope;
    }
  }
  if (PreviousState.dydz ==  -999.9) // Only on initialisation?
  {
    if ( abs(TMS_Kalman::AverageYSlope) > 1.5 )
    {
      std::cerr << "[TMS_Kalman.cpp] Excessive average Y slope = " << TMS_Kalman::AverageYSlope << " of track (first to last hit), setting to 0" << std::endl;
      PreviousState.dydz = 0.0;
    } else {
      PreviousState.dydz = TMS_Kalman::AverageYSlope;
    }
  }

  TVectorD PreviousVec(5);
  PreviousVec[0] = PreviousState.x;
  PreviousVec[1] = PreviousState.y;
  PreviousVec[2] = PreviousState.dxdz;
  PreviousVec[3] = PreviousState.dydz;
  PreviousVec[4] = PreviousState.qp;

  TVectorD UpdateVec = Transfer*(PreviousVec);

  if (Talk) std::cout << "\nPrevious vector: " << std::endl;
  if (Talk) PreviousState.Print();

  // ENERGY LOSS PART
  // Update the energy
  double mom = 1./PreviousState.qp;
  //if (std::isinf(mom)) mom = 4800; // set to 1 GeV
  // Initial energy before energy loss
  double en_initial = sqrt(mom*mom+mass*mass);
  // The energy we'll be changing
  double en = en_initial;

  // Read the position between current point and extrapolated into next bar
  double xval = PreviousState.x;
  double yval = PreviousState.y;
  double zval = PreviousState.z;

  double xval2 = CurrentState.x;
  double yval2 = CurrentState.y;
  double zval2 = PreviousState.z + Transfer(0,2); // Probably a nicer way to do this (:

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
  if (Talk) std::cout << "Energy before " << Materials.size() << " materials: " << en_initial << std::endl;
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

    TotalPathLength += density*thickness;
    TotalLength += thickness;

    // Update the Bethe Bloch calculator to use this material
    try {
      Material matter(density);
      Bethe.fMaterial = matter;
      // Set the material for the multiple scattering
      MSC.fMaterial = matter;
    }
    catch (std::invalid_argument const& ex) {
      std::cout<<"Could not make a material using density "<<density<<", is position within tms bounds?"<<std::endl;
    }

    // Subtract off the energy loss for this material
    if (ForwardFitting) en -= Bethe.Calc_dEdx(en)*density*thickness;
    else                en += Bethe.Calc_dEdx(en)*density*thickness;


    // Variance assuming Gaussian straggling
    double en_var = Bethe.Calc_dEdx_Straggling(en)*density*thickness;
    total_en_var += en_var*en_var;

    // Calculate this before or after the energy subtraction/addition?
    MSC.Calc_MS(en, thickness*density);

    counter++;
  }
  if (Talk) {
    std::cout << " - energy after " << counter << ": " << en << std::endl;
    std::cout << "Total path length (g/cm2): " << TotalPathLength << std::endl;
    std::cout << "Total length (cm): " << TotalLength << std::endl;
  }

  // Updated momentum^2
  double p_2_up = en*en-BetheBloch_Utils::Mm*BetheBloch_Utils::Mm;
  double p_up;
  if (p_2_up > 0) p_up = sqrt(p_2_up);
  else {
    std::cerr << "[TMS_Kalman.cpp] negative momentum squared, setting momentum to 1 MeV" << std::endl;
    //p_up = 1;
    p_up = sqrt(en*en);
  }

  // Update the state's q/p
  if (Talk) std::cout << "setting current p = " << p_up << std::endl;
  CurrentState.qp = 1./p_up;


  double p_var = (2*en/p_up)*(2*en/p_up) * total_en_var;
  double qp_var = 1./(p_up*p_up*p_up*p_up) * p_var;


  // Set pointers to the noise matrix
  TMatrixD &NoiseMatrix = Node.NoiseMatrix;
  TMatrixD &CovarianceMatrix = Node.CovarianceMatrix;
  TMatrixD &UpdatedCovarianceMatrix = Node.UpdatedCovarianceMatrix;

  // 'measurement' matrices
  TMatrixD H   = TMatrixD(KALMAN_DIM,KALMAN_DIM);
  TMatrixD H_T = TMatrixD(KALMAN_DIM,KALMAN_DIM);
  H  .Zero();
  H_T.Zero();
  for (int l=0; l<2; l++) { H(l,l) = 1.0; H_T(l,l) = 1.0; }

  double sigma = MSC.Calc_MS_Sigma();

  if (TotalPathLength >= 1.0) { // If path is 0 we're in the first Node, set initial cov
    Node.FillUpdatedCovarianceMatrix(TotalPathLength, UpdateVec[2], UpdateVec[3], CurrentState.qp, sigma, false); // Fill the matrix for multiple scattering
  } else { // Initialise cov to something 'sane'-ish
    UpdatedCovarianceMatrix.Zero(); // zero this out to be sure
    CovarianceMatrix(0,0) = 200.0; // TODO: enable swapping the x and y for z layers, should rarely happen
    CovarianceMatrix(1,1) = 1.0E3;
    CovarianceMatrix(2,2) = 1.50;
    CovarianceMatrix(3,3) = 2.50;
    CovarianceMatrix(4,4) = 1.0;
    if (Talk) std::cout << "Initialising covariance!" << std::endl;
  }

  Node.FillNoiseMatrix(); // Full the matrix for multiple scattering

  CovarianceMatrix = Transfer*CovarianceMatrix*TransferT;

  CovarianceMatrix += UpdatedCovarianceMatrix;

  TMatrixD GainMatrix = TMatrixD(5,5);

  GainMatrix = CovarianceMatrix + NoiseMatrix;
  for (int l=2; l<KALMAN_DIM; l++) GainMatrix(l,l) = 1.0; // Set diags to 1 for inversion
  GainMatrix = GainMatrix.Invert();

  GainMatrix = CovarianceMatrix*GainMatrix;

  if (Talk) std::cout << "Final cov\n" << std::flush;
  if (Talk) CovarianceMatrix.Print();

  TVectorD FilteredVec = TVectorD(5);
  TVectorD Measurement = TVectorD(5);

  Measurement[0] = CurrentState.x ;//+ NoiseVec[0];
  Measurement[1] = CurrentState.y ;//+ NoiseVec[1];
  Measurement[2] = UpdateVec[2];//0.0;//NoiseVec[2];//CurrentState.dxdz;
  Measurement[3] = UpdateVec[3];//0.0;//NoiseVec[3];//CurrentState.dydz;
  Measurement[4] = 0.0; //CurrentState.qp;

  if (Talk) std::cout << "Gain" << std::flush;
  if (Talk) GainMatrix.Print();

  FilteredVec = UpdateVec + GainMatrix*( Measurement - UpdateVec );

  CurrentState.x    = FilteredVec[0];
  CurrentState.y    = FilteredVec[1];
  CurrentState.dxdz = FilteredVec[2];
  CurrentState.dydz = FilteredVec[3];


  // Calculate chi^2
  Node.rVec[0] = (Measurement[0] - UpdateVec[0]);
  Node.rVec[1] = (Measurement[1] - UpdateVec[1]);

  // Probably a much nicer way to make (sub)matrix from a bigger one, but YOLO
  Node.RMatrix(0,0) = (Node.NoiseMatrix(0,0) - UpdatedCovarianceMatrix(0,0));
  Node.RMatrix(1,0) = (Node.NoiseMatrix(1,0) - UpdatedCovarianceMatrix(1,0));
  Node.RMatrix(0,1) = (Node.NoiseMatrix(0,1) - UpdatedCovarianceMatrix(0,1));
  Node.RMatrix(1,1) = (Node.NoiseMatrix(1,1) - UpdatedCovarianceMatrix(1,1));
  Node.RMatrix.Invert(); // Matrix has to be inverted

  Node.chi2 = Node.rVec*(Node.RMatrix*Node.rVec); // Calc chi^2



  //CovarianceMatrix.Print();
  //GainMatrix.Print();
  //CurrentState.Print();
  if ( (CurrentState.x < TMS_Const::TMS_Start[0]) || (CurrentState.x > TMS_Const::TMS_End[0]) ) // point outside x region
  {
    std::cerr << "[TMS_Kalman.cpp] x value outside TMS: " << CurrentState.y << "\tTMS: [" << TMS_Const::TMS_Start[0] << ", "<< TMS_Const::TMS_End[0] << "]" << std::endl;
    if (Talk)
    {
      Node.PrintTrueReco();
      PreviousState.Print();
      CurrentState.Print();
    }
  }
  if ( (CurrentState.y < TMS_Const::TMS_Start[1]) || (CurrentState.y > TMS_Const::TMS_End[1]) ) // point outside y region
  {
    std::cerr << "[TMS_Kalman.cpp] y value outside TMS: " << CurrentState.y << "\tTMS: [" << TMS_Const::TMS_Start[1] << ", "<< TMS_Const::TMS_End[1] << "]" << std::endl;
    if (Talk)
    {
      Node.PrintTrueReco();
      PreviousState.Print();
      CurrentState.Print();
    }
  }

  Node.SetRecoXY(CurrentState);
  if (Talk) Node.PrintTrueReco();
  if (Talk) CurrentState.Print();
}

// Use the NoiseMatrix from the Kalman state to throw a vector of random noise // Liam: ??????????????
// to add to the state measurement
TVectorD TMS_Kalman::GetNoiseVector(TMS_KalmanNode Node) {
  // Clarence my man... what is this function even for?
  TVectorD rand_vec;
  rand_vec.ResizeTo(KALMAN_DIM);
  for(int i = 0; i < KALMAN_DIM; i++)
    rand_vec[i] = RNG.Gaus();


  TMatrixD &cov = Node.NoiseMatrix;
  TVectorD toy;
  toy.ResizeTo(5);
  toy.Zero();

  toy = cov*rand_vec;
  rand_vec.Print();

  return toy;
}


void TMS_Kalman::SetStartDirection(double ax, double ay)
{
  // Defining a vector by two slopes is a bit odd;
  // it's not possible to define vectors in the x/y plane for example.
  // Consider the dector (dx/dz, dy/dz, dz/dz), the z component is 1 by construction,
  // and thus we normalise this
  double mag  = sqrt(ax*ax + ay*ay + 1);
  double mag2 = ax*ax + ay*ay + 1; // square of mag

  StartDirection[0]=ax/mag;
  StartDirection[1]=ay/mag;
  StartDirection[2]=sqrt(1 - ax*ax/mag2 - ay*ay/mag2);
}

void TMS_Kalman::SetEndDirection(double ax, double ay)
{
  // Defining a vector by two slopes is a bit odd;
  // it's not possible to define vectors in the x/y plane for example.
  // Consider the dector (dx/dz, dy/dz, dz/dz), the z component is 1 by construction,
  // and thus we normalise this
  double mag  = sqrt(ax*ax + ay*ay + 1);
  double mag2 = ax*ax + ay*ay + 1; // square of mag

  StartDirection[0]=ax/mag;
  StartDirection[1]=ay/mag;
  StartDirection[2]=sqrt(1 - ax*ax/mag2 - ay*ay/mag2);
}
