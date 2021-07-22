#ifndef _TMS_GEOM_H_SEEN_
#define _TMS_GEOM_H_SEEN_

#include <iostream>
#include <string>
#include <utility>

#include "TGeoManager.h"

#include "TVector3.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "TMS_Manager.h"
#include "TMS_Constants.h"

// Define the TMS geometry singleton
class TMS_Geom {
  public:

    // Get the instance
    static TMS_Geom& GetInstance() {
      static TMS_Geom Instance;
      return Instance;
    }

    // Get the geometry
    TGeoManager* GetGeometry() {
      // Shout if the geometry isn't set
      if (this->geom == NULL) {
        std::cerr << "Geometry not set!" << std::endl;
        std::cerr << "Call TMS_Geom.SetGeometry(TGeoManager *) first!" << std::endl;
        throw;
      }
      return this->geom;
    }

    const std::string & GetFileName() { 
      return FileName;
    }

    // Set the geometry
    void SetGeometry(TGeoManager *geometry) {
      geom = geometry;
      std::cout << "Global geometry set to " << geometry->GetName() << std::endl;
      geom->LockGeometry();
    }

    void SetFileName(std::string filename) {
      FileName = filename;
    }

    // Largely modlled on TGeoChecker::ShootRay in ROOT (except there's a stopping point, not an infinitely long ray)
    std::vector<std::pair<TGeoMaterial*, double> > GetMaterials(const TVector3 &point1, const TVector3 &point2) {

      // First cd the navigator to the starting point
      geom->FindNode(point1.X(), point1.Y(), point1.Z());

      // Go between the two points
      // Set the current point to be the current point
      geom->SetCurrentPoint(point1.X(), point1.Y(), point1.Z());
      // Set the direction
      geom->SetCurrentDirection((point2-point1).Unit().X(), 
          (point2-point1).Unit().Y(), 
          (point2-point1).Unit().Z());

      // The returned vector of materials
      // Also want how much of the material was passed through
      std::vector<std::pair<TGeoMaterial*,double> > Materials;

      // Count up the total length for debugging
      double total = 0;
      // Walk through until we're in the same volume as our final point
      while (!geom->IsSameLocation(point2.X(), point2.Y(), point2.Z())) {
        // Get the material of the current point
        TGeoMaterial *mat = geom->GetCurrentNode()->GetMedium()->GetMaterial();
        // Step into the next volume
        geom->FindNextBoundaryAndStep();
        // How big was the step
        double snext = geom->GetStep();
        // Push back the information
        std::pair<TGeoMaterial*, double> temp(mat, snext);
        Materials.push_back(temp);
        total += snext;
        // Check the total
        // Detector is roughly 7 meters: 10 times this and something has gone wrong!
        if (total > 7*1000*10) break;
      }

      if (total > 7*1000*10) {
        std::cerr << "Very long distance between points: " << total << std::endl;
        Materials.clear();
        return Materials;
      }

      // Then finally add in the last material too
      // Get the last point
      const Double_t *curpt = geom->GetCurrentPoint();
      TVector3 temp(curpt[0], curpt[1], curpt[2]);
      double extra = (temp-point2).Mag();
      // Change the point to get the material
      geom->SetCurrentPoint(point2.X(), point2.Y(), point2.Z());
      // Update material
      TGeoMaterial *mat = geom->GetCurrentNode()->GetMedium()->GetMaterial();
      std::pair<TGeoMaterial*, double> mypair(mat, extra);
      Materials.push_back(mypair);
      
      if (fabs(total+extra - (point2-point1).Mag()) > 1E-3) {
        std::cout << "Total: " << total << std::endl;
        std::cout << "extra: " << extra << std::endl;
        std::cout << "total+extra: " << total+extra << std::endl;
        std::cout << "Intended: " << (point2-point1).Mag() << std::endl;
        std::cout << "N materials: " << Materials.size() << std::endl;
        throw;
      }

      return Materials;
    }

    // Get the track length between two points by walking through the materials at those points
    double GetTrackLength(const TVector3 &point1, const TVector3 &point2) {
      // First get the collection of materials between point1 and point2
      std::vector<std::pair<TGeoMaterial*, double> > Materials = GetMaterials(point1, point2);

      double TotalPathLength = 0;
      double TotalLength = 0;
      int counter = 0;
      for (auto material : Materials) {

        // Read these directly from a TGeoManager
        double density = material.first->GetDensity()/(CLHEP::g/CLHEP::cm3); // now in g/cm3 (edep-sim geometry is in CLHEP units)
        double thickness = material.second/10.; // in cm (was in mm in geometry)
        material.first->Print();

          std::cout << "Material " << counter << " = " << material.first->GetName() << std::endl;
          std::cout << "  density: " << density << std::endl;
          std::cout << "  thickness: " << thickness << std::endl;
          std::cout << "  thickness*density = " << density*thickness << std::endl;

        // Skip if density or thickness is small
        if (density*thickness < 0.1) {
          std::cout << "  Skipping material, to little path length to bother" << std::endl;
          continue;
        }

        TotalPathLength += density*thickness;
        TotalLength += thickness;

        counter++;
      }

      return TotalPathLength;
    };

    std::vector<std::pair<std::string, const double*> > GetNodes(const TVector3 &point1, const TVector3 &point2) {

      // First cd the navigator to the starting point
      geom->FindNode(point1.X(), point1.Y(), point1.Z());

      // Go between the two points
      // Set the current point to be the current point
      geom->SetCurrentPoint(point1.X(), point1.Y(), point1.Z());
      // Set the direction
      geom->SetCurrentDirection((point2-point1).Unit().X(), 
          (point2-point1).Unit().Y(), 
          (point2-point1).Unit().Z());
      
      std::vector<std::pair<std::string, const double*> > Nodes;

      // Walk through until we're in the same volume as our final point
      while (!geom->IsSameLocation(point2.X(), point2.Y(), point2.Z())) {
        // Get the material of the current point
        std::string nodename = std::string(geom->GetCurrentNode()->GetName());
        const double *pos = geom->GetCurrentPoint();

        // Step into the next volume
        geom->FindNextBoundaryAndStep();

        // Push back the information
        std::pair<std::string, const double*> temp(nodename, pos);
        Nodes.push_back(temp);
      }

      return Nodes;
    }

    std::vector<std::pair<int*, const double*> > GetUniquePlaneBarIdent(const TVector3 &point1, const TVector3 &point2) {

      // First cd the navigator to the starting point
      geom->FindNode(point1.X(), point1.Y(), point1.Z());

      // Go between the two points
      // Set the current point to be the current point
      geom->SetCurrentPoint(point1.X(), point1.Y(), point1.Z());
      // Set the direction
      geom->SetCurrentDirection((point2-point1).Unit().X(), 
          (point2-point1).Unit().Y(), 
          (point2-point1).Unit().Z());
      
      std::vector<std::pair<int*, const double*> > Nodes;

      // Walk through until we're in the same volume as our final point
      while (!geom->IsSameLocation(point2.X(), point2.Y(), point2.Z())) {
        const double *pos = geom->GetCurrentPoint();
        std::string NodeName = std::string(geom->FindNode(pos[0], pos[1], pos[2])->GetName());
        // Plane, bar, global
        int *Plane = new int[3];

        // cd up in the geometry to find the right name
        while (NodeName.find(TMS_Const::TMS_TopLayerName) == std::string::npos &&
               NodeName.find(TMS_Const::TMS_EDepSim_VolumeName) == std::string::npos) {

          // We've found the plane number
          if (NodeName.find(TMS_Const::TMS_ModuleLayerName) != std::string::npos) {
            Plane[0] = geom->GetCurrentNode()->GetNumber();
          } else if (NodeName.find(TMS_Const::TMS_ScintLayerName) != std::string::npos) {
            Plane[1] = geom->GetCurrentNode()->GetNumber();
          } else if (NodeName.find(TMS_Const::TMS_ModuleName) != std::string::npos) {
            Plane[2] = geom->GetCurrentNode()->GetNumber();
          }

          geom->CdUp();
          NodeName = std::string(geom->GetCurrentNode()->GetName());
        }

        // Step into the next volume
        geom->FindNextBoundaryAndStep();

        // Push back the information
        std::pair<int*, const double*> temp(Plane, pos);
        Nodes.push_back(temp);
      }

      return Nodes;
    }


  private:
    // The empty constructor
    TMS_Geom() {
      FileName = TMS_Manager::GetInstance().GetFileName();
      geom = NULL;
    };

    ~TMS_Geom() {};

    TMS_Geom(TMS_Geom const &) = delete;
    void operator=(TMS_Geom const &) = delete;

    // The actual geometry
    TGeoManager *geom;
    std::string FileName;
};

#endif
