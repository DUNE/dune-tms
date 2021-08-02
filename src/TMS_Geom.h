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

#define __GEOM_LARGE_STEP__ 1E10
#define __GEOM_SMALL_STEP__ 1E-5
#define __GEOM_TINY_STEP__ 1E-10

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

      // To step through the volumes, can't use IsSameLocation because the volume between LAr and TMS is same as the volume after/next to TMS and in front/beside LAr.
      // Instead check the step size and current Point, making sure it doesn't pass point2
      double target_dist = (point2-point1).Mag();
      double dist = 0;

      double step = geom->GetStep();

      // Count up the total length for debugging
      double total = 0;
      // Walk through until we're in the same volume as our final point
      while (step < __GEOM_LARGE_STEP__ && target_dist-dist > 0) {
        // Get the material of the current point
        TGeoMaterial *mat = geom->GetCurrentNode()->GetMedium()->GetMaterial();
        // Get the position
        const double *pt = geom->GetCurrentPoint();
        // How far have we gone along the current direction
        TVector3 pt_vec(pt[0], pt[1], pt[2]);
        // How far is this boundary from the starting point
        dist = (pt_vec-point1).Mag();

        // Step into the next volume
        geom->FindNextBoundaryAndStep();
        // Go down to the deepest node
        geom->FindNode();
        // How big was the step in this material
        step = geom->GetStep();
        if (step < __GEOM_TINY_STEP__) {
          geom->SetStep(__GEOM_SMALL_STEP__);
          // Step into the next volume
          geom->Step();
          // Go down to the deepest node
          geom->FindNode();
          // How big was the step in this material
          step = geom->GetStep();
        }

        // This step might take us beyond our target, so modify the step size to be from boundary up until point2
        if (dist+step > target_dist) {
          step = (point2-pt_vec).Mag();
        }

        // Push back the information
        std::pair<TGeoMaterial*, double> temp(mat, step);
        Materials.push_back(temp);
        total += step;
      }

      /*
         if (total > 7*1000*10) {
         std::cerr << "Very long distance between points: " << total << std::endl;
         Materials.clear();
         return Materials;
         }
         */

      /*
      // Then finally add in the last material from the boundary to the current point
      // Get the last point
      geom->FindNode(point2.X(), point2.Y(), point2.Z());
      TVector3 temp(curpt[0], curpt[1], curpt[2]);
      // The distance between the 
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
      */

      return Materials;
    }

    std::vector<std::pair<TGeoMaterial*, TVector3> > GetMaterialsPos(const TVector3 &point1, const TVector3 &point2) {
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
      std::vector<std::pair<TGeoMaterial*,TVector3> > Materials;

      // To step through the volumes, can't use IsSameLocation because the volume between LAr and TMS is same as the volume after/next to TMS and in front/beside LAr.
      // Instead check the step size and current Point, making sure it doesn't pass point2
      double target_dist = (point2-point1).Mag();
      double dist = 0;

      double step = geom->GetStep();
      // Walk through until we're in the same volume as our final point
      //while (!geom->IsSameLocation(point2.X(), point2.Y(), point2.Z()) && step < __GEOM_LARGE_STEP__) {
      while (step < __GEOM_LARGE_STEP__ && target_dist-dist > 0) {
        // Get the material of the current point
        TGeoMaterial *mat = geom->GetCurrentNode()->GetMedium()->GetMaterial();
        // Get the position
        const double *pt = geom->GetCurrentPoint();
        // How far have we gone along the current direction
        TVector3 pt_vec(pt[0], pt[1], pt[2]);
        dist = (pt_vec-point1).Mag();

        // Step into the next volume
        geom->FindNextBoundaryAndStep();
        // Go down to the deepest node
        geom->FindNode();
        // ROOT uses a large step in ShootRay function
        // Check IsEntering?
        // If IsEntering, set step to very small until not entering anymore
        step = geom->GetStep();
        if (step < __GEOM_TINY_STEP__) {
          geom->SetStep(__GEOM_SMALL_STEP__);
          // Step into the next volume
          geom->Step();
          // Go down to the deepest node
          geom->FindNode();
          // How big was the step in this material
          step = geom->GetStep();
        }
        // Check that this step still puts us inside the target
        if (dist+step > target_dist) break;

        // Push back the information
        std::pair<TGeoMaterial*, TVector3> temp(mat, pt_vec);
        Materials.push_back(temp);
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
        //material.first->Print();

        /*
        std::cout << "Material " << counter << " = " << material.first->GetName() << std::endl;
        std::cout << "  density: " << density << std::endl;
        std::cout << "  thickness: " << thickness << std::endl;
        std::cout << "  thickness*density = " << density*thickness << std::endl;
        */

        // Skip if density or thickness is small
        if (density*thickness < 0.1) {
          //std::cout << "  Skipping material, to little path length to bother" << std::endl;
          continue;
        }

        TotalPathLength += density*thickness;
        TotalLength += thickness;

        counter++;
      }

      return TotalPathLength;
    };

    std::vector<std::pair<std::string, TVector3> > GetNodes(const TVector3 &point1, const TVector3 &point2) {

      // First cd the navigator to the starting point
      geom->FindNode(point1.X(), point1.Y(), point1.Z());

      // Go between the two points
      // Set the current point to be the current point
      geom->SetCurrentPoint(point1.X(), point1.Y(), point1.Z());
      // Set the direction
      geom->SetCurrentDirection((point2-point1).Unit().X(), 
          (point2-point1).Unit().Y(), 
          (point2-point1).Unit().Z());

      std::vector<std::pair<std::string, TVector3> > Nodes;

      // To step through the volumes, can't use IsSameLocation because the volume between LAr and TMS is same as the volume after/next to TMS and in front/beside LAr.
      // Instead check the step size and current Point, making sure it doesn't pass point2
      double target_dist = (point2-point1).Mag();
      double dist = 0;

      double step = geom->GetStep();

      // Walk through until we're in the same volume as our final point
      while (step < __GEOM_LARGE_STEP__ && target_dist-dist > 0) {
        // Get the material of the current point
        std::string nodename = std::string(geom->GetCurrentNode()->GetName());
        const double *pt = geom->GetCurrentPoint();
        // How far have we gone along the current direction
        TVector3 pt_vec(pt[0], pt[1], pt[2]);
        // How far is this boundary from the starting point
        dist = (pt_vec-point1).Mag();

        // Step into the next volume
        geom->FindNextBoundaryAndStep();
        // Go down to the deepest node
        geom->FindNode();
        // How big was the step in this material
        step = geom->GetStep();
        if (step < __GEOM_TINY_STEP__) {
          geom->SetStep(__GEOM_SMALL_STEP__);
          // Step into the next volume
          geom->Step();
          // Go down to the deepest node
          geom->FindNode();
          // How big was the step in this material
          step = geom->GetStep();
        }

        // Push back the information
        std::pair<std::string, TVector3> temp(nodename, pt_vec);
        Nodes.push_back(temp);
      }
      return Nodes;
    }

    std::vector<std::pair<int*, TVector3> > GetUniquePlaneBarIdent(const TVector3 &point1, const TVector3 &point2) {

      // First cd the navigator to the starting point
      geom->FindNode(point1.X(), point1.Y(), point1.Z());

      // Go between the two points
      // Set the current point to be the current point
      geom->SetCurrentPoint(point1.X(), point1.Y(), point1.Z());
      // Set the direction
      geom->SetCurrentDirection((point2-point1).Unit().X(), 
          (point2-point1).Unit().Y(), 
          (point2-point1).Unit().Z());

      std::vector<std::pair<int*, TVector3> > Nodes;

      // To step through the volumes, can't use IsSameLocation because the volume between LAr and TMS is same as the volume after/next to TMS and in front/beside LAr.
      // Instead check the step size and current Point, making sure it doesn't pass point2
      double target_dist = (point2-point1).Mag();
      double dist = 0;

      double step = geom->GetStep();

      // Walk through until we're in the same volume as our final point
      while (step < __GEOM_LARGE_STEP__ && target_dist-dist > 0) {
        // Plane, bar, global
        int *Plane = new int[3];
        for (int i = 0; i < 3; ++i) Plane[i] = -999;

        std::string NodeName = std::string(geom->GetCurrentNode()->GetName());
        const double *pt = geom->GetCurrentPoint();
        // How far have we gone along the current direction
        TVector3 pt_vec(pt[0], pt[1], pt[2]);
        // How far is this boundary from the starting point
        dist = (pt_vec-point1).Mag();

        // cd up in the geometry to find the right name
        while (NodeName.find(TMS_Const::TMS_DetEnclosure) == std::string::npos && 
            NodeName.find(TMS_Const::TMS_TopLayerName) == std::string::npos) {
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

        // cd down back into the deepest node
        geom->FindNode();
        // Step into the next volume
        geom->FindNextBoundaryAndStep();
        // How big was the step in this material
        step = geom->GetStep();
        if (step < __GEOM_TINY_STEP__) {
          geom->SetStep(__GEOM_SMALL_STEP__);
          // Step into the next volume
          geom->Step();
          // Go down to the deepest node
          geom->FindNode();
          // How big was the step in this material
          step = geom->GetStep();
        }

        // Don't push back if we're missing info; likely means the volume wasn't a scintillator bar
        if (Plane[0] == -999 || Plane[1] == -999 || Plane[2] == -999) continue;

        // Push back the information
        std::pair<int*, TVector3> temp(Plane, pt_vec);
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
