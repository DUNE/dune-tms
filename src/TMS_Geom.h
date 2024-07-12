#ifndef _TMS_GEOM_H_SEEN_
#define _TMS_GEOM_H_SEEN_

#include <iostream>
#include <string>
#include <utility>

#include "TGeoManager.h"

#include "TVector3.h"
#include "TGeoBBox.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "TMS_Manager.h"
#include "TMS_Constants.h"

#define __GEOM_LARGE_STEP__ 1E10
#define __GEOM_SMALL_STEP__ 1E-4
#define __GEOM_TINY_STEP__ 1E-5

// Define the TMS geometry singleton
class TMS_Geom {
  public:

    // Get the instance
    static TMS_Geom& GetInstance() {
      static TMS_Geom Instance;
      return Instance;
    }
    
    // Positions are fidicual volume
    inline double GetXStartOfLAr() const { return TMS_Const::LAr_Start_Exact[0]; };
    inline double GetYStartOfLAr() const { return TMS_Const::LAr_Start_Exact[1]; };
    inline double GetZStartOfLAr() const { return TMS_Const::LAr_Start_Exact[2]; };
    inline TVector3 GetStartOfLAr() const { return TVector3(GetXStartOfLAr(), GetYStartOfLAr(), GetZStartOfLAr()); };
    inline double GetXEndOfLAr() const { return TMS_Const::LAr_End_Exact[0]; };
    inline double GetYEndOfLAr() const { return TMS_Const::LAr_End_Exact[1]; };
    inline double GetZEndOfLAr() const { return TMS_Const::LAr_End_Exact[2]; };
    inline TVector3 GetEndOfLAr() const { return TVector3(GetXEndOfLAr(), GetYEndOfLAr(), GetZEndOfLAr()); };
    inline double GetXStartOfTMS() const { return TMS_Const::TMS_Start_Bars_Only[0]; };
    inline double GetYStartOfTMS() const { return TMS_Const::TMS_Start_Bars_Only[1]; };
    inline double GetZStartOfTMS() const { return TMS_Const::TMS_Start_Bars_Only[2]; };
    inline TVector3 GetStartOfTMS() const { return TVector3(GetXStartOfTMS(), GetYStartOfTMS(), GetZStartOfTMS()); };
    inline double GetXEndOfTMS() const { return TMS_Const::TMS_End_Bars_Only[0]; };
    inline double GetYEndOfTMS() const { return TMS_Const::TMS_End_Bars_Only[1]; };
    inline double GetZEndOfTMS() const { return TMS_Const::TMS_End_Bars_Only[2]; };
    inline TVector3 GetEndOfTMS() const { return TVector3(GetXEndOfTMS(), GetYEndOfTMS(), GetZEndOfTMS()); };
    inline double GetXStartOfTMSMass() const { return TMS_Const::TMS_Start_Exact[0]; };
    inline double GetYStartOfTMSMass() const { return TMS_Const::TMS_Start_Exact[1]; };
    inline double GetZStartOfTMSMass() const { return TMS_Const::TMS_Start_Exact[2]; };
    inline TVector3 GetStartOfTMSMass() const { return TVector3(GetXStartOfTMS(), GetYStartOfTMSMass(), GetZStartOfTMSMass()); };
    inline double GetXEndOfTMSMass() const { return TMS_Const::TMS_End_Exact[0]; };
    inline double GetYEndOfTMSMass() const { return TMS_Const::TMS_End_Exact[1]; };
    inline double GetZEndOfTMSMass() const { return TMS_Const::TMS_End_Exact[2]; };
    inline TVector3 GetEndOfTMSMass() const { return TVector3(GetXEndOfTMS(), GetYEndOfTMSMass(), GetZEndOfTMSMass()); };
    inline double GetZEndOfTMSThin() const { return TMS_Const::TMS_Thick_Start; };
    inline TVector3 GetEndOfTMSThin() const { return TVector3(GetXEndOfTMS(), GetYEndOfTMS(), GetZEndOfTMSThin()); };
    inline double GetZEndOfTMSFirstTwoModules() const { return GetZStartOfTMS() + 110; }; // module 2 - module 0 = 11cm
    inline TVector3 GetEndOfTMSFirstTwoModules() const { return TVector3(GetXEndOfTMS(), GetYEndOfTMS(), GetZEndOfTMSFirstTwoModules()); };

    inline TVector3 GetStartOfTMSFiducial() const { return TVector3(TMS_Manager::GetInstance().Get_FIDUCIAL_TMS_START_X(), 
        TMS_Manager::GetInstance().Get_FIDUCIAL_TMS_START_Y(), TMS_Manager::GetInstance().Get_FIDUCIAL_TMS_START_Z()); };
    inline TVector3 GetEndOfTMSFiducial() const { return TVector3(TMS_Manager::GetInstance().Get_FIDUCIAL_TMS_END_X(), 
        TMS_Manager::GetInstance().Get_FIDUCIAL_TMS_END_Y(), TMS_Manager::GetInstance().Get_FIDUCIAL_TMS_END_Z()); };
    inline TVector3 GetStartOfLArFiducial() const { return TVector3(TMS_Manager::GetInstance().Get_FIDUCIAL_LAR_START_X(), 
        TMS_Manager::GetInstance().Get_FIDUCIAL_LAR_START_Y(), TMS_Manager::GetInstance().Get_FIDUCIAL_LAR_START_Z()); };
    inline TVector3 GetEndOfLArFiducial() const { return TVector3(TMS_Manager::GetInstance().Get_FIDUCIAL_LAR_END_X(), 
        TMS_Manager::GetInstance().Get_FIDUCIAL_LAR_END_Y(), TMS_Manager::GetInstance().Get_FIDUCIAL_LAR_END_Z()); };
    
    bool IsInsideBox(TVector3 position, TVector3 start, TVector3 end) const {
      if (position.X() < start.X()) return false;
      if (position.Y() < start.Y()) return false;
      if (position.Z() < start.Z()) return false;
      if (position.X() > end.X()) return false;
      if (position.Y() > end.Y()) return false;
      if (position.Z() > end.Z()) return false;
      return true;
    };
    bool IsInsideLAr(TVector3 position) const { return IsInsideBox(position, GetStartOfLArFiducial(), GetEndOfLArFiducial()); };
    static bool StaticIsInsideLAr(TVector3 position) { return TMS_Geom::GetInstance().IsInsideLAr(position); };
    bool IsInsideTMS(TVector3 position) const { return IsInsideBox(position, GetStartOfTMSFiducial(), GetEndOfTMSFiducial()); };
    static bool StaticIsInsideTMS(TVector3 position) { return TMS_Geom::GetInstance().IsInsideTMS(position); };
    bool IsInsideTMSThin(TVector3 position) const { return IsInsideBox(position, GetStartOfTMSFiducial(), GetEndOfTMSThin()); };
    static bool StaticIsInsideTMSThin(TVector3 position) { return TMS_Geom::GetInstance().IsInsideTMSThin(position); };
    bool IsInsideTMSFirstTwoModules(TVector3 position) const { return IsInsideBox(position, GetStartOfTMSFiducial(), GetEndOfTMSFirstTwoModules()); };
    static bool StaticIsInsideTMSFirstTwoModules(TVector3 position) { return TMS_Geom::GetInstance().IsInsideTMSFirstTwoModules(position); };
    bool IsInsideTMSMass(TVector3 position) const { return IsInsideBox(position, GetStartOfTMSMass(), GetEndOfTMSMass()); };
    static bool StaticIsInsideTMSMass(TVector3 position) { return TMS_Geom::GetInstance().IsInsideTMSMass(position); };
    
    bool IsInsideReasonableSize(TVector3 position) const { return IsInsideBox(position, TVector3(-10000, -10000, 3000), TVector3(10000, 10000, 20000)); };
    static bool StaticIsInsideReasonableSize(TVector3 position) { return TMS_Geom::GetInstance().IsInsideReasonableSize(position); };
    
    
    std::string GetNameOfDetector(const TVector3 &point) {
      std::string out = "";
      if (out == "" && IsInsideLAr(point)) out += "LAr";
      if (out == "" && IsInsideTMS(point)) out += "TMS";
      if (out == "" && IsInsideTMSMass(point)) out += "TMS Mass";
      if (out == "" && StaticIsInsideReasonableSize(point)) out += "Reasonable Size";
      if (out == "") out += "Other";
      out += " (";
      out += std::to_string(point.X()) + ", ";
      out += std::to_string(point.Y()) + ", ";
      out += std::to_string(point.Z());
      out += ")";
      return out;
    };
    

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
      geom->LockGeometry();
      
      // There's an overall scale factor depending on if the gmdl was loaded with cm or mm. 
      // So here we're trying to automatically figure out the scale factor
      // If the geometry changes too much, then this would not work and so we'd want to give up.
      TGeoBBox *box = dynamic_cast<TGeoBBox*>(geom->GetTopVolume()->GetShape());
      double dx = 2*box->GetDX();
      if (dx == 600000) ScaleFactor = 1;
      else if (dx == 60000) ScaleFactor = 10;
      else {
       std::cout << "DX: " << box->GetDX() << std::endl;
       std::cerr << "Fatal: Unable to guess geometry's scale factor based on Shape for geometry " << geometry->GetName() << std::endl;
       throw;
      }
      // set ScaleFactor temporarilty to 10
      //ScaleFactor = 10;
      std::cout << "Global geometry set to " << geometry->GetName() << std::endl;
      std::cout << "Geometry scale factor: " << ScaleFactor;
      std::cout << ". Factor is 1 if 1 unit = 1mm (aka edep sim), 10 if 1 unit = 1cm (aka larsoft)." << std::endl;
    }

    void SetFileName(std::string filename) {
      FileName = filename;
    }
    
    inline void Scale(double &x, double &y, double &z) {
      x *= ScaleFactor;
      y *= ScaleFactor;
      z *= ScaleFactor;
    }
    
    inline double Scale(double val) {
      return val * ScaleFactor;
    }
    
    TVector3 Scale(TVector3 point) {
      double x = point.X();
      double y = point.Y();
      double z = point.Z();
      Scale(x, y, z);
      TVector3 out(x, y, z);
      return out;
    }
    
    inline void Unscale(double &x, double &y, double &z) {
      x /= ScaleFactor;
      y /= ScaleFactor;
      z /= ScaleFactor;
    }
    
    inline double Unscale(double val) {
      return val / ScaleFactor;
    }
    
    TVector3 Unscale(TVector3 point) {
      double x = point.X();
      double y = point.Y();
      double z = point.Z();
      Unscale(x, y, z);
      TVector3 out(x, y, z);
      return out;
    }
    
    TGeoNode* FindNode(double x, double y, double z) {
      Unscale(x, y, z);
      return geom->FindNode(x, y, z);
    }

    // Largely modlled on TGeoChecker::ShootRay in ROOT (except there's a stopping point, not an infinitely long ray)
    std::vector<std::pair<TGeoMaterial*, double> > GetMaterials(const TVector3 &point1_temp, const TVector3 &point2_temp) {
      // Make vectors have geometry scale
      TVector3 point1 = Unscale(point1_temp);
      TVector3 point2 = Unscale(point2_temp);

      // The returned vector of materials
      // Also want how much of the material was passed through
      std::vector<std::pair<TGeoMaterial*,double> > Materials;     

      if ((point1 - point2).Mag() == 0) return Materials;
      // First cd the navigator to the starting point
      geom->FindNode(point1.X(), point1.Y(), point1.Z());

      // Go between the two points
      // Set the current point to be the current point
      geom->SetCurrentPoint(point1.X(), point1.Y(), point1.Z());
      // Set the direction
      geom->SetCurrentDirection((point2-point1).Unit().X(), 
          (point2-point1).Unit().Y(), 
          (point2-point1).Unit().Z());

      // To step through the volumes, can't use IsSameLocation because the volume between LAr and TMS is same as the volume after/next to TMS and in front/beside LAr.
      // Instead check the step size and current Point, making sure it doesn't pass point2
      double target_dist = (point2-point1).Mag();
      double dist = 0;

      double step = geom->GetStep();

      // Count up the total length for debugging
      double total = 0;
      // Walk through until we're in the same volume as our final point
      while (step < Unscale(__GEOM_LARGE_STEP__) && target_dist-dist > 0) {
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
        //std::cout<<"Current total: "<<total<<", det "<<GetNameOfDetector(pt_vec)<<", step="<<step<<std::endl;
        if (step < Unscale(__GEOM_TINY_STEP__)) {
          geom->SetStep(Unscale(__GEOM_SMALL_STEP__));
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
        std::pair<TGeoMaterial*, double> temp(mat, Scale(step));
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

    std::vector<std::pair<TGeoMaterial*, TVector3> > GetMaterialsPos(const TVector3 &point1_temp, const TVector3 &point2_temp) {
      // Make vectors have geometry scale
      TVector3 point1 = Unscale(point1_temp);
      TVector3 point2 = Unscale(point2_temp);
      
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
      //while (!geom->IsSameLocation(point2.X(), point2.Y(), point2.Z()) && step < Unscale(__GEOM_LARGE_STEP__)) {
      while (step < Unscale(__GEOM_LARGE_STEP__) && target_dist-dist > 0) {
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
        if (step < Unscale(__GEOM_TINY_STEP__)) {
          geom->SetStep(Unscale(__GEOM_SMALL_STEP__));
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
        std::pair<TGeoMaterial*, TVector3> temp(mat, Scale(pt_vec));
        Materials.push_back(temp);
      }
      return Materials;
    }

    // Get the track length between two points by walking through the materials at those points
    double GetTrackLength(const TVector3 &point1_temp, const TVector3 &point2_temp) {
      // Make vectors have geometry scale
      TVector3 point1 = Unscale(point1_temp);
      TVector3 point2 = Unscale(point2_temp);
      
      // Root's geometry has a bug that if you call a point outside the world (or maybe inside sand?)
      // then from then on we get no materials. Not sure what's going on, but we can avoid it
      // by completely avoiding locations outside "reasonable size"
      // Checked and GetMaterials returns a vector of size zero once the bug is triggered
      if (!IsInsideReasonableSize(point1) || !IsInsideReasonableSize(point2)) return 0;

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

        TotalPathLength += density*thickness;
        TotalLength += thickness;

        counter++;
      }
      TotalPathLength = Scale(TotalPathLength);
      return TotalPathLength;
    };
    
    // Walks through all the pairs of nodes and returns the total track length
    double GetTrackLength(std::vector<TVector3> nodes, bool ignore_y = false) {
      //std::cout<<"Getting track length for "<<nodes.size()<<" nodes"<<std::endl;
      double out = 0;
      if (nodes.size() != 19) return out;
      // for unsigned ints, 0-1 = MAX_UNSIGNED_INT, so need to verify that size > 0
      if (nodes.size() > 1) { 
        // Loop over all pairs of vectors
        for (size_t i = 0; i < nodes.size() - 1; i++) {
          auto p1 = nodes.at(i);
          auto p2 = nodes.at(i+1);
          //std::cout<<"p1 is in "<<GetNameOfDetector(p1)<<" volume, ";
          //std::cout<<"p2 is in "<<GetNameOfDetector(p2)<<" volume"<<std::endl;
          if (ignore_y) {
            TVector3 p1_without_y(p1);
            TVector3 p2_without_y(p2);
            p1_without_y.SetY(-200);
            p2_without_y.SetY(-200);
            p1 = p1_without_y;
            p2 = p2_without_y;
          }
          out += GetTrackLength(p1, p2);
        }
        double distance = (nodes.front() - nodes.back()).Mag();
        if (distance > 0.1 && out < 0.01) {
          std::cout<<"Found track length < 0.01 with distance="<<distance<<", track length="<<out<<std::endl;
          std::vector<std::pair<TGeoMaterial*, double> > Materials = GetMaterials(nodes.front(), nodes.back());
          std::cout<<"Found materials for nodes.front to nodes.back, and I get this many:"<<Materials.size()<<std::endl;
          for (auto& mat : Materials) {
            std::cout<<mat.first<<", "<<mat.second<<std::endl;
          }
        }
      }
      //std::cout<<"Finished track length for "<<nodes.size()<<" nodes and got "<<out<<" g/cm^2"<<std::endl;
      return out;
    }

    std::vector<std::pair<std::string, TVector3> > GetNodes(const TVector3 &point1_temp, const TVector3 &point2_temp) {
      // Make vectors have geometry scale
      TVector3 point1 = Unscale(point1_temp);
      TVector3 point2 = Unscale(point2_temp);

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
      while (step < Unscale(__GEOM_LARGE_STEP__) && target_dist-dist > 0) {
        // Get the material of the current point
        std::string nodename = std::string(geom->GetCurrentNode()->GetName());
        const double *pt = geom->GetCurrentPoint();
        // How far have we gone along the current direction
        TVector3 pt_vec(pt[0], pt[1], pt[2]);
        // How far is this boundary from the starting point
        dist = (pt_vec-point1).Mag();

        std::cout << "new node" << std::endl;
        // Go down to the deepest node
        geom->FindNode();
        std::cout << geom->GetCurrentNode()->GetName() << std::endl;

        // Step into the next volume
        geom->FindNextBoundaryAndStep();
        // Go down to the deepest node
        geom->FindNode();
        std::cout << geom->GetCurrentNode()->GetName() << std::endl;
        // How big was the step in this material
        step = geom->GetStep();
        if (step < Unscale(__GEOM_TINY_STEP__)) {
          geom->SetStep(Unscale(__GEOM_SMALL_STEP__));
          // Step into the next volume
          geom->Step();
          // Go down to the deepest node
          geom->FindNode();
          // How big was the step in this material
          step = geom->GetStep();
        }
        // Try cding around
        std::cout << "***" << nodename << std::endl;
        while (nodename.find(TMS_Const::TMS_DetEnclosure) == std::string::npos) {
          geom->CdUp();
          nodename = std::string(geom->GetCurrentNode()->GetName());
          std::cout << nodename << std::endl;
        }

        // Push back the information
        std::pair<std::string, TVector3> temp(nodename, Scale(pt_vec));
        Nodes.push_back(temp);
      }
      return Nodes;
    }

    std::vector<std::pair<int*, TVector3> > GetUniquePlaneBarIdent(const TVector3 &point1_temp, const TVector3 &point2_temp) {
      // Make vectors have geometry scale
      TVector3 point1 = Unscale(point1_temp);
      TVector3 point2 = Unscale(point2_temp);

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
      while (step < Unscale(__GEOM_LARGE_STEP__) && target_dist-dist > 0) {
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
        if (step < Unscale(__GEOM_TINY_STEP__)) {
          geom->SetStep(Unscale(__GEOM_SMALL_STEP__));
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
        std::pair<int*, TVector3> temp(Plane, Scale(pt_vec));
        Nodes.push_back(temp);
      }

      return Nodes;
    }


      private:
    // The empty constructor
    TMS_Geom() {
      FileName = TMS_Manager::GetInstance().GetFileName();
      geom = NULL;
      ScaleFactor = 1;
    };

    ~TMS_Geom() {};

    TMS_Geom(TMS_Geom const &) = delete;
    void operator=(TMS_Geom const &) = delete;

    // The actual geometry
    TGeoManager *geom;
    std::string FileName;
    double ScaleFactor;
    };

#endif
