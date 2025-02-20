#include <iostream>

#include "BetheBloch.h"

#include "TLegend.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"

int main() {
  const int npoints = 10000;

  std::vector<Material::MaterialType> Materials;
  Materials.push_back(Material::kIron);
  Materials.push_back(Material::kPolyStyrene);
  Materials.push_back(Material::kGraphite);
  Materials.push_back(Material::kGArgon);
  Materials.push_back(Material::kLArgon);
  Materials.push_back(Material::kWater);
  int nMaterials = Materials.size();

  std::vector<TGraph*> graphs;
  for (int i = 0; i < nMaterials; ++i) graphs.push_back(new TGraph(npoints));

  // Just need one calculator which we change the material of
  BetheBloch_Calculator generic(Material::kIron);

  double minimum = BetheBloch_Utils::Mm*1.2;
  for (int i = 0; i < npoints; ++i) {
    double energy = minimum+(i*0.5);
    for (int j = 0; j < nMaterials; ++j) {
      generic.fMaterial = Materials[j];
      double dedx = generic.Calc_dEdx(energy);
      graphs[j]->SetPoint(i, energy, dedx);
    }
  }

  TCanvas *canv = new TCanvas("canv", "canv", 1024, 1024);
  for (int i = 0; i < nMaterials; ++i) {
    generic.fMaterial = Materials[i];
    graphs[i]->SetTitle(generic.fMaterial.MaterialName().c_str());
    graphs[i]->SetLineColor(i+1);
    if (i == 0) {
      graphs[i]->Draw();
      graphs[i]->GetXaxis()->SetTitle("E_{#mu} (MeV)");
      graphs[i]->GetYaxis()->SetTitle("<dE/dX> (MeV/gr/cm^{2})");
    }
    else graphs[i]->Draw("same");
  }
  TLegend *leg = new TLegend(0.2, 0.6, 0.9, 0.9);
  for (int i = 0; i < nMaterials; ++i) {
    leg->AddEntry(graphs[i], graphs[i]->GetTitle(), "l");
  }
  leg->SetNColumns(2);
  leg->SetLineStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw("same");
  canv->Print("dedx.pdf");

}
