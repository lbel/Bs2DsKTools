#include <iostream>

#include <TStyle.h>
#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TPad.h>

#include "oanglePID.hxx"
#include "readDataTree.hxx"
#include "PIDsel.hh"

// #ifdef __CINT__
// gSystem->Load("oanglePID_cxx.so");
// #endif

// Dsπ / DsK = 0.067 or 0.081

using namespace std;

int PIDsel()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetPalette(1); // "rainbow" color palette
  gStyle->SetNumberContours(256); // smooth color palette
  gStyle->SetTitleOffset( 1.2, "xy");

  TFile *fhisto = new TFile( "templates.root", "read");

  TH2D *hDsK  = dynamic_cast<TH2D*> (fhisto->Get("hDsK" )->Clone());
  TH2D *hDspi = dynamic_cast<TH2D*> (fhisto->Get("hDspi")->Clone());

  oanglePID kPID(hDsK, hDspi);
  
  // make and read chain
  TChain K_hypo("K_hypo");
  K_hypo.Add("../../ntuples/data/FitTuple_BsDs1Pi_K_*.root/MyOffSelTree");
  readDataTree ch_K (&K_hypo);

  // TChain pi_hypo("pi_hypo");
  // pi_hypo.Add("../../ntuples/data/FitTuple_BsDs1Pi_Pi_*.root/MyOffSelTree");
  // readDataTree ch_pi(&pi_hypo);

  // invariant mass
  TH1D hdataBsmK ("hdataBsmK", "B_{s} mass with #it{#pi} mass hypothesis", 75, 4500, 6000);
  // TH1D hdataBsmpi  ("hdataBsmpi",  "B_{s} mass with #it{K} mass hypothesis",   75, 4500, 6000);

  hdataBsmK.SetLineColor(kAzure);
  hdataBsmK.SetXTitle("Mass[MeV]");
  hdataBsmK.SetYTitle("Events");
  hdataBsmK.Sumw2();

  // opening angle
  TH2D hdataDsK ("hdataDsK", "#it{B_{s}} (#it{D_{s}#pi} hypothesis) vs #it{#vec{#beta}} #angle #it{h} in #it{B_{s}} rest frame",
		 75, 4500, 6000, 20, -1, 1);
  // TH2D hdataDspi ("hdataDspi", "#it{B_{s}} (#it{D_{s}K} hypothesis) vs #it{#vec{#beta}} #angle #it{h} in #it{B_{s}} rest frame",
  // 		  75, 4500, 6000, 20, -1, 1);

  hdataDsK.SetYTitle("Cosine of the opening angle[deg]");
  hdataDsK.SetXTitle("Mass[MeV]");
  hdataDsK.Sumw2();

  ch_K.Loop(kPID, hdataBsmK, hdataDsK);

  hdataBsmK.Draw("hist");
  gPad->Print("Bs-mass-data-new-PID-K.png");
  gPad->Clear();
  hdataDsK.Draw("COLZ");
  gPad->Print("Bs-opening-angle-data-new-PID-K.png");

  hdataBsmK.Print("base");
  hdataDsK .Print("base");

  return 0;
}


int test_oanglePID()
{
  gStyle->SetPalette(1); // "rainbow" color palette
  gStyle->SetNumberContours(256); // smooth color palette

  TFile *fhisto = new TFile( "templates.root", "read");

  TH2D *hDsK  = dynamic_cast<TH2D*> (fhisto->Get("hDsK" )->Clone());
  TH2D *hDspi = dynamic_cast<TH2D*> (fhisto->Get("hDspi")->Clone());

  oanglePID test(hDsK, hDspi);

  TH2D *hweight = dynamic_cast<TH2D*> (hDsK->Clone());
  hweight->Reset("ICESM");

  // TH2D binning - 150, 4500, 6000, 25, -1, 1);
  int gbin(0);
  for (double mass(4505.); mass < 6000; mass+=10) {
    for (double cosine(-0.98); cosine < 1; cosine+=0.08) {
      gbin = hweight->FindBin(mass, cosine);
      hweight->SetBinContent(gbin, test.GetoangleDLL(mass, cosine));
      // cout << "bin: " << gbin;
      // cout << ", DLL: " << test.GetoangleDLL(mass, cosine) << endl;
    }
  }

  TCanvas c("c", "template compare", 1500, 400);
  c.Divide(3,1);

  c.cd(1);
  hDsK->Draw("COLZ");

  c.cd(2);
  hDspi->Draw("COLZ");

  c.cd(3);
  hweight->Draw("COLZ");

  c.Print("new-PID-test.png");
  // hweight->Divide(hDsK, hDspi, 0.1, 0.9);
  // hweight->Draw("surf");

  return 0;
}
