#include <string>
#include <vector>
#include <iostream>
#include <cassert>
#include <limits>
#include <time.h>

#include <boost/foreach.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

int main (int argc, char **argv)
{
  // define option flags
  po::options_description options("Options");
  options.add_options()
    ("help,h", "Show this help message")
    ;
  // define positional args
  po::options_description pos_args("Hidden options");
  pos_args.add_options()
    ("inputs", po::value<std::vector<std::string> >(),
     "ROOT files with TTree.  Name of the tree is assumed to be `DecayTree'.")
    ;
  // positional args
  po::positional_options_description pos;
  pos.add("inputs", -1);
  // combined options
  po::options_description all;
  all.add(options).add(pos_args);

  // parse arguments
  po::variables_map vm;
  try {
    // FIXME: they throw different exceptions, handle separately
    store(po::command_line_parser(argc, argv).
	  options(all).positional(pos).run(), vm); // po::required_option
    notify(vm);					   // po::error
  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    return 2;
  }

  /* handle all options */
  std::string prog(argv[0]);
  auto usage = [prog, pos_args, options](){
    std::cout << "Usage: $ " << prog << " [options] " << "<files>\n" << std::endl;
    std::cout << "files - " << pos_args.find_nothrow("inputs", false)->description()
              << std::endl << std::endl;
    std::cout << options << std::endl;
  };
  // do nothing if help message requested
  if (vm.count("help")) {
    usage();
    return 0;
  }

  // enforce mandatory arguments
  if (not vm.count("inputs"))
    {
      std::cout << "ERROR: Input file missing (mandatory).\n" << std::endl;
      usage();
      return 1;
    }

  // input files
  std::vector<std::string> filenames(vm["inputs"].as<std::vector<std::string> >());
  BOOST_FOREACH(std::string filename, filenames) {
    std::cout << "Reading from file: " << filename << std::endl;
    TFile ifile(filename.c_str(), "update");
    TTree *itree = dynamic_cast<TTree*>(ifile.Get("DecayTree"));
    assert(itree);

    std::string newFileName = filename.substr(0, filename.length() - 5) + "_extended.root";
    std::cout << "Writing to file: " << newFileName << std::endl;
    TFile ofile(newFileName.c_str(), "create");
    TTree *otree = itree->CloneTree(0);
    otree->SetAutoFlush(10000);

    unsigned long nentries(itree->GetEntries());
    std::cout << "Entries: " << nentries << std::endl;

    double lab1_gp,
      lab3_gp, lab3_pt, lab3_chi2,
      lab4_gp, lab4_pt, lab4_chi2,
      lab5_gp, lab5_pt, lab5_chi2;

    double
      lab3_CosTheta, lab3_IP_OWNPV, lab3_p, lab3_pe, lab3_px, lab3_py, lab3_pz, lab3_M, lab3_PIDe, lab3_PIDmu, lab3_PIDK, lab3_PIDp,
      lab5_CosTheta, lab5_IP_OWNPV, lab5_p, lab5_pe, lab5_px, lab5_py, lab5_pz, lab5_M, lab5_PIDe, lab5_PIDmu, lab5_PIDK, lab5_PIDp;
    int lab3_ID, lab5_ID;

    double
      lab3_ProbNNe, lab3_ProbNNk, lab3_ProbNNp, lab3_ProbNNpi, lab3_ProbNNmu, lab3_ProbNNghost, lab3_TRACK_CHI2NDOF,
      lab5_ProbNNe, lab5_ProbNNk, lab5_ProbNNp, lab5_ProbNNpi, lab5_ProbNNmu, lab5_ProbNNghost, lab5_TRACK_CHI2NDOF;
    bool lab3_isMuon, lab5_isMuon;

    unsigned int ndswitchvars = 22;
    double* dswitchvars3[] = {&lab3_CosTheta, &lab3_IP_OWNPV, &lab3_p, &lab3_pe, &lab3_px, &lab3_py, &lab3_pz, &lab3_M, &lab3_PIDe, &lab3_PIDmu, &lab3_PIDK, &lab3_PIDp, &lab3_ProbNNe, &lab3_ProbNNk, &lab3_ProbNNp, &lab3_ProbNNpi, &lab3_ProbNNmu, &lab3_ProbNNghost, &lab3_TRACK_CHI2NDOF, &lab3_gp, &lab3_pt, &lab3_chi2};
    double* dswitchvars5[] = {&lab5_CosTheta, &lab5_IP_OWNPV, &lab5_p, &lab5_pe, &lab5_px, &lab5_py, &lab5_pz, &lab5_M, &lab5_PIDe, &lab5_PIDmu, &lab5_PIDK, &lab5_PIDp, &lab5_ProbNNe, &lab5_ProbNNk, &lab5_ProbNNp, &lab5_ProbNNpi, &lab5_ProbNNmu, &lab5_ProbNNghost, &lab5_TRACK_CHI2NDOF, &lab5_gp, &lab5_pt, &lab5_chi2};

    unsigned int niswitchvars = 1;
    int* iswitchvars3[] = {&lab3_ID};
    int* iswitchvars5[] = {&lab5_ID};
    
    unsigned int nbswitchvars = 1;
    bool* bswitchvars3[] = {&lab3_isMuon};
    bool* bswitchvars5[] = {&lab5_isMuon};

    double lab0_ENDVERTEX_X, lab0_ENDVERTEX_Y, lab0_OWNPV_X, lab0_OWNPV_Y, lab0_RFD;
    double lab2_ENDVERTEX_X, lab2_ENDVERTEX_Y, lab2_OWNPV_X, lab2_OWNPV_Y, lab2_RFD;
    double lab345_MIN_PT, lab345_MIN_TRACK_CHI2, lab1345_TRACK_GhostProb;
    double lab0_ENDVERTEX_CHI2, lab0_VCHI2NDOF, lab2_ENDVERTEX_CHI2, lab2_VCHI2NDOF;
    int lab0_ENDVERTEX_NDOF, lab2_ENDVERTEX_NDOF;

    itree->SetBranchAddress("lab1_TRACK_GhostProb", &lab1_gp);
    itree->SetBranchAddress("lab3_TRACK_GhostProb", &lab3_gp);
    itree->SetBranchAddress("lab4_TRACK_GhostProb", &lab4_gp);
    itree->SetBranchAddress("lab5_TRACK_GhostProb", &lab5_gp);

    itree->SetBranchAddress("lab3_PT", &lab3_pt);
    itree->SetBranchAddress("lab4_PT", &lab4_pt);
    itree->SetBranchAddress("lab5_PT", &lab5_pt);

    itree->SetBranchAddress("lab3_IPCHI2_OWNPV", &lab3_chi2);
    itree->SetBranchAddress("lab4_IPCHI2_OWNPV", &lab4_chi2);
    itree->SetBranchAddress("lab5_IPCHI2_OWNPV", &lab5_chi2);

    itree->SetBranchAddress("lab3_CosTheta", &lab3_CosTheta); itree->SetBranchAddress("lab5_CosTheta", &lab5_CosTheta);
    itree->SetBranchAddress("lab3_IP_OWNPV", &lab3_IP_OWNPV); itree->SetBranchAddress("lab5_IP_OWNPV", &lab5_IP_OWNPV);
    itree->SetBranchAddress("lab3_P",        &lab3_p       ); itree->SetBranchAddress("lab5_P",        &lab5_p);
    itree->SetBranchAddress("lab3_PE",       &lab3_pe      ); itree->SetBranchAddress("lab5_PE",       &lab5_pe);
    itree->SetBranchAddress("lab3_PX",       &lab3_px      ); itree->SetBranchAddress("lab5_PX",       &lab5_px);
    itree->SetBranchAddress("lab3_PY",       &lab3_py      ); itree->SetBranchAddress("lab5_PY",       &lab5_py);
    itree->SetBranchAddress("lab3_PZ",       &lab3_pz      ); itree->SetBranchAddress("lab5_PZ",       &lab5_pz);
    itree->SetBranchAddress("lab3_M",        &lab3_M       ); itree->SetBranchAddress("lab5_M",        &lab5_M);
    itree->SetBranchAddress("lab3_ID",       &lab3_ID      ); itree->SetBranchAddress("lab5_ID",       &lab5_ID);
    itree->SetBranchAddress("lab3_PIDe",     &lab3_PIDe    ); itree->SetBranchAddress("lab5_PIDe",     &lab5_PIDe);
    itree->SetBranchAddress("lab3_PIDmu",    &lab3_PIDmu   ); itree->SetBranchAddress("lab5_PIDmu",    &lab5_PIDmu);
    itree->SetBranchAddress("lab3_PIDK",     &lab3_PIDK    ); itree->SetBranchAddress("lab5_PIDK",     &lab5_PIDK);
    itree->SetBranchAddress("lab3_PIDp",     &lab3_PIDp    ); itree->SetBranchAddress("lab5_PIDp",     &lab5_PIDp);

    itree->SetBranchAddress("lab3_ProbNNe",        &lab3_ProbNNe       ); itree->SetBranchAddress("lab5_ProbNNe",        &lab5_ProbNNe       );
    itree->SetBranchAddress("lab3_ProbNNk",        &lab3_ProbNNk       ); itree->SetBranchAddress("lab5_ProbNNk",        &lab5_ProbNNk       );
    itree->SetBranchAddress("lab3_ProbNNp",        &lab3_ProbNNp       ); itree->SetBranchAddress("lab5_ProbNNp",        &lab5_ProbNNp       );
    itree->SetBranchAddress("lab3_ProbNNpi",       &lab3_ProbNNpi      ); itree->SetBranchAddress("lab5_ProbNNpi",       &lab5_ProbNNpi      );
    itree->SetBranchAddress("lab3_ProbNNmu",       &lab3_ProbNNmu      ); itree->SetBranchAddress("lab5_ProbNNmu",       &lab5_ProbNNmu      );
    itree->SetBranchAddress("lab3_ProbNNghost",    &lab3_ProbNNghost   ); itree->SetBranchAddress("lab5_ProbNNghost",    &lab5_ProbNNghost   );
    itree->SetBranchAddress("lab3_isMuon",         &lab3_isMuon        ); itree->SetBranchAddress("lab5_isMuon",         &lab5_isMuon        );
    itree->SetBranchAddress("lab3_TRACK_CHI2NDOF", &lab3_TRACK_CHI2NDOF); itree->SetBranchAddress("lab5_TRACK_CHI2NDOF", &lab5_TRACK_CHI2NDOF);

    itree->SetBranchAddress("lab0_ENDVERTEX_X", &lab0_ENDVERTEX_X);
    itree->SetBranchAddress("lab0_ENDVERTEX_Y", &lab0_ENDVERTEX_Y);
    itree->SetBranchAddress("lab0_OWNPV_X", &lab0_OWNPV_X);
    itree->SetBranchAddress("lab0_OWNPV_Y", &lab0_OWNPV_Y);
    itree->SetBranchAddress("lab0_ENDVERTEX_CHI2", &lab0_ENDVERTEX_CHI2);
    itree->SetBranchAddress("lab0_ENDVERTEX_NDOF", &lab0_ENDVERTEX_NDOF);

    itree->SetBranchAddress("lab2_ENDVERTEX_X", &lab2_ENDVERTEX_X);
    itree->SetBranchAddress("lab2_ENDVERTEX_Y", &lab2_ENDVERTEX_Y);
    itree->SetBranchAddress("lab2_OWNPV_X", &lab2_OWNPV_X);
    itree->SetBranchAddress("lab2_OWNPV_Y", &lab2_OWNPV_Y);
    itree->SetBranchAddress("lab2_ENDVERTEX_CHI2", &lab2_ENDVERTEX_CHI2);
    itree->SetBranchAddress("lab2_ENDVERTEX_NDOF", &lab2_ENDVERTEX_NDOF);

    std::vector<double> ghost_prob, trk_pt, trk_chi2;
    std::vector<double> *pghost_prob(&ghost_prob), *ptrk_pt(&trk_pt), *ptrk_chi2(&trk_chi2);
    ghost_prob.reserve(4);
    trk_pt.reserve(3);
    trk_chi2.reserve(3);

    otree->Branch("ghost_prob", "std::vector<double>", &pghost_prob);
    otree->Branch("trk_pt", "std::vector<double>", &ptrk_pt);
    otree->Branch("trk_chi2", "std::vector<double>", &ptrk_chi2);

    otree->Branch("lab0_RFD", &lab0_RFD);
    otree->Branch("lab2_RFD", &lab2_RFD);
    otree->Branch("lab345_MIN_PT", &lab345_MIN_PT);
    otree->Branch("lab345_MIN_TRACK_CHI2", &lab345_MIN_TRACK_CHI2);
    otree->Branch("lab1345_TRACK_GhostProb", &lab1345_TRACK_GhostProb);
    otree->Branch("lab0_VCHI2NDOF", &lab0_VCHI2NDOF);
    otree->Branch("lab2_VCHI2NDOF", &lab2_VCHI2NDOF);

    time_t start_time = time(NULL);
    unsigned long bytes(0L);
    for(unsigned long i = 0; i < nentries; ++i) {
      bytes += itree->GetEntry(i);
      if (i%10000 == 0) std::cout << "@entry: " << i << "/" << nentries
				  << ", bytes read so far: " << bytes << ", ETA: "
                  << difftime(time(NULL), start_time) / ((1. * i) / nentries)
                  << " seconds" << std::endl;
      assert(lab1_gp > 10*std::numeric_limits<double>::epsilon()
	     or lab1_gp <= 0);

      // swap lab3 and lab5, if necessary
      if ((lab3_ID == 321 || lab3_ID == -321) && (lab5_ID == 211 || lab5_ID == -211)) {
          for (unsigned int si = 0; si < ndswitchvars; si++) {
              double tmp = *dswitchvars3[si];
              *dswitchvars3[si] = *dswitchvars5[si];
              *dswitchvars5[si] = tmp;
          }
          for (unsigned int si = 0; si < niswitchvars; si++) {
              int tmp = *iswitchvars3[si];
              *iswitchvars3[si] = *iswitchvars5[si];
              *iswitchvars5[si] = tmp;
          }
          for (unsigned int si = 0; si < nbswitchvars; si++) {
              bool tmp = *bswitchvars3[si];
              *bswitchvars3[si] = *bswitchvars5[si];
              *bswitchvars5[si] = tmp;
          }
      }

      ghost_prob.push_back(lab1_gp);
      ghost_prob.push_back(lab3_gp);
      ghost_prob.push_back(lab4_gp);
      ghost_prob.push_back(lab5_gp);

      trk_pt.push_back(lab3_pt);
      trk_pt.push_back(lab4_pt);
      trk_pt.push_back(lab5_pt);

      trk_chi2.push_back(lab3_chi2);
      trk_chi2.push_back(lab4_chi2);
      trk_chi2.push_back(lab5_chi2);

      lab0_RFD = TMath::Sqrt((lab0_ENDVERTEX_X - lab0_OWNPV_X) * (lab0_ENDVERTEX_X - lab0_OWNPV_X) + (lab0_ENDVERTEX_Y - lab0_OWNPV_Y) * (lab0_ENDVERTEX_Y - lab0_OWNPV_Y));
      lab2_RFD = TMath::Sqrt((lab2_ENDVERTEX_X - lab2_OWNPV_X) * (lab2_ENDVERTEX_X - lab2_OWNPV_X) + (lab2_ENDVERTEX_Y - lab2_OWNPV_Y) * (lab2_ENDVERTEX_Y - lab2_OWNPV_Y));
      lab345_MIN_PT = TMath::Min(TMath::Min(lab3_pt, lab4_pt), lab5_pt);
      lab345_MIN_TRACK_CHI2 = TMath::Min(TMath::Min(lab3_chi2, lab4_chi2), lab5_chi2);
      lab1345_TRACK_GhostProb = TMath::Max(TMath::Max(lab1_gp, lab3_gp), TMath::Max(lab4_gp, lab5_gp));
      lab0_VCHI2NDOF = lab0_ENDVERTEX_CHI2 / lab0_ENDVERTEX_NDOF;
      lab2_VCHI2NDOF = lab2_ENDVERTEX_CHI2 / lab2_ENDVERTEX_NDOF;

      otree->Fill();

      // cleanup
      ghost_prob.clear();
      trk_pt.clear();
      trk_chi2.clear();
    }

    // cleanup
    otree->Write(); delete otree; ofile.Close();
    delete itree; ifile.Close();
  }
  return 0;
}
