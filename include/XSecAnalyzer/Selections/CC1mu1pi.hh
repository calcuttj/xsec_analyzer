#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/SelectionBase.hh"

class CC1mu1pi : public SelectionBase {
 public:
  CC1mu1pi();
  ~CC1mu1pi(){};

  int CategorizeEvent(AnalysisEvent* Event);
  bool Selection(AnalysisEvent* Event);
  bool DefineSignal(AnalysisEvent* Event);
  void ComputeRecoObservables(AnalysisEvent* Event);
  void ComputeTrueObservables(AnalysisEvent* Event);
  void DefineOutputBranches();
  void DefineConstants();
  void DefineCategoryMap();

 private:

    float muon_chi2_cutval_ = 30.,
          proton_chi2_cutval_ = 60.,
          promu_chi2_ratio_cutval_ = 7.;

    bool sig_isNuMu_;
    bool sig_inFV_;
    bool sig_has_fs_muon_;
    bool sig_is_signal_;
    size_t sig_nPion_;

    bool sel_reco_vertex_in_FV_;
    bool sel_pfp_starts_in_PCV_;
    bool sel_has_muon_candidate_;
    bool sel_topo_cut_passed_;
    bool sel_nu_mu_cc_;
    bool sel_muon_contained_;

    int muon_candidate_idx_;
    int pion_candidate_idx_;

    MyPointer< TVector3 > reco_p3mu_;
    MyPointer< TVector3 > mc_p3mu_;
    MyPointer< TVector3 > reco_p3pi_;
    MyPointer< TVector3 > mc_p3pi_;

};
