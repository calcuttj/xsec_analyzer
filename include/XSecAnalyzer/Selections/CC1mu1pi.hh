#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/SelectionBase.hh"

class CC1mu1pi : public SelectionBase {
 public:
  CC1mu1pi();
  ~CC1mu1pi(){};

  int categorize_event(AnalysisEvent* Event);
  bool selection(AnalysisEvent* Event);
  bool define_signal(AnalysisEvent* Event);
  void compute_reco_observables(AnalysisEvent* Event);
  void compute_true_observables(AnalysisEvent* Event);
  void define_output_branches();
  void define_constants();
  void define_category_map();
  void reset();

 private:

    float muon_chi2_cutval_ = 30.,
          proton_chi2_cutval_ = 60.,
          promu_chi2_ratio_cutval_ = 7.,
          presel_topo_cutval_ = .25,
          cc1pi_topo_cutval_ = .67;

    bool sig_isNuMu_;
    bool sig_inFV_;
    bool sig_has_fs_muon_;
    bool sig_is_signal_;
    size_t sig_nPion_;
    bool min_2_tracks_;
    bool max_1_escape_;
    bool num_non_protons_good_;

    bool sel_reco_vertex_in_FV_;
    bool sel_pfp_starts_in_PCV_;
    bool sel_has_muon_candidate_;
    bool sel_cc1pi_topo_cut_passed_;
    bool muon_not_in_gap_;
    bool pion_not_in_gap_;
    bool sel_presel_topo_cut_passed_;
    bool sel_good_topo_score_;
    bool sel_good_opening_angle_;
    bool sel_nu_mu_cc_;
    bool sel_nu_mu_cc_1pi_;
    bool sel_muon_contained_;

    float sel_opening_angle_;
    int muon_candidate_idx_;
    int pion_candidate_idx_;

    MyPointer< TVector3 > reco_p3mu_;
    MyPointer< TVector3 > mc_p3mu_;
    MyPointer< TVector3 > reco_p3pi_;
    MyPointer< TVector3 > mc_p3pi_;

};
