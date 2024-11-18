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
          cc1pi_topo_cutval_ = .67,
          proton_thresh_ = .300;

    // "Proton containment volume" from https://arxiv.org/abs/2403.19574. We keep
    // this definition so that this toy selection matches the CC inclusive
    // portion of the selection used in Gardiner's CC0piNp analysis.
    FiducialVolume PCV_;

    // Muon Containment Volume  -- 5cm from all edges
    FiducialVolume MCV_;

    bool sig_isNuMu_;
    bool sig_is_CC_;
    bool sig_inFV_;
    bool sig_has_fs_muon_;
    bool sig_is_signal_;
    int sig_nPion_;
    bool sel_min_2_tracks_;
    bool sel_max_1_escape_;
    int sel_n_escape_;
    bool sel_num_non_protons_good_;
    bool sel_good_bdt_scores_;
    bool sel_reselected_muon_;

    bool sel_reco_vertex_in_FV_;
    bool sel_pfp_starts_in_PCV_;
    bool sel_has_muon_candidate_;
    bool sel_cc1pi_topo_cut_passed_;
    bool sel_muon_not_in_gap_;
    bool sel_pion_not_in_gap_;
    bool sel_presel_topo_cut_passed_;
    bool sel_good_topo_score_;
    bool sel_phase_space_opening_angle_;
    bool sel_phase_space_pion_mom_;
    bool sel_phase_space_muon_mom_;
    bool sel_nu_mu_cc_;
    bool sel_nu_mu_cc_1pi_;
    bool sel_muon_contained_;
    bool sel_all_near_vertex_;
    float sel_pion_dEdx_;
    bool sel_pion_dEdx_good_;

    float sel_opening_angle_;
    int muon_candidate_idx_;
    int pion_candidate_idx_;

    //MyPointer< TVector3 > reco_p3mu_;
    MyPointer< TVector3 > mc_p3mu_;
    //MyPointer< TVector3 > reco_p3pi_;
    MyPointer< TVector3 > mc_p3pi_;

    float reco_mu_mcs_mom_ = BOGUS;
    float reco_mu_range_mom_ = BOGUS;
    float reco_mu_dirx_ = BOGUS;
    float reco_mu_diry_ = BOGUS;
    float reco_mu_dirz_ = BOGUS;

    float reco_pi_mcs_mom_ = BOGUS;
    float reco_pi_range_mom_ = BOGUS;
    float reco_pi_dirx_ = BOGUS;
    float reco_pi_diry_ = BOGUS;
    float reco_pi_dirz_ = BOGUS;

    float reco_pi_theta_ = BOGUS;
    float reco_pi_phi_ = BOGUS;
    float reco_mu_theta_ = BOGUS;
    float reco_mu_phi_ = BOGUS;

    float true_pi_mom_ = BOGUS;
    float true_pi_dirx_ = BOGUS;
    float true_pi_diry_ = BOGUS;
    float true_pi_dirz_ = BOGUS;
    float true_mu_mom_ = BOGUS;
    float true_mu_dirx_ = BOGUS;
    float true_mu_diry_ = BOGUS;
    float true_mu_dirz_ = BOGUS;

    int true_n_proton_ = BOGUS_INT;
    int true_n_proton_above_thresh_ = BOGUS_INT;

    float max_vertex_distance_ = -1.*BOGUS;
};
