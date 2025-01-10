#include "XSecAnalyzer/Selections/EventCategories1pi.hh"
#include "XSecAnalyzer/Selections/CC1mu1pi.hh"

CC1mu1pi::CC1mu1pi() : SelectionBase("CC1mu1pi") {}

void CC1mu1pi::define_category_map() {
  categ_map_ = CC1mu1pi_MAP;
};

void CC1mu1pi::define_output_branches() {
  // Save any additional variables to the output TTree
  this->set_branch( &sig_isNuMu_, "mc_is_numu");
  this->set_branch( &sig_is_CC_, "mc_is_cc");
  this->set_branch( &sig_inFV_, "mc_vertex_in_FV");
  this->set_branch( &sig_nPion_, "mc_npion");
  this->set_branch( &sig_has_fs_muon_, "mc_has_fs_muon");
  this->set_branch( &sig_is_signal_, "mc_is_signal");

  this->set_branch( &sel_reco_vertex_in_FV_, "sel_reco_vertex_in_FV");
  this->set_branch( &sel_pfp_starts_in_PCV_, "sel_pfp_starts_in_PCV");
  this->set_branch( &sel_has_muon_candidate_, "sel_has_muon_candidate");
  this->set_branch( &sel_presel_topo_cut_passed_, "sel_presel_topo_cut_passed");
  this->set_branch( &sel_cc1pi_topo_cut_passed_, "sel_cc1pi_topo_cut_passed");
  this->set_branch( &sel_phase_space_opening_angle_, "sel_phase_space_opening_angle");
  this->set_branch( &sel_phase_space_pion_mom_, "sel_phase_space_pion_mom");
  this->set_branch( &sel_phase_space_muon_mom_, "sel_phase_space_muon_mom");
  this->set_branch( &sel_good_bdt_scores_, "sel_good_bdt_scores");
  this->set_branch( &sel_reselected_muon_, "sel_reselected_muon_");

  this->set_branch( &sel_opening_angle_, "sel_opening_angle");
  this->set_branch( &sel_nu_mu_cc_, "sel_nu_mu_cc");
  this->set_branch( &sel_nu_mu_cc_1pi_, "sel_nu_mu_cc_1pi");
  this->set_branch( &sel_muon_contained_, "sel_muon_contained");
  this->set_branch( &sel_all_near_vertex_, "sel_all_near_vertex");
  this->set_branch( &max_vertex_distance_, "max_vertex_distance");
  this->set_branch( &sel_pion_dEdx_, "sel_pion_dEdx");
  this->set_branch( &sel_pion_dEdx_good_, "sel_pion_dEdx_good");
  this->set_branch( &sel_min_2_tracks_, "sel_min_2_tracks");
  this->set_branch( &sel_num_non_protons_good_, "sel_num_non_protons_good");
  this->set_branch( &sel_max_1_escape_, "sel_max_1_escape");
  this->set_branch( &sel_n_escape_, "sel_n_escape");
  this->set_branch( &sel_pion_not_in_gap_, "sel_pion_not_in_gap");
  this->set_branch( &sel_muon_not_in_gap_, "sel_muon_not_in_gap");

  this->set_branch( &muon_candidate_idx_, "muon_candidate_idx");
  this->set_branch( &pion_candidate_idx_, "pion_candidate_idx");

  //this->set_branch( reco_p3mu_, "reco_p3mu");
  //this->set_branch( reco_p3pi_, "reco_p3pi");
  this->set_branch( &reco_pi_mcs_mom_, "reco_pi_mcs_mom_");
  this->set_branch( &reco_pi_range_mom_, "reco_pi_range_mom_");
  this->set_branch( &reco_mu_mcs_mom_, "reco_mu_mcs_mom_");
  this->set_branch( &reco_mu_range_mom_, "reco_mu_range_mom_");

  this->set_branch( &reco_pi_dirx_ , "reco_pi_dirx");
  this->set_branch( &reco_pi_diry_ , "reco_pi_diry");
  this->set_branch( &reco_pi_dirz_ , "reco_pi_dirz");

  this->set_branch( &reco_pi_costheta_ , "reco_pi_costheta");
  this->set_branch( &reco_mu_costheta_ , "reco_mu_costheta");
  this->set_branch( &reco_pi_phi_ , "reco_pi_phi");
  this->set_branch( &reco_mu_phi_ , "reco_mu_phi");

  this->set_branch( &true_pi_costheta_ , "true_pi_costheta");
  this->set_branch( &true_mu_costheta_ , "true_mu_costheta");
  this->set_branch( &true_pi_phi_ , "true_pi_phi");
  this->set_branch( &true_mu_phi_ , "true_mu_phi");

  this->set_branch( &reco_mu_dirx_ , "reco_mu_dirx");
  this->set_branch( &reco_mu_diry_ , "reco_mu_diry");
  this->set_branch( &reco_mu_dirz_ , "reco_mu_dirz");

  this->set_branch( &true_pi_dirx_ , "true_pi_dirx");
  this->set_branch( &true_pi_diry_ , "true_pi_diry");
  this->set_branch( &true_pi_dirz_ , "true_pi_dirz");
  this->set_branch( &true_mu_dirx_ , "true_mu_dirx");
  this->set_branch( &true_mu_diry_ , "true_mu_diry");
  this->set_branch( &true_mu_dirz_ , "true_mu_dirz");
  this->set_branch( &true_mu_mom_ , "true_mu_mom");
  this->set_branch( &true_pi_mom_ , "true_pi_mom");

  this->set_branch( &true_n_proton_, "true_n_proton");
  this->set_branch( &true_n_proton_above_thresh_, "true_n_proton_above_thresh");

  //this->set_branch( &mc_p3mu_, "true_p3mu");
  //this->set_branch( &mc_p3pi_, "true_p3pi");

};

void CC1mu1pi::compute_reco_observables(AnalysisEvent* Event) {
  if (pion_candidate_idx_ != BOGUS_INDEX) {
    reco_pi_dirx_ = Event->track_dirx_->at(pion_candidate_idx_);
    reco_pi_diry_ = Event->track_diry_->at(pion_candidate_idx_);
    reco_pi_dirz_ = Event->track_dirz_->at(pion_candidate_idx_);

    reco_pi_costheta_ = reco_pi_dirz_;
    reco_pi_phi_ = std::acos(
        reco_pi_dirx_/
        sqrt(reco_pi_dirx_*reco_pi_dirx_ + reco_pi_diry_*reco_pi_diry_)
    );

    reco_pi_mcs_mom_ = Event->track_mcs_mom_mu_->at(pion_candidate_idx_);
    reco_pi_range_mom_ = Event->track_range_mom_mu_->at(pion_candidate_idx_);
  }

  if (muon_candidate_idx_ != BOGUS_INDEX) {
    reco_mu_dirx_ = Event->track_dirx_->at(muon_candidate_idx_);
    reco_mu_diry_ = Event->track_diry_->at(muon_candidate_idx_);
    reco_mu_dirz_ = Event->track_dirz_->at(muon_candidate_idx_);

    reco_mu_costheta_ = reco_mu_dirz_;
    reco_mu_phi_ = std::acos(
        reco_mu_dirx_/
        sqrt(reco_mu_dirx_*reco_mu_dirx_ + reco_mu_diry_*reco_mu_diry_)
    );

    reco_mu_mcs_mom_ = Event->track_mcs_mom_mu_->at(muon_candidate_idx_);
    reco_mu_range_mom_ = Event->track_range_mom_mu_->at(muon_candidate_idx_);
  }
};

void CC1mu1pi::compute_true_observables(AnalysisEvent* Event) {
  // Evaluate the true kinematic variables of interest

  // If there isn't one, we don't need to do anything else
  if (!sig_is_signal_) return;

  // Loop over the true final-state particles of the input event
  size_t num_fs_particles = Event->mc_nu_daughter_pdg_->size();

  // Find the final-state muon with the highest momentum, and store its true
  // 3-momentum in the owned TVector3 object.
  mc_p3mu_->SetXYZ(0., 0., 0.);
  mc_p3pi_->SetXYZ(0., 0., 0.);

  bool found_muon = false;
  bool found_pion = false;
  for (size_t f = 0; f < num_fs_particles; ++f) {

    int pdg = abs(Event->mc_nu_daughter_pdg_->at(f));
    float px = Event->mc_nu_daughter_px_->at(f);
    float py = Event->mc_nu_daughter_py_->at(f);
    float pz = Event->mc_nu_daughter_pz_->at(f);

    if (pdg == MUON) {

      found_muon = true;
      //float px = Event->mc_nu_daughter_px_->at(f);
      //float py = Event->mc_nu_daughter_py_->at(f);
      //float pz = Event->mc_nu_daughter_pz_->at(f);

      // TODO CHECK THIS
      // Replace the stored 3-momentum with the one for the current final-state
      // muon if the latter is larger
      TVector3 temp_p3mu(px, py, pz);
      if (temp_p3mu.Mag() > mc_p3mu_->Mag()) {
        *mc_p3mu_ = temp_p3mu;

        true_mu_mom_  = mc_p3mu_->Mag();
        true_mu_dirx_ = px/true_mu_mom_;
        true_mu_diry_ = py/true_mu_mom_;
        true_mu_dirz_ = pz/true_mu_mom_;

        true_mu_costheta_ = true_mu_dirz_;
        true_mu_phi_ = std::acos(
            true_mu_dirx_/
            sqrt(true_mu_dirx_*true_mu_dirx_ + true_mu_diry_*true_mu_diry_)
        );


      }

    } // final-state muon
    if (pdg == PI_PLUS) {

      found_pion = true;

      // Replace the stored 3-momentum with the one for the current final-state
      // pion if the latter is larger
      TVector3 temp_p3pi(px, py, pz);
      if (temp_p3pi.Mag() > mc_p3pi_->Mag()) {
        *mc_p3pi_ = temp_p3pi;
        true_pi_mom_  = mc_p3pi_->Mag();
        true_pi_dirx_ = px/true_pi_mom_;
        true_pi_diry_ = py/true_pi_mom_;
        true_pi_dirz_ = pz/true_pi_mom_;

        true_pi_costheta_ = true_pi_dirz_;
        true_pi_phi_ = std::acos(
            true_pi_dirx_/
            sqrt(true_pi_dirx_*true_pi_dirx_ + true_pi_diry_*true_pi_diry_)
        );

      }

    } // FS pion

    if (pdg == PROTON) {
      ++true_n_proton_;
      double p = sqrt(px*px + py*py + pz*pz);

      if (p > proton_thresh_) {
        ++true_n_proton_above_thresh_;
      }
    }
  } // loop over final-state particles

  if (!found_muon) {
    std::cout << "WARNING: Missing muon in MC signal event!\n";
  }
  if (!found_pion) {
    std::cout << "WARNING: Missing pion in MC signal event!\n";
  }
};

int  CC1mu1pi::categorize_event(AnalysisEvent* Event) {
  // Identify the event category of the selected event

  // Real data has a bogus true neutrino PDG code that is not one of the
  // allowed values (±12, ±14, ±16)
  int abs_mc_nu_pdg = std::abs(Event->mc_nu_pdg_);
  Event->is_mc_ = (abs_mc_nu_pdg == ELECTRON_NEUTRINO ||
    abs_mc_nu_pdg == MUON_NEUTRINO || abs_mc_nu_pdg == TAU_NEUTRINO);
  if (!Event->is_mc_) {
    return kUnknown;
  }

  // All events outside of the true fiducial volume should be categorized
  // as "out of fiducial volume"
  bool mcVertexInFV = point_inside_FV(this->true_FV(),
    Event->mc_nu_vx_, Event->mc_nu_vy_, Event->mc_nu_vz_);
  if (!mcVertexInFV) {
    return kOOFV;
  }

  bool isNC = (Event->mc_nu_ccnc_ == NEUTRAL_CURRENT);
  if (isNC) return kNC;

  //abs? TODO
  if (abs_mc_nu_pdg == ELECTRON_NEUTRINO) {
    return kNuECC;
  }
  if (abs_mc_nu_pdg != MUON_NEUTRINO) {
    return kOther;
  }

  int interaction_type = Event->mc_nu_interaction_type_;

  if (interaction_type == 0) {
    if (sig_nPion_ == 0) {
      return kNuMuCC0pi_CCQE;
    }
    else if (sig_nPion_ == 1) {
      return kNuMuCC1pi_CCQE;
    }
    else {
      return kNuMuCCNpi_CCQE;
    }
  }
  else if (interaction_type == 10) {
    if (sig_nPion_ == 0) {
      return kNuMuCC0pi_CCMEC;
    }
    else if (sig_nPion_ == 1) {
      return kNuMuCC1pi_CCMEC;
    }
    else {
      return kNuMuCCNpi_CCMEC;
    }
  }
  else if (interaction_type == 1) {
    if (sig_nPion_ == 0) {
      return kNuMuCC0pi_CCRES;
    }
    else if (sig_nPion_ == 1) {
      return kNuMuCC1pi_CCRES;
    }
    else {
      return kNuMuCCNpi_CCRES;
    }
  }
  else {
    if (sig_nPion_ == 0) {
      return kNuMuCC0pi_Other;
    }
    else if (sig_nPion_ == 1) {
      return kNuMuCC1pi_Other;
    }
    else {
      return kNuMuCCNpi_Other;
    }
  }

  // We shouldn't ever get here, but return "unknown" just in case
  return kUnknown;
};

void CC1mu1pi::reset() {

  sig_nPion_ = 0;
  sig_isNuMu_ = false;
  sig_is_CC_ = false;
  sig_is_signal_ = false;
  sig_has_fs_muon_ = false;
  sig_inFV_ = false;

  sel_pfp_starts_in_PCV_ = false;
  sel_reco_vertex_in_FV_ = false;
  sel_presel_topo_cut_passed_ = false;
  sel_has_muon_candidate_ = false;
  sel_nu_mu_cc_ = false;
  sel_muon_contained_ = false;
  sel_opening_angle_ = BOGUS;
  sel_phase_space_opening_angle_ = false;
  sel_phase_space_pion_mom_ = false;
  sel_phase_space_muon_mom_ = false;
  sel_cc1pi_topo_cut_passed_ = false;
  sel_all_near_vertex_ = false;
  sel_pion_dEdx_ = BOGUS;
  sel_pion_dEdx_good_ = false;
  sel_nu_mu_cc_1pi_ = false;
  sel_min_2_tracks_ = false;
  sel_num_non_protons_good_ = false;
  sel_max_1_escape_ = false;
  sel_n_escape_ = 0;
  sel_pion_not_in_gap_ = false;
  sel_muon_not_in_gap_ = false;

  muon_candidate_idx_ = BOGUS_INDEX;
  pion_candidate_idx_ = BOGUS_INDEX;

  mc_p3mu_->SetXYZ(BOGUS, BOGUS, BOGUS);
  mc_p3pi_->SetXYZ(BOGUS, BOGUS, BOGUS);

  reco_pi_mcs_mom_ = BOGUS;
  reco_pi_range_mom_ = BOGUS;
  reco_pi_dirx_ = BOGUS;
  reco_pi_diry_ = BOGUS;
  reco_pi_dirz_ = BOGUS;

  reco_pi_costheta_ = BOGUS;
  reco_mu_costheta_ = BOGUS;
  reco_pi_phi_ = BOGUS;
  reco_mu_phi_ = BOGUS;

  true_pi_costheta_ = BOGUS;
  true_mu_costheta_ = BOGUS;
  true_pi_phi_ = BOGUS;
  true_mu_phi_ = BOGUS;

  reco_mu_mcs_mom_ = BOGUS;
  reco_mu_range_mom_ = BOGUS;
  reco_mu_dirx_ = BOGUS;
  reco_mu_diry_ = BOGUS;
  reco_mu_dirz_ = BOGUS;
  max_vertex_distance_ = -1.*BOGUS;

  true_n_proton_ = 0;
  true_n_proton_above_thresh_ = 0;

  sel_good_bdt_scores_ = true;
  sel_reselected_muon_ = false;
}

bool CC1mu1pi::selection(AnalysisEvent* Event) {

  //std::cout << "Doing selection" << std::endl;


  // Apply the containment cut to the starting positions of all
  // reconstructed tracks and showers. Pass this cut by default.
  // L218 -- doc41264
  sel_pfp_starts_in_PCV_ = true;
  //Also check here that all are nearby
  sel_all_near_vertex_ = true;
  max_vertex_distance_ = -1.*BOGUS;

  // Loop over each PFParticle in the event
  for ( int p = 0; p < Event->num_pf_particles_; ++p ) {

    // Only check direct neutrino daughters (generation == 2)
    unsigned int generation = Event->pfp_generation_->at( p );
    if ( generation != 2 ) continue;

    // Use the track reconstruction results to get the start point for every
    // PFParticle for the purpose of verifying containment. We could in
    // principle differentiate between tracks and showers here, but track
    // information was used in the CC0piNp analysis (which rejected all showers
    // downstream of this cut anyway). For consistency, we therefore apply the
    // track reconstruction here unconditionally.
    float x = Event->track_startx_->at( p );
    float y = Event->track_starty_->at( p );
    float z = Event->track_startz_->at( p );

    // Verify that the start of the PFParticle lies within the containment
    // volume. See https://stackoverflow.com/a/2488507 for an explanation of
    // the use of &= here. Don't worry, it's type-safe since both operands are
    // bool.
    sel_pfp_starts_in_PCV_ &= point_inside_FV( PCV_, x, y, z );

    //Look at all of the gen-2 particles and see if they're within a 
    //good distance 
    float vertex_distance = Event->track_start_distance_->at(p);
    sel_all_near_vertex_ &= (Event->track_start_distance_->at(p) < 9.5);
    if (vertex_distance > max_vertex_distance_) max_vertex_distance_ = vertex_distance;
  }

  // Require selected events to have a reco vertex within the reco fiducial
  // volume
  // L220 -- doc41264
  sel_reco_vertex_in_FV_ = point_inside_FV( this->reco_FV(),
    Event->nu_vx_, Event->nu_vy_, Event->nu_vz_ );

  // Require the event to have a topological score that passes the cut
  // L221 -- doc41264
  sel_presel_topo_cut_passed_ = Event->topological_score_ > TOPO_SCORE_CUT;

  // Sets the sel_has_muon_candidate_ flag as appropriate. The threshold check
  // is handled later.
  std::vector<int> muon_candidate_indices;
  std::vector<float> muon_pid_scores;
  for ( int p = 0; p < Event->num_pf_particles_; ++p ) {
    // Only direct neutrino daughters (generation == 2) will be considered as
    // possible muon candidates
    unsigned int generation = Event->pfp_generation_->at( p );
    if (generation != 2) continue;

    float track_score = Event->pfp_track_score_->at( p );
    float start_dist = Event->track_start_distance_->at( p );
    float track_length = Event->track_length_->at( p );
    //float pid_score = Event->track_llr_pid_score_->at( p );
    float proton_chi2 = Event->track_chi2_proton_->at(p);
    float muon_chi2 = Event->track_chi2_muon_->at(p);

    //From values from Phil's docdb 41264
    if ((track_score > MUON_TRACK_SCORE_CUT) && // > .85
        (start_dist < MUON_VTX_DISTANCE_CUT) && // < 4cm
        (track_length > MUON_LENGTH_CUT) && // > 20cm
        (proton_chi2 > proton_chi2_cutval_) && // > 60
        (muon_chi2 < muon_chi2_cutval_) && // < 30
        (proton_chi2/muon_chi2 > promu_chi2_ratio_cutval_) // > 7
    ) { // TODO -- Check this
      muon_candidate_indices.push_back( p );
    }
  }

  //Flip this if there's a muon candidate
  size_t num_candidates = muon_candidate_indices.size();
  sel_has_muon_candidate_ = (num_candidates > 0);

  sel_nu_mu_cc_ = sel_reco_vertex_in_FV_ && sel_pfp_starts_in_PCV_
    && sel_has_muon_candidate_ && sel_presel_topo_cut_passed_;

  std::vector<std::pair<int, float>> contained_ids_bdt_scores;
  std::vector<int> escaping_ids;
  //Reselect the muon at this point
  //See doc 33809 p.47 figure 34
  for (int candidate_id = 0; candidate_id < Event->num_pf_particles_;
        ++candidate_id) {

    unsigned int generation = Event->pfp_generation_->at(candidate_id);
    if (generation != 2) continue;

    float endx = Event->track_endx_->at( candidate_id );
    float endy = Event->track_endy_->at( candidate_id );
    float endz = Event->track_endz_->at( candidate_id );
    bool end_contained = point_inside_FV( MCV_, endx, endy, endz );
    //save muon bdt score to check contained
    float bdt_score = Event->pfp_muon_bdt_responses_->at(candidate_id);
    if (!end_contained) {
      escaping_ids.push_back(candidate_id);
    }
    else if (bdt_score < BOGUS) {
      contained_ids_bdt_scores.push_back(
        {candidate_id, bdt_score}
      );
    }
  }

  sel_n_escape_ = escaping_ids.size();
  if (sel_n_escape_ == 1) {
    //Just use the escaper
    muon_candidate_idx_ = escaping_ids[0];
    sel_muon_contained_ = false;
  }
  //If no good BDT vals skip event
  else if (sel_n_escape_ == 0 && (contained_ids_bdt_scores.size() > 0)) {
    //All contained? Use highest muon bdt score
    std::sort(contained_ids_bdt_scores.begin(), contained_ids_bdt_scores.end(),
              [](auto & a, auto & b){return (a.second > b.second);});
    muon_candidate_idx_ = contained_ids_bdt_scores[0].first;
    sel_muon_contained_ = true;
  }
  else if (contained_ids_bdt_scores.size() == 0) {
    sel_good_bdt_scores_ = false;
  }

  sel_reselected_muon_ = (muon_candidate_idx_ != BOGUS_INDEX);

  //More than 1? Will get cut out
  sel_max_1_escape_ = (sel_n_escape_ <= 1);

  //will be used later
  sel_min_2_tracks_ = (Event->num_tracks_ >= 2);

  //size_t n_escape = 0;
  std::vector<int> proton_candidates, non_proton_candidates;

  //std::cout << "Doing pi sel" << std::endl;

  //Now do the Generic Pion Selection
  for (int p = 0; p < Event->num_pf_particles_; ++p) {

    //Skip the already-selected muon
    if (p == muon_candidate_idx_) continue;

    unsigned int generation = Event->pfp_generation_->at( p );
    if (generation != 2) continue;

    float proton_bdt_score = Event->pfp_proton_bdt_responses_->at(p);
    //Turn into config/parameter TODO 
    if (proton_bdt_score > -.06) proton_candidates.push_back(p);
    else non_proton_candidates.push_back(p);
  }

  //FIGURE OUT WHAT TO DO HERE. WE NEED MAX 2 NON PROTONS
  //What if there's a muon candidate that has proton bdt > cut???

  //Total is num pfp? or num tracks? TODO ASK
  size_t num_non_protons = (non_proton_candidates.size());
  //std::cout << "Got " << num_non_protons << " non protons and " <<
  //             proton_candidates.size() << " protons" << std::endl;
  if (num_non_protons == 1) {
    pion_candidate_idx_ = non_proton_candidates[0];
    sel_pion_not_in_gap_ = (
      (Event->pfp_hitsV_->at(pion_candidate_idx_) > 0) &&
      (Event->pfp_hitsU_->at(pion_candidate_idx_) > 0) &&
      (Event->pfp_hitsY_->at(pion_candidate_idx_) > 0)
    );
    //Set this good here. Technically it's called 2NonProtons in Andy's TN
    //but we already have one as the muon
    sel_num_non_protons_good_ = true;

    //Only u here?
    sel_pion_dEdx_ = Event->track_trunk_dEdx_u_->at(pion_candidate_idx_);
    sel_pion_dEdx_good_ = (sel_pion_dEdx_ > 1.0); //Make this a member TODO
  }

  //Check if the muon has hits on all three planes
  if (muon_candidate_idx_ != BOGUS_INDEX) {
    sel_muon_not_in_gap_ = (
      (Event->pfp_hitsV_->at(muon_candidate_idx_) > 0) &&
      (Event->pfp_hitsU_->at(muon_candidate_idx_) > 0) &&
      (Event->pfp_hitsY_->at(muon_candidate_idx_) > 0)
    );
  }

  if ((muon_candidate_idx_ != BOGUS_INDEX) &&
      (pion_candidate_idx_ != BOGUS_INDEX)) {
    float mu_dirx = Event->track_dirx_->at(muon_candidate_idx_);
    float mu_diry = Event->track_diry_->at(muon_candidate_idx_);
    float mu_dirz = Event->track_dirz_->at(muon_candidate_idx_);
    float pi_dirx = Event->track_dirx_->at(pion_candidate_idx_);
    float pi_diry = Event->track_diry_->at(pion_candidate_idx_);
    float pi_dirz = Event->track_dirz_->at(pion_candidate_idx_);
    sel_opening_angle_ = std::acos(
      mu_dirx*pi_dirx + mu_diry*pi_diry + mu_dirz*pi_dirz
    );
  }

  sel_phase_space_opening_angle_ = ((sel_opening_angle_ < 2.65) && //TODO -- make this a member
                             (sel_opening_angle_ < BOGUS));
  sel_cc1pi_topo_cut_passed_ = (Event->topological_score_ > cc1pi_topo_cutval_); // >.67


  sel_nu_mu_cc_1pi_ = (
    sel_nu_mu_cc_ &&
    sel_min_2_tracks_ && sel_max_1_escape_ &&
    sel_num_non_protons_good_ &&
    //sel_pion_dEdx_good_ && sel_phase_space_opening_angle_ &&
    sel_pion_not_in_gap_ && sel_muon_not_in_gap_ &&
    sel_cc1pi_topo_cut_passed_
  );
  return sel_nu_mu_cc_1pi_;
};

void CC1mu1pi::define_constants() {
  // Define reco & true fiducial volumes, alongside any other constants used
  // within selection cuts

  // FV definitions as in P.Detje's BNB CC1mu1pi analysis
  // (see docdb 41264) expanding upon Andy Smith's work
  // (docdb 33809)
  // x_min, x_max, y_min, y_max, z_min, z_max
  this->define_true_FV(10., 246.35, -106.5, 106.5, 10, 986.8);
  this->define_reco_FV(10., 246.35, -106.5, 106.5, 10, 986.8);

  PCV_.X_Min = 10.;
  PCV_.X_Max = 246.35;
  PCV_.Y_Min = -106.5;
  PCV_.Y_Max = 106.5;
  PCV_.Z_Min = 10.;
  PCV_.Z_Max = 1026.8;

  MCV_.X_Min = 5.;
  MCV_.X_Max = 251.35;
  MCV_.Y_Min = -111.5;
  MCV_.Y_Max = 111.5;
  MCV_.Z_Min = 5.;
  MCV_.Z_Max = 1031.8;

}

bool CC1mu1pi::define_signal(AnalysisEvent* Event) {
  // Determine whether or not the input event matches the signal definition
  // for this selection. This determination should be done only with MC truth
  // information.

  // Require signal events to be inside the true fiducial volume
  sig_inFV_ = point_inside_FV(
      this->true_FV(),
      Event->mc_nu_vx_,
      Event->mc_nu_vy_,
      Event->mc_nu_vz_);

  // Require an incident muon neutrino
  sig_isNuMu_ = (Event->mc_nu_pdg_ == MUON_NEUTRINO);
  sig_is_CC_ = (Event->mc_nu_ccnc_ == CHARGED_CURRENT);
  sig_nPion_ = 0;

  // Require a final-state muon
  // And 1 pion
  // ALSO ADD IN KINEMATICS
  for (size_t p = 0; p < Event->mc_nu_daughter_pdg_->size(); ++p) {
    
    //Use abs because we are nu/antinu agnostic
    int pdg = std::abs(Event->mc_nu_daughter_pdg_->at(p));
    if (pdg == MUON) {
      sig_has_fs_muon_ = true;
    }
    else if (pdg == PI_PLUS) {
      ++sig_nPion_;
    }
  }

  sig_is_signal_ = (
    sig_inFV_ && sig_isNuMu_ &&
    sig_is_CC_ &&
    sig_has_fs_muon_ && (sig_nPion_ == 1)
 );

  return sig_is_signal_;
}
