#include "XSecAnalyzer/Selections/EventCategories1pi.hh"
#include "XSecAnalyzer/Selections/CC1mu1pi.hh"

CC1mu1pi::CC1mu1pi() : SelectionBase("CC1mu1pi") {}

void CC1mu1pi::define_category_map() {
  categ_map_ = CC1mu1pi_MAP;
};

void CC1mu1pi::define_output_branches() {
  // Save any additional variables to the output TTree
  this->set_branch( &sig_isNuMu_, "mc_is_numu");
  this->set_branch( &sig_inFV_, "mc_vertex_in_FV");
  this->set_branch( &sig_nPion_, "mc_npion");
  this->set_branch( &sig_has_fs_muon_, "mc_has_fs_muon");
  this->set_branch( &sig_is_signal_, "mc_is_signal");
  this->set_branch( &sel_reco_vertex_in_FV_, "sel_reco_vertex_in_FV");
  this->set_branch( &sel_pfp_starts_in_PCV_, "sel_pfp_starts_in_PCV");
  this->set_branch( &sel_has_muon_candidate_, "sel_has_muon_candidate");
  this->set_branch( &sel_topo_cut_passed_, "sel_topo_cut_passed");
  this->set_branch( &sel_nu_mu_cc_, "sel_nu_mu_cc");
  this->set_branch( &sel_muon_contained_, "sel_muon_contained");
  this->set_branch( &muon_candidate_idx_, "muon_candidate_idx");
  this->set_branch( &pion_candidate_idx_, "pion_candidate_idx");
  this->set_branch( reco_p3mu_, "reco_p3_mu");
  this->set_branch( reco_p3pi_, "reco_p3_pi");
  this->set_branch( mc_p3mu_, "true_p3_mu");
  this->set_branch( mc_p3pi_, "true_p3_pi");
};

void CC1mu1pi::compute_reco_observables(AnalysisEvent* Event) {

};

void CC1mu1pi::compute_true_observables(AnalysisEvent* Event) {
  // Evaluate the true kinematic variables of interest

  // Check if there is a true final-state muon in this event
  bool is_signal = (
    sig_isNuMu_ && (Event->mc_nu_ccnc_ == CHARGED_CURRENT) &&
    (sig_nPion_ == 1)
  );

  // If there isn't one, we don't need to do anything else
  if (!is_signal) return;

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
    if (pdg == MUON) {

      found_muon = true;
      float px = Event->mc_nu_daughter_px_->at(f);
      float py = Event->mc_nu_daughter_py_->at(f);
      float pz = Event->mc_nu_daughter_pz_->at(f);

      // Replace the stored 3-momentum with the one for the current final-state
      // muon if the latter is larger
      TVector3 temp_p3mu(px, py, pz);
      if (temp_p3mu.Mag() > mc_p3mu_->Mag()) {
        *mc_p3mu_ = temp_p3mu;
      }

    } // final-state muon
    if (pdg == PI_PLUS) {

      found_pion = true;
      float px = Event->mc_nu_daughter_px_->at(f);
      float py = Event->mc_nu_daughter_py_->at(f);
      float pz = Event->mc_nu_daughter_pz_->at(f);

      // Replace the stored 3-momentum with the one for the current final-state
      // pion if the latter is larger
      TVector3 temp_p3pi(px, py, pz);
      if (temp_p3pi.Mag() > mc_p3pi_->Mag()) {
        *mc_p3pi_ = temp_p3pi;
      }

    } // FS pion
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

  if (Event->mc_nu_pdg_ == ELECTRON_NEUTRINO) {
    return kNuECC;
  }
  if (!(Event->mc_nu_pdg_ == MUON_NEUTRINO)) {
    return kOther;
  }

  if (this->is_event_mc_signal()) {
    // Categorize all numuCC events as "CC other" since we don't look at the
    // hadronic content in this selection
    // TODO: revisit this!
    return kNuMuCCOther;
  }

  // We shouldn't ever get here, but return "unknown" just in case
  return kUnknown;
};
bool CC1mu1pi::selection(AnalysisEvent* Event) {



  // "Proton containment volume" from https://arxiv.org/abs/2403.19574. We keep
  // this definition so that this toy selection matches the CC inclusive
  // portion of the selection used in Gardiner's CC0piNp analysis.
  FiducialVolume PCV;
  PCV.X_Min = 10.;
  PCV.X_Max = 246.35;
  PCV.Y_Min = -106.5;
  PCV.Y_Max = 106.5;
  PCV.Z_Min = 10.;
  PCV.Z_Max = 1026.8;


  // Apply the containment cut to the starting positions of all
  // reconstructed tracks and showers. Pass this cut by default.
  sel_pfp_starts_in_PCV_ = true;

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
    sel_pfp_starts_in_PCV_ &= point_inside_FV( PCV, x, y, z );
  }

  // Require selected events to have a reco vertex within the reco fiducial
  // volume
  sel_reco_vertex_in_FV_ = point_inside_FV( this->reco_FV(),
    Event->nu_vx_, Event->nu_vy_, Event->nu_vz_ );

  // Require the event to have a topological score that passes the cut
  sel_topo_cut_passed_ = Event->topological_score_ > TOPO_SCORE_CUT;

  // TODO -- PCV? https://github.com/uboone/xsec_analyzer/blob/tutorial-umn/src/selections/TutorialCC1mu.cxx#L177

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
    float proton_bdt_score = Event->track_llr_pid_score_->at( p );
    //float proton_chi2 = Event->track_chi2_proton_->at(p);
    //float muon_chi2 = Event->track_chi2_muon_->at(p);

    //From values from Phil's docdb 41264
    if ((track_score > MUON_TRACK_SCORE_CUT) && // > .85
        (start_dist < MUON_VTX_DISTANCE_CUT) && // < 4cm
        (track_length > MUON_LENGTH_CUT) //&& // > 20cm
        //(proton_chi2 > proton_chi2_cutval_) && // > 60
        //(muon_chi2 < muon_chi2_cutval_) && // < 30
        //(proton_chi2/muon_chi2 > promu_chi2_ratio_cutval_) // > 7
    ) { // TODO -- Check this
      muon_candidate_indices.push_back( p );
      //muon_pid_scores.push_back( pid_score );

      ///if (pid_score > )

      // Check whether the muon candidate is contained. Use the PCV as the
      // containment volume.
      sel_muon_contained_ = false;
      float endx = Event->track_endx_->at( muon_candidate_idx_ );
      float endy = Event->track_endy_->at( muon_candidate_idx_ );
      float endz = Event->track_endz_->at( muon_candidate_idx_ );
      bool end_contained = point_inside_FV( PCV, endx, endy, endz );
      /*if (end_contained) {
        contained_muon_candidates.push_back(p);
      }
      else {
        escaping_muon_candidates.push_back(p);
      }*/

      //muon_proton_chi2s.push_back(proton_chi2);
      //muon_muon_chi2s.push_back(muon_chi2);
    }
  }

  size_t num_candidates = muon_candidate_indices.size();
  sel_has_muon_candidate_ = (num_candidates > 0);

  if ( num_candidates == 1 ) {
    muon_candidate_idx_ = muon_candidate_indices.front();
  }
  else if ( num_candidates > 1 ) {
    // In the case of multiple muon candidates, choose the one with the highest
    // PID score (most muon-like) as the one to use
    float highest_score = LOW_FLOAT;
    int chosen_index = BOGUS_INDEX;
    for ( size_t c = 0; c < num_candidates; ++c ) {
      float score = muon_pid_scores.at( c );
      if ( highest_score < score ) {
        highest_score = score;
        chosen_index = muon_candidate_indices.at( c );
      }
    }
    muon_candidate_idx_ = chosen_index;
  }
  else {
    muon_candidate_idx_ = BOGUS_INDEX;
  }


  // Check whether the muon candidate is contained. Use the PCV as the
  // containment volume.
  sel_muon_contained_ = false;
  if ( muon_candidate_idx_ != BOGUS_INDEX ) {
    float endx = Event->track_endx_->at( muon_candidate_idx_ );
    float endy = Event->track_endy_->at( muon_candidate_idx_ );
    float endz = Event->track_endz_->at( muon_candidate_idx_ );
    bool end_contained = point_inside_FV( PCV, endx, endy, endz );
    sel_muon_contained_ = end_contained;
  }

  sel_nu_mu_cc_ = sel_reco_vertex_in_FV_ && sel_pfp_starts_in_PCV_
    && sel_has_muon_candidate_ && sel_topo_cut_passed_;

  return sel_nu_mu_cc_;
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
    sig_inFV_ && sig_isNuMu_ && sig_has_fs_muon_ && (sig_nPion_ == 1)
 );

  return sig_is_signal_;
}
