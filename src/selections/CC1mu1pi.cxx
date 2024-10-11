#include "XSecAnalyzer/Selections/EventCategories1pi.hh"
#include "XSecAnalyzer/Selections/CC1mu1pi.hh"

CC1mu1pi::CC1mu1pi() : SelectionBase("CC1mu1pi") {}

void CC1mu1pi::DefineCategoryMap() {
  categ_map_ = CC1mu1pi_MAP;
};

void CC1mu1pi::DefineOutputBranches() {};

void CC1mu1pi::ComputeRecoObservables(AnalysisEvent* Event) {};
void CC1mu1pi::ComputeTrueObservables(AnalysisEvent* Event) {
  // Evaluate the true kinematic variables of interest

  // Check if there is a true final-state muon in this event
  bool has_true_muon = (sig_isNuMu_ && Event->mc_nu_ccnc_ == CHARGED_CURRENT);

  // If there isn't one, we don't need to do anything else
  if (!has_true_muon) return;

  // Loop over the true final-state particles of the input event
  size_t num_fs_particles = Event->mc_nu_daughter_pdg_->size();

  // Find the final-state muon with the highest momentum, and store its true
  // 3-momentum in the owned TVector3 object.
  mc_p3mu_->SetXYZ(0., 0., 0.);
  mc_p3pi_->SetXYZ(0., 0., 0.);

  bool found_muon = false;
  bool found_pion = false;
  for (size_t f = 0u; f < num_fs_particles; ++f) {
    //Abs??
    int pdg = Event->mc_nu_daughter_pdg_->at(f);
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

int  CC1mu1pi::CategorizeEvent(AnalysisEvent* Event) {
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
  bool mcVertexInFV = point_inside_FV(this->ReturnTrueFV(),
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

  if (this->IsEventMCSignal()) {
    // Categorize all numuCC events as "CC other" since we don't look at the
    // hadronic content in this selection
    // TODO: revisit this!
    return kNuMuCCOther;
  }

  // We shouldn't ever get here, but return "unknown" just in case
  return kUnknown;
};
bool CC1mu1pi::Selection(AnalysisEvent* Event) {return true;};

void CC1mu1pi::DefineConstants() {
  // Define reco & true fiducial volumes, alongside any other constants used
  // within selection cuts

  // FV definitions as in P.Detje's BNB CC1mu1pi analysis
  // (see docdb 41264) expanding upon Andy Smith's work
  // (docdb 33809)
  // x_min, x_max, y_min, y_max, z_min, z_max
  this->DefineTrueFV(10., 246.35, -106.5, 106.5, 10, 986.8);
  this->DefineRecoFV(10., 246.35, -106.5, 106.5, 10, 986.8);
}

bool CC1mu1pi::DefineSignal(AnalysisEvent* Event) {
  // Determine whether or not the input event matches the signal definition
  // for this selection. This determination should be done only with MC truth
  // information.

  // Require signal events to be inside the true fiducial volume
  sig_inFV_ = point_inside_FV(
      this->ReturnTrueFV(),
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
    //TODO CHECK THAT THIS NEEDS TO BE ABS OR IF ANTIMU IS ABS'D
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
