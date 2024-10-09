#include "XSecAnalyzer/Selections/EventCategoriesXp.hh"
#include "XSecAnalyzer/Selections/CC1mu1pi.hh"

CC1mu1pi::CC1mu1pi() : SelectionBase("CC1mu1pi") {}

void CC1mu1pi::DefineCategoryMap() {};
void CC1mu1pi::DefineOutputBranches() {};
void CC1mu1pi::ComputeRecoObservables(AnalysisEvent* Event) {};
void CC1mu1pi::ComputeTrueObservables(AnalysisEvent* Event) {};
int  CC1mu1pi::CategorizeEvent(AnalysisEvent* Event) {return 1;};
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
