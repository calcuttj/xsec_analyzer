// ROOT includes
#include "TH1.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/EventCategories1pi.hh"

std::map< int, std::pair< std::string, int > > CC1mu1pi_MAP = {
  { kUnknown, { "Unknown", kGray } },
  { kNuMuCC1pi_CCQE, { "CCmu1pi (CCQE)", kBlue - 2 } },
  { kNuMuCC1pi_CCMEC, { "CCmu1pi (CCMEC)", kBlue - 6 } },
  { kNuMuCC1pi_CCRES, { "CCmu1pi (CCRES)", kBlue - 9 } },
  { kNuMuCC1pi_Other, { "CCmu1pi (Other)", kBlue - 10 } },

  { kNuMuCC0pi_CCQE, { "CCmu0pi (CCQE)", kBlue - 2 } },
  { kNuMuCC0pi_CCMEC, { "CCmu0pi (CCMEC)", kBlue - 6 } },
  { kNuMuCC0pi_CCRES, { "CCmu0pi (CCRES)", kBlue - 9 } },
  { kNuMuCC0pi_Other, { "CCmu0pi (Other)", kBlue - 10 } },

  { kNuMuCCNpi_CCQE, { "CCmuNpi (CCQE)", kBlue - 2 } },
  { kNuMuCCNpi_CCMEC, { "CCmuNpi (CCMEC)", kBlue - 6 } },
  { kNuMuCCNpi_CCRES, { "CCmuNpi (CCRES)", kBlue - 9 } },
  { kNuMuCCNpi_Other, { "CCmuNpi (Other)", kBlue - 10 } },

  //{ kNuMuCCNpi, { "#nu_{#mu} CCN#pi", kAzure - 2 } },
  { kNuMuCCOther, { "Other #nu_{#mu} CC", kAzure } },
  { kNuECC, { "#nu_{e} CC", kViolet } },
  { kNC, { "NC", kOrange } },
  { kOOFV, {"Out FV", kRed + 3 } },
  { kOther, { "Other", kRed + 1 } }
};
