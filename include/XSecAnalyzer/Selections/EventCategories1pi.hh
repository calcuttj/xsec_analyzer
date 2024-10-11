#pragma once

// Standard library includes
#include <map>
#include <string>

// Enum used to label event categories of interest for analysis plots in
// the CC1pi analyses
enum EventCategory1pi {

  // Unable to categorize (e.g., because the event is real data and thus
  // has no MC truth information)
  kUnknown = 0,

  // Signal events broken down by underlying reaction mode
  kNuMuCC1pi_CCQE = 1,
  kNuMuCC1pi_CCMEC = 2,
  kNuMuCC1pi_CCRES = 3,
  kNuMuCC1pi_Other = 4,

  // True numu CC event with at least one final-state pion above threshold
  kNuMuCCNpi = 5,

  // Any true numu CC event which does not satisfy the criteria for inclusion
  // in one of the other categories above
  kNuMuCCOther = 6,

  // True nue CC event
  kNuECC = 7,

  // True neutral current event for any neutrino flavor
  kNC = 8,

  // True neutrino vertex (any reaction mode and flavor combination) is outside
  // of the fiducial volume
  kOOFV = 9,

  // All events that do not fall within any of the other categories (e.g.,
  // numubar CC)
  kOther = 10,
};

extern std::map< int, std::pair< std::string, int > > CC1mu1pi_MAP;
