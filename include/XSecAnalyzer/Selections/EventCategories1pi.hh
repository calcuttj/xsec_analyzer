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

  // True numu CC event with more than one final-state pion above threshold
  kNuMuCCNpi_CCQE = 5,
  kNuMuCCNpi_CCMEC = 6,
  kNuMuCCNpi_CCRES = 7,
  kNuMuCCNpi_Other = 8,

  // True numu CC event with 0 final-state pion above threshold
  kNuMuCC0pi_CCQE = 9,
  kNuMuCC0pi_CCMEC = 10,
  kNuMuCC0pi_CCRES = 11,
  kNuMuCC0pi_Other = 12,

  // Any true numu CC event which does not satisfy the criteria for inclusion
  // in one of the other categories above
  kNuMuCCOther = 13,

  // True nue CC event
  kNuECC = 14,

  // True neutral current event for any neutrino flavor
  kNC = 15,

  // True neutrino vertex (any reaction mode and flavor combination) is outside
  // of the fiducial volume
  kOOFV = 16,

  // All events that do not fall within any of the other categories (e.g.,
  // numubar CC)
  kOther = 17,
};

extern std::map< int, std::pair< std::string, int > > CC1mu1pi_MAP;
