
// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/CC1mu1p0pi.hh"
#include "XSecAnalyzer/Selections/CC1mu2p0pi.hh"
#include "XSecAnalyzer/Selections/CC1muNp0pi.hh"
#include "XSecAnalyzer/Selections/CC1mu1pi.hh"
#include "XSecAnalyzer/Selections/CC1mu1piNG2.hh"
#include "XSecAnalyzer/Selections/DummySelection.hh"
#include "XSecAnalyzer/Selections/SelectionFactory.hh"

SelectionFactory::SelectionFactory() {
}

SelectionBase* SelectionFactory::CreateSelection(
  const std::string& selection_name )
{
  SelectionBase* sel;
  if ( selection_name == "CC1mu1p0pi" ) {
    sel = new CC1mu1p0pi;
  }
  else if ( selection_name == "CC1mu2p0pi" ) {
    sel = new CC1mu2p0pi;
  }
  else if ( selection_name == "CC1muNp0pi" ) {
    sel = new CC1muNp0pi;
  }
  else if ( selection_name == "CC1mu1pi" ) {
    sel = new CC1mu1pi;
  }
  else if ( selection_name == "CC1mu1piNG2" ) {
    sel = new CC1mu1piNG2;
  }
  else if ( selection_name == "Dummy" ) {
    sel = new DummySelection;
  }
  else {
    std::cerr << "Selection name requested: " << selection_name
      << " is not implemented in " << __FILE__ << '\n';
    throw;
  }

  return sel;
}
