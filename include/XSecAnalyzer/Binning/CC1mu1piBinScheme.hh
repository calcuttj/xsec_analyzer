#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Binning/BinSchemeBase.hh"

class CC1mu1piBinScheme : public BinSchemeBase {

  public:

    CC1mu1piBinScheme();
    virtual void DefineBlocks() override;
};
