#include "Analysis/MssmHbb/interface/MssmHbbAnalyser.h"

using namespace analysis;
using namespace analysis::tools;
using namespace analysis::mssmhbb;

int main(int argc, char ** argv)
{
    TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms

    MssmHbbAnalyser mssmhbb(argc,argv);

    for ( int i = 0 ; i < mssmhbb.nEvents() ; ++i )
    {
        if ( ! mssmhbb.event(i)             )  continue;    // read event, run selection/json
        if ( ! mssmhbb.triggerSelection()   )  continue;    // trigger
        if ( ! mssmhbb.jetSelector()        )  continue;    // jet selector (remove fake jets)
        if ( ! mssmhbb.jetCorrections()     )  continue;    // jet corrections
        if ( ! mssmhbb.jetSelection()       )  continue;    // jets
        if ( ! mssmhbb.muonSelector()       )  continue;    // muon selector (remove fake muons)
        if ( ! mssmhbb.muonSelection()      )  continue;    // muon selection
        if ( ! mssmhbb.muonJetSelection()   )  continue;    // muon-jet association
        if ( ! mssmhbb.btagSelection()      )  continue;    // btagging
        if (   mssmhbb.muonVeto()           )  continue;    // muon veto (for full hadronic)
        mssmhbb.endSelection();
    }
} // end main
     



