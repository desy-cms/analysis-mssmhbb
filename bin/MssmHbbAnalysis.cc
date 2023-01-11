#include "Analysis/MssmHbb/interface/MssmHbbAnalyser.h"

using namespace analysis;
using namespace analysis::tools;
using namespace analysis::mssmhbb;

int main(int argc, char ** argv)
{
    TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms

    MssmHbbAnalyser mssmhbb(argc,argv);
    
    auto config = mssmhbb.config();
    if ( config -> muonsVeto() )
    {
       mssmhbb.jetHistograms(config->nJetsMin(),"final_noveto");
       mssmhbb.jetHistograms(config->nJetsMin(),"final_muon_veto");
    }

    for ( int i = 0 ; i < mssmhbb.nEvents() ; ++i )
    {
        if ( ! mssmhbb.event(i)             )  continue;    // read event, run selection/json
        if ( ! mssmhbb.triggerSelection()   )  continue;    // trigger
        if ( ! mssmhbb.jetSelector()        )  continue;    // jet selector (remove fake jets)
        mssmhbb.HEMCorrection();    // HEM modules correction
        if ( ! mssmhbb.selectionPrimaryVertex()       )  continue;  // PV selection
        if ( ! mssmhbb.jetCorrections()     )  continue;    // jet corrections
        if ( ! mssmhbb.muonCorrections()    )  continue;    // muon corrections
        mssmhbb.fsrCorrections(mssmhbb.mainJets(),mssmhbb.fsrCandidates());
        if ( ! mssmhbb.jetSelection()       )  continue;    // jets
        if ( ! mssmhbb.muonSelector()       )  continue;    // muon selector (remove fake muons)
        if ( ! mssmhbb.muonSelection()      )  continue;    // muon selection
        if ( ! mssmhbb.muonJetSelection()   )  continue;    // muon-jet association
        if ( ! mssmhbb.btagSelection()      )  continue;    // btagging
        if ( config -> muonsVeto() ) mssmhbb.fillJetHistograms("final_noveto");
        if (   mssmhbb.muonVeto()           )  continue;    // muon veto (for full hadronic)
        if ( config -> muonsVeto() ) mssmhbb.fillJetHistograms("final_muon_veto");
        mssmhbb.actionApplyPrefiringWeight();
        mssmhbb.endSelection();
    }
} // end main
     



