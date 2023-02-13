#include "Analysis/MssmHbb/interface/MssmHbbAnalyser.h"

using namespace analysis;
using namespace analysis::tools;
using namespace analysis::mssmhbb;

int main(int argc, char ** argv)
{
    TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms

    MssmHbbAnalyser mssmhbb(argc,argv);
    
    auto config = mssmhbb.config();
    auto workflow = config->workflow();

    if ( config->muonsVeto() )
    {
        mssmhbb.jetHistograms(config->nJetsMin(),"NoMuonVeto");
        mssmhbb.jetHistograms(config->nJetsMin(),"MuonVeto");
    }
    // Jets in the analysis
    static constexpr int jet1 = 1;
    static constexpr int jet2 = 2;
    static constexpr int jet3 = 3;

    // Muon in the analysis
    static constexpr int muon1 = 1;

    // List of jets for btag matching
    std::vector<int> bjets_matched = {jet1, jet2, jet3};
    int max_bjets_matched = 2;

    // TODO: make two workflows, each should also have the semileptonic workflow (will run depending on the configuration file)
    // TODO: the semileptonic workflow will also be used for the muon veto in the full hadronic case
    // TODO: make two workflows: one with btag selection and another with btag weights
    if (workflow == 1) // nominal MSSM Hbb workflow with btag selection
    {
        for ( int i = 0 ; i < mssmhbb.nEvents() ; ++i )
        {
            if ( ! mssmhbb.event(i)                       )  continue;    // read event, run selection/json
            if ( ! mssmhbb.triggerSelection()             )  continue;
            if ( ! mssmhbb.selectionPrimaryVertex()       )  continue;
            if ( ! mssmhbb.selectionJetId()               )  continue;
            if ( ! mssmhbb.selectionJetPileupId()         )  continue;
            if ( ! mssmhbb.selectionNJets()               )  continue;
            mssmhbb.actionApplyBjetRegression();                          // b jet regression applied
            mssmhbb.actionApplyJER();                                     // jet energy resolution applied
            if ( ! mssmhbb.selectionJet(jet1)             )  continue;    // jet selection
            if ( ! mssmhbb.selectionJet(jet2)             )  continue;    // jet selection
            if ( ! mssmhbb.selectionDiJetMass(jet1,jet2)  )  continue;
            if ( ! mssmhbb.selectionJet(jet3)             )  continue;    // jet selection
            if ( ! mssmhbb.selectionJetDeta(jet1,jet2)    )  continue;    // jet deltaEta selection
            if ( ! mssmhbb.selectionJetDr(jet1,jet2)      )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionJetDr(jet1,jet3)      )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionJetDr(jet2,jet3)      )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionBJet(jet1)            )  continue;    mssmhbb.actionApplyBtagSF(jet1); // btag selection // btag scale factor
            if ( ! mssmhbb.selectionBJet(jet2)            )  continue;    mssmhbb.actionApplyBtagSF(jet2); // btag selection // btag scale factor                             
            if ( ! mssmhbb.selectionBJet(jet3)            )  continue;    mssmhbb.actionApplyBtagSF(jet3); // btag selection // btag scale factor                            
            if ( ! mssmhbb.onlineJetMatching(jet1)        )  continue;    // jet trg matching
            if ( ! mssmhbb.onlineJetMatching(jet2)        )  continue;    // jet trg matching
            mssmhbb.actionApplyJetOnlineSF(jet1);
            mssmhbb.actionApplyJetOnlineSF(jet2);
            if ( ! mssmhbb.onlineBJetMatching(bjets_matched, max_bjets_matched)) continue;    // bjet trigger matching
            // * Semileptonic - if muons not defined in the configuration, none of these have any effect
            if ( ! mssmhbb.selectionMuonId()              ) continue;      // muon id selection
            if ( ! mssmhbb.selectionNMuons()              ) continue;
            if ( ! mssmhbb.selectionMuons()               ) continue;
            // ! the next three lines order matters!!!
            if ( ! mssmhbb.onlineMuonMatching()           ) continue;      // if there is a trigger, keeps all muons matching to online
            if ( ! mssmhbb.muonJet()                      ) continue;      // keep only only one muon that is within a jet
            mssmhbb.actionApplyMuonOnlineSF(muon1);                        // if there is a trigger, applies the muon scale factor
            // ! ----
            // * Semileptonic - end
            mssmhbb.fsrCorrections(mssmhbb.mainJets(), mssmhbb.fsrCandidates()); // FSR better at the end, for it may bias the matching
            mssmhbb.actionApplyScaleCorrection("L1 prefiring");
            if ( config->muonsVeto() ) mssmhbb.fillJetHistograms("NoMuonVeto");
            if (   mssmhbb.muonVeto()                     )  continue;    // muon veto (for full hadronic)
            if ( config->muonsVeto() ) mssmhbb.fillJetHistograms("MuonVeto");
            mssmhbb.endSelection();
        } // end loop
    } // end workflow 1

   if (workflow == 2) // nominal MSSM Hbb workflow with btag weight
    {
        for ( int i = 0 ; i < mssmhbb.nEvents() ; ++i )
        {
            if ( ! mssmhbb.event(i)                       )  continue;    // read event, run selection/json
            if ( ! mssmhbb.triggerSelection()             )  continue;
            if ( ! mssmhbb.selectionPrimaryVertex()       )  continue;
            if ( ! mssmhbb.selectionJetId()               )  continue;
            if ( ! mssmhbb.selectionJetPileupId()         )  continue;
            if ( ! mssmhbb.selectionNJets()               )  continue;
            mssmhbb.actionApplyBjetRegression();                          // b jet regression applied
            mssmhbb.actionApplyJER();                                     // jet energy resolution applied
            if ( ! mssmhbb.selectionJet(jet1)             )  continue;    // jet selection
            if ( ! mssmhbb.selectionJet(jet2)             )  continue;    // jet selection
            if ( ! mssmhbb.selectionDiJetMass(jet1,jet2)  )  continue;
            if ( ! mssmhbb.selectionJet(jet3)             )  continue;    // jet selection
            if ( ! mssmhbb.selectionJetDeta(jet1,jet2)    )  continue;    // jet deltaEta selection
            if ( ! mssmhbb.selectionJetDr(jet1,jet2)      )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionJetDr(jet1,jet3)      )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionJetDr(jet2,jet3)      )  continue;    // jet deltaR selection
            mssmhbb.btagEfficiencyWeight();  // btag weight on jets 1,2
            if ( ! mssmhbb.selectionBJet(jet3)            )  continue;
            if ( ! mssmhbb.onlineJetMatching(jet1)        )  continue;    // jet trg matching
            if ( ! mssmhbb.onlineJetMatching(jet2)        )  continue;    // jet trg matching
            mssmhbb.actionApplyJetOnlineSF(jet1);
            mssmhbb.actionApplyJetOnlineSF(jet2);
            if ( ! mssmhbb.onlineBJetMatching(bjets_matched, max_bjets_matched)) continue;    // bjet trigger matching
            // * Semileptonic - if muons not defined in the configuration, none of these have any effect
            if ( ! mssmhbb.selectionMuonId()              ) continue;      // muon id selection
            if ( ! mssmhbb.selectionNMuons()              ) continue;
            if ( ! mssmhbb.selectionMuons()               ) continue;
            // ! the next three lines order matters!!!
            if ( ! mssmhbb.onlineMuonMatching()           ) continue;      // if there is a trigger, keeps all muons matching to online
            if ( ! mssmhbb.muonJet()                      ) continue;      // keep only only one muon that is within a jet
            mssmhbb.actionApplyMuonOnlineSF(muon1);                        // if there is a trigger, applies the muon scale factor
            // ! ----
            // * Semileptonic - end
            mssmhbb.fsrCorrections(mssmhbb.mainJets(), mssmhbb.fsrCandidates()); // FSR better at the end, for it may bias the matching
            mssmhbb.actionApplyScaleCorrection("L1 prefiring");
            if ( config->muonsVeto() ) mssmhbb.fillJetHistograms("NoMuonVeto");
            if (   mssmhbb.muonVeto()                     )  continue;    // muon veto (for full hadronic)
            if ( config->muonsVeto() ) mssmhbb.fillJetHistograms("MuonVeto");
            mssmhbb.endSelection();
        } // end loop
    } // end workflow 2


    // for ( int i = 0 ; i < mssmhbb.nEvents() ; ++i )
    // {
    //     if ( ! mssmhbb.event(i)             )  continue;    // read event, run selection/json
    //     if ( ! mssmhbb.triggerSelection()   )  continue;    // trigger
    //     if ( ! mssmhbb.jetSelector()        )  continue;    // jet selector (remove fake jets)
    //     if ( ! mssmhbb.jetCorrections()     )  continue;    // jet corrections
    //     if ( ! mssmhbb.jetSelection()       )  continue;    // jets
    //     if ( ! mssmhbb.muonSelector()       )  continue;    // muon selector (remove fake muons)
    //     if ( ! mssmhbb.muonSelection()      )  continue;    // muon selection
    //     if ( ! mssmhbb.muonJetSelection()   )  continue;    // muon-jet association
    //     if ( ! mssmhbb.btagSelection()      )  continue;    // btagging
    //     if ( config -> muonsVeto() ) mssmhbb.fillJetHistograms("final_noveto");
    //     if (   mssmhbb.muonVeto()           )  continue;    // muon veto (for full hadronic)
    //     if ( config -> muonsVeto() ) mssmhbb.fillJetHistograms("final_muon_veto");
    //     mssmhbb.actionApplyPrefiringWeight();
    //     mssmhbb.endSelection();
    // }
} // end main
     



