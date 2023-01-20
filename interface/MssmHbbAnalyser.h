#ifndef Analysis_MssmHbb_MssmHbbAnalyser_h
#define Analysis_MssmHbb_MssmHbbAnalyser_h 1

// -*- C++ -*-
//
// Package:    Analysis/MssmHbbAnalyser
// Class:      MssmHbbAnalyser
// 
/**\class Analysis MssmHbbAnalyser.cc Analysis/MssmHbb/src/MssmHbbAnalyser.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Roberval Walsh Bastos Rangel
//         Created:  Mon, 20 Oct 2014 14:24:08 GMT
//
//

// system include files
#include <memory>
#include <vector>
#include <string>
// 
// user include files

#include "Analysis/Tools/interface/Analyser.h"
#include "Analysis/Tools/interface/Jet.h"

//
// class declaration
//

namespace analysis {
   namespace mssmhbb {

      class MssmHbbAnalyser : public analysis::tools::Analyser {
         
         public:
            MssmHbbAnalyser();
            MssmHbbAnalyser(int argc, char ** argv);
           ~MssmHbbAnalyser();
           
            virtual bool event(const int &i);
//            virtual bool muonJet(const bool & swap = false);
            virtual bool muonJet();
            virtual bool onlineMuonMatching();
            void fillMssmHbbTree();
            void mssmHbbTree();
            void mssmHbbHistograms(const std::string & label="");
            void fillMssmHbbHistograms(const std::string & label="");
            bool triggerSelection();
            bool jetSelector();
            bool jetSelection();
            bool muonSelector();
            bool muonSelection();
            bool muonJetSelection();
            bool btagSelection();
            void btagEfficiencyWeight();
            bool muonVeto();
            bool endSelection();
            void sortedBTagScorePlus3Jets();
            std::vector< std::shared_ptr<analysis::tools::Jet> > mainJets();
            std::vector< std::shared_ptr<analysis::tools::Jet> > fsrCandidates();
            // ----------member data ---------------------------
         protected:
            bool do_tree_;
            std::shared_ptr<TTree> mssmhbb_tree_;
            double mbb_;
            double mbbw_;
               
         private:
               

      };
   }
}

#endif  // Analysis_MssmHbb_MssmHbbAnalyser_h
