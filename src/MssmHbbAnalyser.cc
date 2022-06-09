/**\class MssmHbb MssmHbbAnalyser.cc Analysis/Tools/src/MssmHbbAnalyser.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Roberval Walsh
//         Created:  Mon, 20 Oct 2014 14:24:08 GMT
//
//

// system include files
#include <iostream>
// 
// user include files
#include "Analysis/MssmHbb/interface/MssmHbbAnalyser.h"
#include "Analysis/Tools/interface/Composite.h"


//
// class declaration
//

using namespace analysis;
using namespace analysis::tools;
using namespace analysis::mssmhbb;

//
// constructors and destructor
//
MssmHbbAnalyser::MssmHbbAnalyser()
{
}

MssmHbbAnalyser::MssmHbbAnalyser(int argc, char ** argv) : BaseAnalyser(argc,argv), Analyser(argc,argv)
{
//   histograms("cutflow");
//   histograms("jet",config_->nJetsMin());
   do_tree_ = config_->doTree();
   if ( do_tree_ ) this -> mssmHbbTree();
   mbb_ = -1.;
   this->jetHistograms(config_->nJetsMin(),"final_selection");
   this->mssmHbbHistograms();
}

MssmHbbAnalyser::~MssmHbbAnalyser()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   
   if ( do_tree_ ) mssmhbb_tree_ -> Write();
}


//
// member functions
//
// ------------ method called for each event  ------------
bool MssmHbbAnalyser::event(const int & i)
{
   // parent function checks only json and run range validity
   if ( ! Analyser::event(i) ) return false;
      
   return true;
}bool MssmHbbAnalyser::triggerSelection()
{
   return ( this->selectionL1() && this->selectionHLT() );
}

bool MssmHbbAnalyser::jetSelection()
{
      if ( ! this->selectionNJets()         )   return false;
      // kinematics selection
      for ( int i = 1; i <= config_->nJetsMin(); ++i )
      {
         if ( ! this->selectionJet(i)          )   return false;
         if ( i <= 2  && (! this->onlineJetMatching(i))     )   return false;
      }
      // delta R selection
      for ( int i = 1; i <= config_->nJetsMin(); ++i )
      {
         for ( int j = i+1; j <= config_->nJetsMin(); ++j )
            if ( ! this->selectionJetDr(i,j)      )   return false;
      }
      
      // delta eta jet selection
      if ( ! this->selectionJetDeta(1,2)    )   return false;
      
      return true;
      
}
bool MssmHbbAnalyser::btagSelection()
{
   if ( ! this->onlineBJetMatching(1)    )   return false;
   if ( ! this->onlineBJetMatching(2)    )   return false;
   if ( config_->btagEfficiencies(1) == "" )
   {
      if ( ! this->selectionBJet(1)      )   return false;
      if ( ! this->selectionBJet(2)      )   return false;
      if ( ! this->selectionBJet(3)      )   return false;
      this->actionApplyBtagSF(1);
      this->actionApplyBtagSF(2);
      this->actionApplyBtagSF(3);
   }
   else
   {
      this->actionApplyBtagEfficiency(1,1);
      this->actionApplyBtagEfficiency(2,1);
      if ( config_->signalRegion() )
      {
         this->actionApplyBtagEfficiency(3,2);
      }
      else
      {
         if ( ! this->selectionBJet(3)      )   return false;
         this->actionApplyBtagSF(3);
      }
   }
   return true;
      
}
bool MssmHbbAnalyser::muonSelection()
{
      if ( ! this->selectionNMuons()        )   return false;
      if ( ! this->selectionMuons()         )   return false;
      if ( ! this->onlineMuonMatching()     )   return false;
      return true;
      
}
bool MssmHbbAnalyser::endSelection()
{
      this->fillJetHistograms("final_selection");
      this->fillMssmHbbHistograms();
      this->fillMssmHbbTree();
      return true;
      
}

bool MssmHbbAnalyser::muonJet(const bool & swap)
{
   // jet 1 is the muon jet, swap with jet 2 in case jet1 does not have a muon 
   // Muon jet association, a muon can be in either of the two leading jets.
   // If argument is false(default, the jet label follows the jet ranking, regardless of where the muon is
   // If argument is true, the jet label 1 refers to the muon-jet: *** the order of the workflow may matter!!! ***
   if ( ! muonsanalysis_ || ! jetsanalysis_ ) return true;  // will skip this selection
   
   int r1 = 1;
   int r2 = 2;
   int j1 = r1-1;
   int j2 = r2-1;
   ++ cutflow_;
   if ( std::string(h1_["cutflow"] -> GetXaxis()-> GetBinLabel(cutflow_+1)) == "" ) 
   {
      std::string label = Form("MSSMHbb: Jet-muon association (delta_R < %4.2f)", config_->jetsMuonsDRMax());;
      std::string label_swap = Form("%s | Muon-Jet index 1 (swap)", label.c_str());;
      if ( swap ) h1_["cutflow"] -> GetXaxis()-> SetBinLabel(cutflow_+1,label_swap.c_str());
      else        h1_["cutflow"] -> GetXaxis()-> SetBinLabel(cutflow_+1,label.c_str());
   }
   
   if ( selectedJets_.size() < 2 ||  selectedMuons_.size() < 1 ) return false;
   
   auto jet1 = selectedJets_[j1];
   jet1 -> addMuon(selectedMuons_,config_->jetsMuonsDRMax());
   auto jet2 = selectedJets_[j2];
   jet2 -> addMuon(selectedMuons_,config_->jetsMuonsDRMax());

   // Only muons within those 2 jets are kept
   std::vector< std::shared_ptr<Muon> > sel_muons;
   sel_muons.clear();
   for ( auto & sm : selectedMuons_ )
   {
      if ( sm == jet1->muon() || sm == jet2->muon() )
         sel_muons.push_back(sm);
   }
   selectedMuons_ = sel_muons;
   
      
   if ( selectedMuons_.size() == 0 ) return false;
   if ( !  jet1 -> muon() && swap ) this->jetSwap(r1,r2);
   
   h1_["cutflow"] -> Fill(cutflow_,weight_);
   return true;
   
}

bool MssmHbbAnalyser::onlineMuonMatching()
{
   if ( ! MuonAnalyser::onlineMuonMatching() ) return false;
   return true;
      
}

void MssmHbbAnalyser::fillMssmHbbTree()
{
   if ( ! do_tree_ ) return;
   ++ cutflow_;
   if ( std::string(h1_["cutflow"] -> GetXaxis()-> GetBinLabel(cutflow_+1)) == "" ) 
      h1_["cutflow"] -> GetXaxis()-> SetBinLabel(cutflow_+1,"Fill MssmHbb tree");
   
   Composite<Jet,Jet> c_12(*(selectedJets_[0]),*(selectedJets_[1]));
   if ( config_->isMC() || !config_->signalRegion() ) mbb_ = c_12.m();
   
   mbbw_ = weight_;
   
   mssmhbb_tree_ -> Fill();
   h1_["cutflow"] -> Fill(cutflow_,weight_);

}

void MssmHbbAnalyser::mssmHbbTree()
{
//   do_tree_ = true;
   this->output()->cd();
   mssmhbb_tree_ = std::make_shared<TTree>("mssmhbb","TTree with mbb and weight for FitModel");
   mssmhbb_tree_ -> Branch("mbb",&mbb_,"mbb/D");
   mssmhbb_tree_ -> Branch("weight",&mbbw_,"weight/D");
}

void MssmHbbAnalyser::fillMssmHbbHistograms(const std::string & label)
{
   this->output()->cd();
   std::string cf_label = "Fill MssmHbb Histograms";
   std::string prefix = "mssmhbb_";
   if ( label != "" )
   {
      cf_label = Form("Fill MssmHbb Histograms - dir = %s", label.c_str());
      this->output()->cd(label.c_str());
      prefix = "";
   }
   
   auto njets = config_->nJetsMin();
   std::string name;
   std::string pname;
   for ( int r = 1 ; r <= njets ; ++r )
   {
      int j1 = r-1;
      auto jet_j1 = selectedJets_[j1];
      
      name  = Form("pt_jet%d",r);  pname = Form("%s%s",prefix.c_str(),name.c_str());
      h1_[pname] -> Fill(jet_j1->pt(),weight_);
      name  = Form("eta_jet%d",r);  pname = Form("%s%s",prefix.c_str(),name.c_str());
      h1_[pname] -> Fill(jet_j1->eta(),weight_);
      name  = Form("btag_jet%d",r);  pname = Form("%s%s",prefix.c_str(),name.c_str());
      h1_[pname] -> Fill(JetAnalyser::btag(*jet_j1,config_->btagAlgorithm()),weight_);
      for ( int r2 = r+1 ; r2 <= njets ; ++r2 )
      {
         if ( r2 == r ) continue;
         int j2 = r2-1;
         auto jet_j2 = selectedJets_[j2];
         float mxy = 0.;
         Composite<Jet,Jet> c_xy(*(jet_j1),*(jet_j2));
         if ( config_->isMC() || !config_->signalRegion() ) mxy = c_xy.m();
         name  = Form("m_jet%d%d",r,r2);  pname = Form("%s%s",prefix.c_str(),name.c_str());
         h1_[pname] -> Fill(mxy,weight_);
         
      } 
   }
   // deprecated
   Composite<Jet,Jet> c_12(*(selectedJets_[0]),*(selectedJets_[1]));
   float mbb = 0.;
   if ( config_->isMC() || !config_->signalRegion() ) mbb = c_12.m();
   name  = "mbb";  pname = Form("%s%s",prefix.c_str(),name.c_str());
   h1_[pname] -> Fill(mbb,weight_);
   //
   
   cutflow(cf_label);

}


void MssmHbbAnalyser::mssmHbbHistograms(const std::string & label)
{
   this->output()->cd();
   if ( label != "" )  this->output()->cd(label.c_str());
   
   auto njets = config_->nJetsMin();
   std::string prefix = "mssmhbb_";
   if ( label != "" ) prefix = "";
   std::string name;
   std::string pname;
   for ( int r = 1 ; r <= njets ; ++r )
   {
      name  = Form("pt_jet%d",r);  pname = Form("%s%s",prefix.c_str(),name.c_str());
      h1_[pname] = std::make_shared<TH1F>(name.c_str(),name.c_str(), 2000,0,2000);
      name  = Form("eta_jet%d",r);  pname = Form("%s%s",prefix.c_str(),name.c_str());
      h1_[pname] = std::make_shared<TH1F>(name.c_str(),name.c_str(), 100,-2.5,2.5);
      name  = Form("btag_jet%d",r);  pname = Form("%s%s",prefix.c_str(),name.c_str());
      h1_[pname] = std::make_shared<TH1F>(name.c_str(),name.c_str(), 1000,0,1.);
      for ( int r2 = r+1 ; r2 <= njets ; ++r2 )
      {
         if ( r2 == r ) continue;
         name  = Form("m_jet%d%d",r,r2);  pname = Form("%s%s",prefix.c_str(),name.c_str());
         h1_[pname] = std::make_shared<TH1F>(name.c_str(),name.c_str(), 3000,0,3000);
         
      } 
   }
   name  = "mbb";  pname = Form("%s%s",prefix.c_str(),name.c_str());
   h1_[pname] = std::make_shared<TH1F>(name.c_str(),"MSSM Hbb mbb (deprecated?)", 3000,0,3000);
   
}
