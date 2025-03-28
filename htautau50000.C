//#include "/home/nfyst/anaconda3/include/TLorentzVector.h"
// include-filer som trengs i "htautau50000.h"
#include <cmath>
//#include <Math/GenVector/Boost.h>
//#include <Math/GenVector/LorentzRotation.h>
//#include <Math/GenVector/LorentzVector.h>
//#include <TLorentzRotation.h>
#include <TLorentzVector.h>
#include <TVector.h>
#include <TRef.h>
#include <TRefArray.h>



#define htautau50000_cxx
#include "htautau50000.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


void htautau50000::Loop()
{
//   In a ROOT session, you can do:
//      root> .L htautau50000.C
//      root> htautau50000 t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
   
     TLorentzVector mu1,mu2,ele1,ele2,partmu1,partmu2,jet1,jet2,pi1,pi2;
     TLorentzVector higgs4P,higgs4Pcm;
     TVector3  thbeta;

     cout << "Booking " << endl;
     TH1 *histNJet = new TH1F("n_Jet","Number of jets",10,-.5,9.5);
     TH1 *histJetPT = new TH1F("jet_pt", "jet P_{T}", 100, 0.0, 100.0);
     TH1 *histJetEta = new TH1F("jet_eta", "jet eta", 100, -5.,5.);
     TH1 *histtauJetPT = new TH1F("taujet_pt", "taujet P_{T}", 100, 0.0, 100.0);
     TH1 *histtauJetEta = new TH1F("taujet_eta", "taujet eta", 100, -5.,5.);
     TH1 *histpipiM = new TH1F("pipi_M","pion pair mass",75,0.,150.);
     TH1 *histcharge = new TH1F("charge","Taujet charge",10,-4.5,5.5);
     
     TH1 *histNGenJet = new TH1F("n_GenJet","Number of generated jets",10,-.5,9.5);
     TH1 *histGenJetPT = new TH1F("genjet_pt", "jet P_{T}", 100, 0.,100.);     
     TH1 *histGenJetEta = new TH1F("genjet_eta", "jet eta", 100, -5.,5.);
     TH1 *histtauGenJetPT = new TH1F("taugenjet_pt", "taujet P_{T}", 100, 0.0, 100.0);
     TH1 *histtauGenJetEta = new TH1F("taugenjet_eta", "taujet eta", 100, -5.,5.);    

     TH1 *histNMU = new TH1F("n_muons","Number of muons",10,-.5,9.5);    
     TH1 *histMuonEta = new TH1F("muon_eta", "muon eta", 100, -5.,5.);
     TH1 *histMuonPT = new TH1F("muon_pt", "muon P_{T}", 100, 0.0, 100.0);
     TH1 *histMumuM = new TH1F("Mumu_M","muon pair mass",75,0.,150.);

     TH1 *histPIDPART = new TH1F("pidPart","Particle ID",220,-.5,219.5);
     TH1 *histNgenMU = new TH1F("n_gen_muons","Generated Number of muons",10,-.5,9.5);   
     TH1 *histNgenPART = new TH1F("n_gen_PART","Generated Number of Particles",100,-.5,999.5); 
     TH1 *histGenMuonEta = new TH1F("muon_gen_eta", "Generated muon eta", 100, -5.,5.);
     TH1 *histGenMuonPT = new TH1F("muon_gen__pt", "Generated muon P_{T}", 100, 0.0, 100.0);
     TH1 *histGenMumuM = new TH1F("Mumu_gen_M","Generated muon pair mass",75,0.,150.);  
     TH1 *histbetaHiggs = new TH1F("betaHiggs","Beta of higgs",50,0.,1.); 
     TH1 *histHiggsPcm = new TH1F("higgsPcm","Higgs Pcm",50,0.,150.);
     
     TH2 *histGenPmu1vsPmu2 = new TH2F("GenPmu1vsPmu2","Generated pmu1 vs pmu2",50,0.,100,50,0.,100);     
     TH2 *histJetpi1vspi2 = new TH2F("Jetpi1vspi2","Jetpi1 vs Jetpi2",50,0.,100,50,0.,100);
     cout << "finished booking " << endl;
     TFile myOutputFile("hist_htautau200000.root","RECREATE");
cout << "finished outputfile " << endl;

   Long64_t nentries = fChain->GetEntriesFast();
   // speed!
   fChain->SetBranchStatus("*",0);
   fChain->SetBranchStatus("Muon",1);
   fChain->SetBranchStatus("Muon.PT",1);  
   fChain->SetBranchStatus("Muon.Eta",1);
   fChain->SetBranchStatus("Muon.Phi",1);
   fChain->SetBranchStatus("Jet",1);
   fChain->SetBranchStatus("Jet.PT",1);
   fChain->SetBranchStatus("Jet.Eta",1);
   fChain->SetBranchStatus("Jet.Phi",1);
   fChain->SetBranchStatus("Jet.TauTag",1);
   fChain->SetBranchStatus("Jet.Charge",1);
   // What is this?
   fChain->SetBranchStatus("GenJet",1);
   fChain->SetBranchStatus("GenJet.PT",1);
   fChain->SetBranchStatus("GenJet.Eta",1);
   fChain->SetBranchStatus("GenJet.TauTag",1);
   
   // Generator level info
   fChain->SetBranchStatus("Particle",1);
   fChain->SetBranchStatus("Particle.PT",1);  
   fChain->SetBranchStatus("Particle.Eta",1);
   fChain->SetBranchStatus("Particle.Phi",1);
   fChain->SetBranchStatus("Particle.PID",1);
   fChain->SetBranchStatus("Particle.Mass",1);
   cout << "finished branches " << endl;
   
   Long64_t nbytes = 0, nb = 0;   //nentries
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
       // se etter generert informasjon
      int ngenMU = 0; //or tau,e, or pi   (11,13,15,211) 
      int isel1 = 0;
      int isel2 = 0;
      Double_t high1 = -1.;
      Double_t high2 = -2.;
      histNgenPART->Fill(Particle_); 
          
      // generator info
      for (int i=0; i<Particle_;i++) {
         histPIDPART->Fill(Particle_PID[i]);
         // Higgs is particle 25
         if((Particle_PID[i] == 25) && (Particle_PT[i] > .001)) {  
            if (jentry < 10) {
            cout << " pt "<< Particle_PT[i] <<" eta " << Particle_Eta[i] << " phi " << Particle_Phi[i] << " M  " << Particle_Mass[i] << endl;}      
            higgs4P.SetPtEtaPhiM(Particle_PT[i],Particle_Eta[i],Particle_Phi[i],Particle_Mass[i]);
            //PtEtaPhiMVector higgs4P(Particle_PT[i],Particle_Eta[i],Particle_Phi[i],Particle_Mass[i]);
            Double_t hP = higgs4P.P();
            Double_t hE = higgs4P.E();
            Double_t hBetaAbs = hP/hE;
            histbetaHiggs->Fill(hBetaAbs);
            
            
            thbeta = higgs4P.BoostVector();
            //TVector3 vhbeta(hbeta[0],hbeta[1],hbeta[2]);
        
            higgs4Pcm = higgs4P;
            higgs4Pcm.Boost(-thbeta); 
           
          
            if (jentry < 10){
            cout << " hP " << hP << " hE " << hE << " higgs4P.Px " << higgs4P.Px() << endl;
            cout << " HBETA " << thbeta.X() << "  " << thbeta.Y() << " " << thbeta.Z() << endl;
            cout << " boostapx "<< higgs4Pcm.Px() <<" mass " << higgs4Pcm.M() << endl;
           
            }
            
            
            histHiggsPcm->Fill(higgs4Pcm.P());
            
            
         }               
         
         if(Particle_PID[i] == 13) {
            
            ngenMU++;
            histGenMuonPT->Fill(Particle_PT[i]);
            histGenMuonEta->Fill(Particle_Eta[i]);
            
            if (Particle_PT[i] > high1){                                 
               isel1 = i;
               high1 = Particle_PT[i];
            }
                     
         }   
         if(Particle_PID[i] == -13) {
            
            ngenMU++;
            histGenMuonPT->Fill(Particle_PT[i]);
            histGenMuonEta->Fill(Particle_Eta[i]);
            
            if (Particle_PT[i] > high2){                                 
               isel2 = i;
               high2 = Particle_PT[i];
            }
                     
         }  
      }
      if (ngenMU >1){
      
        partmu1.SetPtEtaPhiM(Particle_PT[isel1],Particle_Eta[isel1],Particle_Phi[isel1],0.105);
        partmu2.SetPtEtaPhiM(Particle_PT[isel2],Particle_Eta[isel2],Particle_Phi[isel2],0.105);
            
        Double_t Partminv = (partmu1+partmu2).M();       
        histGenMumuM->Fill(Partminv);
        histNgenMU->Fill(ngenMU);   
        histGenPmu1vsPmu2->Fill(partmu1.E(),partmu2.E());
      }
      
      histNMU->Fill(Muon_);
      mu1.SetPtEtaPhiM(0.,0.,0.,.105);
      mu2.SetPtEtaPhiM(0.,0.,0.,.105);
      partmu1.SetPtEtaPhiM(0.,0.,0.,.105);
      partmu2.SetPtEtaPhiM(0.,0.,0.,.105);
      for (int i=0; i<Muon_;i++) {
         histMuonPT->Fill(Muon_PT[i]);
         histMuonEta->Fill(Muon_Eta[i]);
      }
      if (Muon_ == 2) {
         mu1.SetPtEtaPhiM(Muon_PT[0],Muon_Eta[0],Muon_Phi[0],0.105);
         mu2.SetPtEtaPhiM(Muon_PT[1],Muon_Eta[1],Muon_Phi[1],0.105);
         Double_t minv = (mu1+mu2).M();
         histMumuM->Fill(minv);
      }
      histNJet->Fill(Jet_);
      int itau1 = 0;
      int itau2 = 0;
      Double_t tauE1 = -2.;
      Double_t tauE2 = -2.;
      for (int i=0; i<Jet_;i++) {
         histJetPT->Fill(Jet_PT[i]);
         histJetEta->Fill(Jet_Eta[i]);
         if (Jet_TauTag[i] == 1){
            histcharge->Fill(Jet_Charge[i]);
            histtauJetPT->Fill(Jet_PT[i]);
            histtauJetEta->Fill(Jet_Eta[i]);
            if (Jet_Charge[i] > 0) {
               if (Jet_PT[i] > tauE1) {
                  itau1 = i;
                  tauE1 = Jet_PT[i];
               }
            }
            if (Jet_Charge[i] < 0) {
               if (Jet_PT[i] > tauE2) {
                  itau2 = i;
                  tauE2 = Jet_PT[i];
               }
            }
         }
      }
      if (Jet_Charge[itau1]*Jet_Charge[itau2] <0){
         pi1.SetPtEtaPhiM(Jet_PT[itau1],Jet_Eta[itau1],Jet_Phi[itau1],0.135);
         pi2.SetPtEtaPhiM(Jet_PT[itau2],Jet_Eta[itau2],Jet_Phi[itau2],0.135);
         Double_t minvpi = (pi1+pi2).M();
         histpipiM->Fill(minvpi);
         histJetpi1vspi2->Fill(Jet_PT[itau1],Jet_PT[itau2]);
      }        
      histNGenJet->Fill(GenJet_);
      for (int i=0; i<GenJet_;i++) {
         histGenJetPT->Fill(GenJet_PT[i]);
         histGenJetEta->Fill(GenJet_Eta[i]);
         if (GenJet_TauTag[i] == 1){
            histtauGenJetPT->Fill(GenJet_PT[i]);
            histtauGenJetEta->Fill(GenJet_Eta[i]);
         }
      }    
     
      
      
      
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
   cout << "eventloop finished" << endl;
   histNJet->Write();
   
   histJetPT->Write();
   histJetEta->Write();
   histtauJetPT->Write();
   histtauJetEta->Write();
   histpipiM->Write();
   histJetpi1vspi2->Write();
   histcharge->Write();
   
   histNGenJet->Write();
   histGenJetPT->Write();
   histGenJetEta->Write();
   histtauGenJetPT->Write();
   histtauGenJetEta->Write();  

   histMuonEta->Write();
   histMuonPT->Write();
   histNMU->Write();
   histMumuM->Write();
   histNgenPART->Write();
   histGenMuonPT->Write();
   histGenMuonEta->Write();
   histGenMumuM->Write();
   histNgenMU->Write();
   histGenPmu1vsPmu2->Write();
   histPIDPART->Write();
   histbetaHiggs->Write();
   histHiggsPcm->Write();
   
   
   myOutputFile.Close();
}
