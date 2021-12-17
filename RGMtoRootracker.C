//Look at hipo file and calculate variables from scattered electron
//Output to rootracker format tree
//skeleton version


//#include "DC_Fiducial_Cuts_CLAS12.cxx"
#include <iostream>
#include <cstdlib>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <vector>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TPaletteAxis.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include <stdlib.h>
#include "Riostream.h"
#include "TLine.h"
#include "TVirtualPad.h"
#include "TClass.h"
#include "TVirtualX.h"
#include "TMath.h"
#include "TStyle.h"
#include "Constants.h"

using namespace clas12;

 


void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
  rp->par()->getPz(),p4.M());

}


void RGMtoRootracker(TString inFileName){

  TString inputFile = inFileName;

  // Creating a TChain of all the input files
  TChain fake("hipo");
  // Adding the different input files to the TChain
  fake.Add(inputFile.Data());
  // fake.Add(inputFile2.Data());


  auto db=TDatabasePDG::Instance();

  double en_beam = 10.6;  //nominal RGA value, RGM will run at 1, 2, 4 and 6 GeV
  //en_beam["1161"]=1.161;
  //en_beam["2261"]=2.261;
  //en_beam["4461"]=4.461;
  //en_beam["6661"]=6.661;
  //try and figure out how to read these in from epics variables

  TLorentzVector beam(0,0,en_beam,en_beam); // beam Px,Py,Pz,E  //CURRENTLY SET TO RGA NOMINAL VALUES
  TLorentzVector target(0,0,0,db->GetParticle(2212)->Mass()); // target Px,Py,Pz,E
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass()); // scattered e^- Px,Py,Pz,E
  TLorentzVector pr(0,0,0,db->GetParticle(2212)->Mass()); // proton Px,Py,Pz,E
  TLorentzVector pip(0,0,0,db->GetParticle(211)->Mass()); // pi^+ Px,Py,Pz,E
  TLorentzVector pim(0,0,0,db->GetParticle(-211)->Mass()); // pi^- Px,Py,Pz,E
  TLorentzVector kpl(0,0,0,db->GetParticle(321)->Mass()); //k^+ Px,Py,Pz,E
  TLorentzVector kmi(0,0,0,db->GetParticle(-321)->Mass()); //k^- Px,Py,Pz,E
  //pi0's, others?

  std::vector<TLorentzVector> loc_proton_v4;
  //std::vector<TLorentzVector> loc_pion_v4;
  std::vector<TLorentzVector> loc_pipl_v4;
  std::vector<TLorentzVector> loc_pimi_v4;
  std::vector<TLorentzVector> loc_kpl_v4;
  std::vector<TLorentzVector> loc_kmi_v4;
  
  
  
  Int_t negative, positive, nonelectron, nonproton, nonpion, npip, npim, nprot, nkpl, nkmi;
  
  
  // Retrieving list of files
  auto files=fake.GetListOfFiles();
  // Gets total events in all files for run dependence binning
  Int_t Bins = files->GetEntries();
  // Output file location and name
  
  
  TString outFileName( inFileName(0,inFileName.Length()-5) + "_rootracker_filtered.root"); //trim '.hipo' and add '.root' for output file name
  
  
  
  TFile fileOutput1(outFileName,"recreate");
  
  
  
  int TargetPdgCode, TargetZ, TargetA;

  //RGA tests, hydrogen target  
  TargetZ = 1;
  TargetA = 1;

  //   if (ftarget=="3He") { TargetPdgCode = 1000020030; TargetZ = 2; TargetA = 3; }
  //   if (ftarget=="4He") { TargetPdgCode = 1000020040; TargetZ = 2; TargetA = 4; }
  //   if (ftarget=="C12") { TargetPdgCode = 1000060120; TargetZ = 6; TargetA = 12; }
  //   if (ftarget=="56Fe") { TargetPdgCode = 1000260560; TargetZ = 26; TargetA = 56; }
  
  //setup and declaration of genie tree
  TTree* mytree;

  //--------
  //Tree variables declared here


  //--------

  //--------
  //Tree set up here
  mytree = new TTree("rootracker","rootracker");

  //Branches declared here
  //mytree->Branch("iev", &genie_iev, "iev/I");  //example from RGMtoGenie

  //--------

  // Looping over data
  
  //Loop over files
  for(Int_t i=0;i<files->GetEntries();i++){
    
    //create the event reader
    clas12reader c12(files->At(i)->GetTitle());
    
    //Binno++; // Count the number of files, therefore the number of x bins
    
    c12.setEntries(1E5);
    while(c12.next()==true){
      
      //Initialise variables

      
      int NEventsTotal = 0;
      

      
      
      
      
      auto particles = c12.getDetParticles();
      
      negative = 0;
      positive = 0;
      nonelectron = 0;
      nonproton = 0;
      nonpion = 0;
      nprot = 0;
      npip = 0;
      npim = 0;      
      nkpl = 0;
      nkmi = 0;      


      auto electrons=c12.getByID(11);
      auto protons=c12.getByID(2212);
      auto pips=c12.getByID(211);
      auto pims=c12.getByID(-211);
      auto kpls=c12.getByID(321);
      auto kmis=c12.getByID(-321);
      
      
      //----------------------------
      //-----Electron Selection-----
      //----------------------------
      if(electrons.size()<1) continue; //skip if no electron
      
      SetLorentzVector(el,electrons[0]);
      
      //----------------------
      //we can add some basic electron selection (vertex cuts, etc) as we go
      
      //----------------------

       //some kinematic quantities from the scattered electron (beam energy must be appropriately defined)
       double nu = -(el-beam).E();
       double Q2 = -(el-beam).Mag2();
       double x_bjk = Q2/(2*m_prot*nu);
       TVector3 V3_q = (beam-el).Vect();
       double W_var = TMath::Sqrt((m_prot+nu)*(m_prot+nu)-V3_q*V3_q);
       
//        genie_x = x_bjk;
//        genie_Q2reco = Q2;
//        genie_W = W_var;
//        genie_El = el.E();
//        genie_pxl = el.Px();
//        genie_pyl = el.Py();
//        genie_pzl = el.Pz();
//        genie_pl = el.Rho();
//        genie_cthl = cos(el.Theta());
//        genie_vtxx = 0;
//        genie_vtxy = 0;
//        genie_vtxz = 0;


      //-----------------------------
      //--End of Electron Selection--
      //-----------------------------
       
       //---------------------------------------------
       //---------START OF HADRON SELECTION-----------
       //---------------------------------------------
       loc_proton_v4.clear();
       //loc_pion_v4.clear();
       loc_pipl_v4.clear();
       loc_pimi_v4.clear();
       
       int IndexInFinalStateParticleArray = 0;
       
       for(auto& p : particles){
	 
	 if(p->par()->getCharge() > 0 ){
	   positive++;
	   
	   if(p->par()->getPid()==2212){
	     //pr.SetXYZM(p->par()->getPx(),p->par()->getPy(),p->par()->getPz(),db->GetParticle(2212)->Mass());
	     SetLorentzVector(pr,protons[nprot]);

//  	     genie_pdgf[IndexInFinalStateParticleArray] = 2212;
//  	     genie_Ef[IndexInFinalStateParticleArray] = pr.E();
//  	     genie_pxf[IndexInFinalStateParticleArray] = pr.Px();
//  	     genie_pyf[IndexInFinalStateParticleArray] = pr.Py();
//  	     genie_pzf[IndexInFinalStateParticleArray] = pr.Pz();
//  	     genie_pf[IndexInFinalStateParticleArray] = pr.Rho();
//  	     genie_cthf[IndexInFinalStateParticleArray] = pr.CosTheta();
	     
 	     //double CorrectedProtonMomentum = ProtonMomCorrection_He3_4Cell(ftarget,V4_uncorrprot,p_vert_corr);
 	     //V4_uncorrprot.SetRho(CorrectedProtonMomentum);
 	     //double CorrectedProtonE = sqrt(CorrectedProtonMomentum*CorrectedProtonMomentum + m_prot*m_prot);
 	     //V4_uncorrprot.SetE(CorrectedProtonE);
	     
//  	     genie_pdgf[IndexInFinalStateParticleArray+shift] = 2212;
//  	     genie_Ef[IndexInFinalStateParticleArray+shift] = pr.E();
//  	     genie_pxf[IndexInFinalStateParticleArray+shift] = pr.Px();
//  	     genie_pyf[IndexInFinalStateParticleArray+shift] = pr.Py();
//  	     genie_pzf[IndexInFinalStateParticleArray+shift] = pr.Pz();
//  	     genie_pf[IndexInFinalStateParticleArray+shift] = pr.Rho();
//  	     genie_cthf[IndexInFinalStateParticleArray+shift] = pr.CosTheta();
//  	     IndexInFinalStateParticleArray++;
	     

	     nprot++;
	     loc_proton_v4.push_back(pr);
	   }
	   
 	   if(p->par()->getPid()==211){
 	     SetLorentzVector(pip,pips[npip]);

 	     double PiPlusP = pip.Vect().Mag();
 	     double PiPlusE = TMath::Sqrt(PiPlusP*PiPlusP + m_pipl * m_pipl);
	     
//  	     genie_pdgf[IndexInFinalStateParticleArray] = 211;
//  	     genie_Ef[IndexInFinalStateParticleArray] = PiPlusE;
//  	     genie_pxf[IndexInFinalStateParticleArray] = pip.Vect().X();
//  	     genie_pyf[IndexInFinalStateParticleArray] = pip.Vect().Y();
//  	     genie_pzf[IndexInFinalStateParticleArray] = pip.Vect().Z();
//  	     genie_pf[IndexInFinalStateParticleArray] = PiPlusP;
//  	     genie_cthf[IndexInFinalStateParticleArray] = pip.Vect().CosTheta();
//  	     IndexInFinalStateParticleArray++;
	     
 	     npip++;
 	     loc_pipl_v4.push_back(pip);
 	   }

	   if(p->par()->getPid()==321){
	     SetLorentzVector(kpl,kpls[nkpl]);

//  	     genie_pdgf[IndexInFinalStateParticleArray] = 321;
//  	     genie_Ef[IndexInFinalStateParticleArray] = kpl.E();
//  	     genie_pxf[IndexInFinalStateParticleArray] = kpl.Px();
//  	     genie_pyf[IndexInFinalStateParticleArray] = kpl.Py();
//  	     genie_pzf[IndexInFinalStateParticleArray] = kpl.Pz();
//  	     genie_pf[IndexInFinalStateParticleArray] = kpl.Rho();
//  	     genie_cthf[IndexInFinalStateParticleArray] = kpl.CosTheta();
//  	     IndexInFinalStateParticleArray++;
	     

	     nkpl++;
	     loc_kpl_v4.push_back(kpl);
	   }
	 }
	 
 	 else if(p->par()->getCharge() < 0 ){
 	   negative++;
	   
 	   if(p->par()->getPid()==-211){
 	     SetLorentzVector(pim,pims[npim]);

 	     double PiMinusP = pim.Vect().Mag();
 	     double PiMinusE = TMath::Sqrt(PiMinusP*PiMinusP + m_pimi * m_pimi);
	     
//  	     genie_pdgf[IndexInFinalStateParticleArray] = -211;
//  	     genie_Ef[IndexInFinalStateParticleArray] = PiMinusE;
//  	     genie_pxf[IndexInFinalStateParticleArray] = pim.Vect().X();
//  	     genie_pyf[IndexInFinalStateParticleArray] = pim.Vect().Y();
//  	     genie_pzf[IndexInFinalStateParticleArray] = pim.Vect().Z();
//  	     genie_pf[IndexInFinalStateParticleArray] = PiMinusP;
//  	     genie_cthf[IndexInFinalStateParticleArray] = pim.Vect().CosTheta();
//  	     IndexInFinalStateParticleArray++;


 	     npim++;
 	     loc_pimi_v4.push_back(pim);
 	   }

	   if(p->par()->getPid()==-321){
	     SetLorentzVector(kmi,kmis[nkmi]);

//  	     genie_pdgf[IndexInFinalStateParticleArray] = -321;
//  	     genie_Ef[IndexInFinalStateParticleArray] = kmi.E();
//  	     genie_pxf[IndexInFinalStateParticleArray] = kmi.Px();
//  	     genie_pyf[IndexInFinalStateParticleArray] = kmi.Py();
//  	     genie_pzf[IndexInFinalStateParticleArray] = kmi.Pz();
//  	     genie_pf[IndexInFinalStateParticleArray] = kmi.Rho();
//  	     genie_cthf[IndexInFinalStateParticleArray] = kmi.CosTheta();
//  	     IndexInFinalStateParticleArray++;
	     

	     nkmi++;
	     loc_kmi_v4.push_back(kpl);
	   }
	 }
	 
	 
	 //what about neutrals???
	 
       }


    //---------------------------------------------
    //----------END OF HADRON SELECTION------------
    //---------------------------------------------

//        genie_nfp = nprot;
//        genie_nfn = 0;
//        genie_nfpip = npip;
//        genie_nfpim = npim;

//        genie_nfkp = nkpl;
//        genie_nfkm = nkmi;

//        genie_nfpi0 = 0;
//        //genie_nfpi0 = ec_num_n;

//        //genie_nf = genie_nfp + genie_nfn + genie_nfpip + genie_nfpim + genie_nfpi0;
//        genie_nf = genie_nfp + genie_nfn + genie_nfpip + genie_nfpim + genie_nfkp + genie_nfkm + genie_nfpi0;
       
//        //genie_iev = NEventsTotal;
//        //NEventsTotal++;


       mytree->Fill();

      }
      
    //}
    


//       //v_Runno.push_back(to_string(runno)); // Converting runno integer to a string
//   //}

//mytree->Write();

  //saving the file
  fileOutput1.Write();
  fileOutput1.Close();
  }
}


//Function wrapper for hard coded filename
void RGMtoRootracker(){

  // Data files to process
  //TString inFile("/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v0/dst/train/skim4/skim4_005*.hipo");
  //TString inFile("/volatile/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v0/dst/train/dst/train/skim4/*.hipo");
  TString inFile("/home/stuart/CLAS/Data/skim4_00503*.hipo");
  // TString inputFile2("/volatile/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train/inc/*.hipo");

  RGMtoRootracker(inFile); //call the analysis function with this filename 

}

