//Look at hipo file and calculate variables from scattered electron
//Output to GENIE format tree
//No hadrons yet


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


//void RGMtoGenie(TString inFileName){
void RGMtoGenie(TString inFileName, TString outFileName, const std::string databaseF){

  TString inputFile = inFileName;

  // Creating a TChain of all the input files
  TChain fake("hipo");
  // Adding the different input files to the TChain
  fake.Add(inputFile.Data());
  // fake.Add(inputFile2.Data());


  auto db=TDatabasePDG::Instance();

  //open comunications with RCDB DBfile must be name of prepared database root file
  const std::string DBfile = databaseF;
  clas12databases::SetRCDBRootConnection(DBfile);
  
  double en_beam = 6.0;  //This needs to come from RCDB
  //en_beam["1161"]=1.161;
  //en_beam["2261"]=2.261;
  //en_beam["4461"]=4.461;
  //en_beam["6661"]=6.661;
  //try and figure out how to read these in from epics variables
  Double_t beam_E;

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
  
  
  ////TString outFileName( inFileName(0,inFileName.Length()-5) + "_genie_filtered.root"); //trim '.hipo' and add '.root' for output file name
  //TString outFileName( "ca40_qe_genie_test1_torus-1.0_sol-1.0_genie_filtered.root"); //test output
  
  
  
  TFile fileOutput1(outFileName,"recreate");
  
  
  
  int TargetPdgCode, TargetZ, TargetA;

  //RGM tests, ca-40 target  
  TargetZ = 20;
  TargetA = 20;

  //   if (ftarget=="3He") { TargetPdgCode = 1000020030; TargetZ = 2; TargetA = 3; }
  //   if (ftarget=="4He") { TargetPdgCode = 1000020040; TargetZ = 2; TargetA = 4; }
  //   if (ftarget=="C12") { TargetPdgCode = 1000060120; TargetZ = 6; TargetA = 12; }
  //   if (ftarget=="56Fe") { TargetPdgCode = 1000260560; TargetZ = 26; TargetA = 56; }
  
  //setup and declaration of genie tree
  TTree* mytree;
  
  Int_t           genie_iev;
  Int_t           genie_neu;
  Int_t           genie_fspl;
  Int_t           genie_tgt;
  Int_t           genie_Z;
  Int_t           genie_A;
  Int_t           genie_hitnuc;
  Int_t           genie_hitqrk;
  Int_t           genie_resid;
  Bool_t          genie_sea;
  Bool_t          genie_qel;
  Bool_t          genie_mec;
  Bool_t          genie_res;
  Bool_t          genie_dis;
  Bool_t          genie_coh;
  Bool_t          genie_dfr;
  Bool_t          genie_imd;
  Bool_t          genie_imdanh;
  Bool_t          genie_singlek;
  Bool_t          genie_nuel;
  Bool_t          genie_em;
  Bool_t          genie_CC;
  Bool_t          genie_nc;
  Bool_t          genie_charm;
  Int_t           genie_neut_code;
  Int_t           genie_nuance_code;
  Double_t        genie_wght;
  Double_t        genie_xs;
  Double_t        genie_ys;
  Double_t        genie_ts;
  Double_t        genie_Q2s;
  Double_t        genie_Ws;
  Double_t        genie_x;
  Double_t        genie_y;
  Double_t        genie_t;
  Double_t        genie_Q2reco;
  Double_t        genie_W;
  Double_t        genie_EvRF;
  Double_t        genie_Ev;
  Double_t        genie_pxv;
  Double_t        genie_pyv;
  Double_t        genie_pzv;
  Double_t        genie_En;
  Double_t        genie_pxn;
  Double_t        genie_pyn;
  Double_t        genie_pzn;
  Double_t        genie_El;
  Double_t        genie_pxl;
  Double_t        genie_pyl;
  Double_t        genie_pzl;
  Double_t        genie_pl;
  Double_t        genie_cthl;
  Int_t           genie_nfp;
  Int_t           genie_nfn;
  Int_t           genie_nfpip;
  Int_t           genie_nfpim;
  Int_t           genie_nfpi0;
  Int_t           genie_nfkp;
  Int_t           genie_nfkm;
  Int_t           genie_nfk0;
  Int_t           genie_nfem;
  Int_t           genie_nfother;
  Int_t           genie_nip;
  Int_t           genie_nin;
  Int_t           genie_nipip;
  Int_t           genie_nipim;
  Int_t           genie_nipi0;
  Int_t           genie_nikp;
  Int_t           genie_nikm;
  Int_t           genie_nik0;
  Int_t           genie_niem;
  Int_t           genie_niother;
  Int_t           genie_ni;
  int             InitialStateParticles = 2;
  Int_t           genie_pdgi[InitialStateParticles];   //[ni]
  Int_t           genie_resc[InitialStateParticles];   //[ni]
  Double_t        genie_Ei[InitialStateParticles];   //[ni]
  Double_t        genie_pxi[InitialStateParticles];   //[ni]
  Double_t        genie_pyi[InitialStateParticles];   //[ni]
  Double_t        genie_pzi[InitialStateParticles];   //[ni]
  Int_t           genie_nf;
  const int FinalStateParticles = 120;
  Int_t           genie_pdgf[FinalStateParticles];   //[nf]
  Double_t        genie_Ef[FinalStateParticles];   //[nf]
  Double_t        genie_pxf[FinalStateParticles];   //[nf]
  Double_t        genie_pyf[FinalStateParticles];   //[nf]
  Double_t        genie_pzf[FinalStateParticles];   //[nf]
  Double_t        genie_pf[FinalStateParticles];   //[nf]
  Double_t        genie_cthf[FinalStateParticles];   //[nf]
  Double_t        genie_vtxx;
  Double_t        genie_vtxy;
  Double_t        genie_vtxz;
  Double_t        genie_vtxt;
  Double_t        genie_sumKEf;
  Double_t        genie_calresp0;




  mytree = new TTree("gst","gst");  
  
  mytree->Branch("iev", &genie_iev, "iev/I");
  mytree->Branch("neu", &genie_neu, "neu/I");
  mytree->Branch("fspl", &genie_fspl, "fspl/I");
  mytree->Branch("tgt", &genie_tgt, "tgt/I");
  mytree->Branch("Z", &genie_Z, "Z/I");
  mytree->Branch("A", &genie_A, "A/I");
  mytree->Branch("hitnuc", &genie_hitnuc, "hitnuc/I");
  mytree->Branch("hitqrk", &genie_hitqrk, "hitqrk/I");
  mytree->Branch("resid", &genie_resid, "resid/I");
  mytree->Branch("sea", &genie_sea, "sea/O");
  mytree->Branch("qel", &genie_qel, "qel/O");
  mytree->Branch("mec", &genie_mec, "mec/O");
  mytree->Branch("res", &genie_res, "res/O");
  mytree->Branch("dis", &genie_dis, "dis/O");
  mytree->Branch("coh", &genie_coh, "coh/O");
  mytree->Branch("dfr", &genie_dfr, "dfr/O");
  mytree->Branch("imd", &genie_imd, "imd/O");
  mytree->Branch("imdanh", &genie_imdanh, "imdanh/O");
  mytree->Branch("singlek", &genie_singlek, "singlek/O");
  mytree->Branch("nuel", &genie_nuel, "nuel/O");
  mytree->Branch("em", &genie_em, "em/O");
  mytree->Branch("cc", &genie_CC, "cc/O");
  mytree->Branch("nc", &genie_nc, "nc/O");
  mytree->Branch("charm", &genie_charm, "charm/O");
  mytree->Branch("neut_code", &genie_neut_code, "neut_code/I");
  mytree->Branch("nuance_code", &genie_nuance_code, "nuance_code/I");
  mytree->Branch("wght", &genie_wght, "wght/D");
  mytree->Branch("xs", &genie_xs, "xs/D");
  mytree->Branch("ys", &genie_ys, "ys/D");
  mytree->Branch("ts", &genie_ts, "ts/D");
  mytree->Branch("Q2s", &genie_Q2s, "Q2s/D");
  mytree->Branch("Ws", &genie_Ws, "Ws/D");
  mytree->Branch("x", &genie_x, "x/D");
  mytree->Branch("y", &genie_y, "x/D");
  mytree->Branch("t", &genie_t, "t/D");
  mytree->Branch("Q2", &genie_Q2reco, "Q2/D");
  mytree->Branch("W", &genie_W, "W/D");
  mytree->Branch("EvRF", &genie_EvRF, "EvRF/D");
  mytree->Branch("Ev", &genie_Ev, "Ev/D");
  mytree->Branch("pxv", &genie_pxv, "pxv/D");
  mytree->Branch("pyv", &genie_pyv, "pyv/D");
  mytree->Branch("pzv", &genie_pzv, "pzv/D");
  mytree->Branch("En", &genie_En, "En/D");
  mytree->Branch("pxn", &genie_pxn, "pxn/D");
  mytree->Branch("pyn", &genie_pyn, "pyn/D");
  mytree->Branch("pzn", &genie_pzn, "pzn/D");
  mytree->Branch("El", &genie_El, "El/D");
  mytree->Branch("pxl", &genie_pxl, "pxl/D");
  mytree->Branch("pyl", &genie_pyl, "pyl/D");
  mytree->Branch("pzl", &genie_pzl, "pzl/D");
  mytree->Branch("pl", &genie_pl, "pl/D");
  mytree->Branch("cthl", &genie_cthl, "cthl/D");
  mytree->Branch("nfp", &genie_nfp, "nfp/I");
  mytree->Branch("nfn", &genie_nfn, "nfn/I");
  mytree->Branch("nfpip", &genie_nfpip, "nfpip/I");
  mytree->Branch("nfpim", &genie_nfpim, "nfpim/I");
  mytree->Branch("nfpi0", &genie_nfpi0, "nfpi0/I");
  mytree->Branch("nfkp", &genie_nfkp, "nfkp/I");
  mytree->Branch("nfkm", &genie_nfkm, "nfkm/I");
  mytree->Branch("nfk0", &genie_nfk0, "nfk0/I");
  mytree->Branch("nfem", &genie_nfem, "nfem/I");
  mytree->Branch("nfother", &genie_nfother, "nfother/I");
  mytree->Branch("nip", &genie_nip, "nip/I");
  mytree->Branch("nin", &genie_nin, "nin/I");
  mytree->Branch("nipip", &genie_nipip, "nipip/I");
  mytree->Branch("nipim", &genie_nipim, "nipim/I");
  mytree->Branch("nipi0", &genie_nipi0, "nipi0/I");
  mytree->Branch("nikp", &genie_nikp, "nikp/I");
  mytree->Branch("nikm", &genie_nikm, "nikm/I");
  mytree->Branch("nik0", &genie_nik0, "nik0/I");
  mytree->Branch("niem", &genie_niem, "niem/I");
  mytree->Branch("niother", &genie_niother, "niother/I");
  mytree->Branch("ni", &genie_ni, "ni/I");
  mytree->Branch("pdgi", &genie_pdgi,"pdgi[2]/I");
  mytree->Branch("resc", &genie_resc, "resc[2]/I");
  mytree->Branch("Ei", &genie_Ei, "Ei[2]/D");
  mytree->Branch("pxi", &genie_pxi, "pxi[2]/D");
  mytree->Branch("pyi", &genie_pyi, "pyi[2]/D");
  mytree->Branch("pzi", &genie_pzi, "pzi[2]/D");
  mytree->Branch("nf", &genie_nf, "nf/I");
  
  mytree->Branch("pdgf", &genie_pdgf, "pdgf[120]/I");
  mytree->Branch("Ef", &genie_Ef, "Ef[120]/D");
  mytree->Branch("pxf", &genie_pxf, "pxf[120]/D");
  mytree->Branch("pyf", &genie_pyf, "pyf[120]/D");
  mytree->Branch("pzf", &genie_pzf, "pzf[120]/D");
  mytree->Branch("pf", &genie_pf, "pf[120]/D");
  mytree->Branch("cthf", &genie_cthf, "cthf[120]/D");
  mytree->Branch("vtxx", &genie_vtxx, "vtxx/D");
  mytree->Branch("vtxy", &genie_vtxy, "vtxy/D");
  mytree->Branch("vtxz", &genie_vtxz, "vtxz/D");
  mytree->Branch("vtxt", &genie_vtxt, "vtxt/D");
  mytree->Branch("sumKEf", &genie_sumKEf, "sumKEf/D");
  mytree->Branch("calresp0", &genie_calresp0, "calresp0/D");
  
  int shift = 60;
  


  // Looping over data
  
  //Loop over files
  for(Int_t i=0;i<files->GetEntries();i++){
    
    //create the event reader
    clas12reader c12(files->At(i)->GetTitle());

 //Creating a chain for the data from different filesc12->queryRcdb();
   clas12root::HipoChain chain;
   chain.Add(inFileName);
   chain.SetReaderTags({0});  //create clas12reader with just tag 0 events
   auto config_c12=chain.GetC12Reader();

   //dont talk to QADB
    chain.db()->turnOffQADB();

    //RCDB data is for the current run I think
    auto& rcdbData= config_c12->rcdb()->current();
    
    //Binno++; // Count the number of files, therefore the number of x bins
    
    c12.setEntries(2E5);
    while(c12.next()==true){
      
       //GeV it up
      beam_E = rcdbData.beam_energy*0.001;
      //beam_E = rcdbData.beam_energy;
      
      beam.SetXYZM(0, 0, beam_E, 0.000511);

      en_beam = beam_E;

      //cout << "beam energy " << beam_E << endl;

     
      genie_neu = 11;
      genie_fspl = 11;
      genie_tgt = TargetPdgCode;
      genie_Z = TargetZ;
      genie_A = TargetA;
      genie_hitnuc = -99.;
      genie_hitqrk = -99.;
      genie_resid = -99.;
      genie_sea = false;
      genie_qel = false;
      genie_mec = false;
      genie_res = false;
      genie_dis = false;
      genie_coh = false;
      genie_dfr = false;
      genie_imd = false;
      genie_imdanh = false;
      genie_singlek = false;
      genie_nuel = false;
      genie_em = true;
      genie_CC = false;
      genie_nc = false;
      genie_charm = false;
      genie_neut_code = -99;
      genie_nuance_code = -99;
      genie_wght = 1.;
      genie_xs = -99;
      genie_ys = -99;
      genie_ts = -99;
      genie_Q2s = -99;
      genie_Ws = -99;
      genie_Ev = en_beam;
      genie_EvRF = -99.;
      genie_pxv = 0.;
      genie_pyv = 0.;
      genie_pzv = sqrt(en_beam * en_beam - e_mass*e_mass); 
      genie_En = -99.;
      genie_pxn = -99.;
      genie_pyn = -99.;
      genie_pzn = -99.;
      genie_nfkp = -99;
      genie_nfkm = -99;
      genie_nfk0 = -99;
      genie_nfem = -99;
      genie_nfother = -99;
      genie_nip = -99;
      genie_nin = -99;
      genie_nipip = -99;
      genie_nipim = -99;
      genie_nipi0 = -99;
      genie_nikp = -99;
      genie_nikm = -99;
      genie_nik0 = -99;
      genie_niem = -99;
      genie_niother = -99;
      genie_ni = -99;
      genie_calresp0 = -99;
      genie_y = -99.;
      genie_t = -99;
      
      int NEventsTotal = 0;
      
      for (int WhichInitialStateParticle = 0; WhichInitialStateParticle < InitialStateParticles; WhichInitialStateParticle ++) {
	
	genie_pdgi[WhichInitialStateParticle] = -99;
	genie_resc[WhichInitialStateParticle] = -99;
	genie_Ei[WhichInitialStateParticle] = -99.;
	genie_pxi[WhichInitialStateParticle] = -99.;
	genie_pyi[WhichInitialStateParticle] = -99.;
	genie_pzi[WhichInitialStateParticle] = -99.;
	
      }
      
      for (int WhichFinalStateParticle = 0; WhichFinalStateParticle < FinalStateParticles; WhichFinalStateParticle ++) {
	
	genie_pdgf[WhichFinalStateParticle] = -99;
	genie_Ef[WhichFinalStateParticle] = -99.;
	genie_pxf[WhichFinalStateParticle] = -99.;
	genie_pyf[WhichFinalStateParticle] = -99.;
	genie_pzf[WhichFinalStateParticle] = -99.;
	genie_pf[WhichFinalStateParticle] = -99.;
	genie_cthf[WhichFinalStateParticle] = -99.;
	
      }
      
      
      
      
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
      
      //copy-pasted from CTOF, can be made more efficient, all we're after is the particle 4-vectors      
      //       for(auto& p : particles){
      
      //         // Looking at negative particles
      //         if(p->par()->getCharge() < 0){
      // 	  negative++;
      //           // negative particles not electron are set to pi^-
      //           if(p->par()->getPid() != 11){
      //             nonelectron++;
      //             pim.SetXYZM(p->par()->getPx(),p->par()->getPy(),p->par()->getPz(),db->GetParticle(-211)->Mass());
      //           }
      //         }
      //         // Looking at positive particles
      //         else if(p->par()->getCharge() > 0){
      // 	  positive++;
      // 	  //positive particles not proton are set to pi+
      //           if(p->par()->getPid() != 2212){
      //             nonproton++;
      //             pip.SetXYZM(p->par()->getPx(),p->par()->getPy(),p->par()->getPz(),db->GetParticle(211)->Mass());
      //           }
      // 	  //positive particles not pion are set to proton
      //           if(p->par()->getPid() != 211){
      //             nonpion++;
      //             pr.SetXYZM(p->par()->getPx(),p->par()->getPy(),p->par()->getPz(),db->GetParticle(2212)->Mass());
      //           }
      // 	}
      //       }
      
      //       //if(nonpion == 1 && electrons.size() == 1 && pims.size() == 1 && pips.size() == 1){
      
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
       
       genie_x = x_bjk;
       genie_Q2reco = Q2;
       genie_W = W_var;
       genie_El = el.E();
       genie_pxl = el.Px();
       genie_pyl = el.Py();
       genie_pzl = el.Pz();
       genie_pl = el.Rho();
       genie_cthl = cos(el.Theta());
       genie_vtxx = 0;
       genie_vtxy = 0;
       genie_vtxz = 0;


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

 	     genie_pdgf[IndexInFinalStateParticleArray] = 2212;
 	     genie_Ef[IndexInFinalStateParticleArray] = pr.E();
 	     genie_pxf[IndexInFinalStateParticleArray] = pr.Px();
 	     genie_pyf[IndexInFinalStateParticleArray] = pr.Py();
 	     genie_pzf[IndexInFinalStateParticleArray] = pr.Pz();
 	     genie_pf[IndexInFinalStateParticleArray] = pr.Rho();
 	     genie_cthf[IndexInFinalStateParticleArray] = pr.CosTheta();
	     
 	     //double CorrectedProtonMomentum = ProtonMomCorrection_He3_4Cell(ftarget,V4_uncorrprot,p_vert_corr);
 	     //V4_uncorrprot.SetRho(CorrectedProtonMomentum);
 	     //double CorrectedProtonE = sqrt(CorrectedProtonMomentum*CorrectedProtonMomentum + m_prot*m_prot);
 	     //V4_uncorrprot.SetE(CorrectedProtonE);
	     
 	     genie_pdgf[IndexInFinalStateParticleArray+shift] = 2212;
 	     genie_Ef[IndexInFinalStateParticleArray+shift] = pr.E();
 	     genie_pxf[IndexInFinalStateParticleArray+shift] = pr.Px();
 	     genie_pyf[IndexInFinalStateParticleArray+shift] = pr.Py();
 	     genie_pzf[IndexInFinalStateParticleArray+shift] = pr.Pz();
 	     genie_pf[IndexInFinalStateParticleArray+shift] = pr.Rho();
 	     genie_cthf[IndexInFinalStateParticleArray+shift] = pr.CosTheta();
 	     IndexInFinalStateParticleArray++;
	     

	     nprot++;
	     loc_proton_v4.push_back(pr);
	   }
	   
 	   if(p->par()->getPid()==211){
 	     SetLorentzVector(pip,pips[npip]);

 	     double PiPlusP = pip.Vect().Mag();
 	     double PiPlusE = TMath::Sqrt(PiPlusP*PiPlusP + m_pipl * m_pipl);
	     
 	     genie_pdgf[IndexInFinalStateParticleArray] = 211;
 	     genie_Ef[IndexInFinalStateParticleArray] = PiPlusE;
 	     genie_pxf[IndexInFinalStateParticleArray] = pip.Vect().X();
 	     genie_pyf[IndexInFinalStateParticleArray] = pip.Vect().Y();
 	     genie_pzf[IndexInFinalStateParticleArray] = pip.Vect().Z();
 	     genie_pf[IndexInFinalStateParticleArray] = PiPlusP;
 	     genie_cthf[IndexInFinalStateParticleArray] = pip.Vect().CosTheta();
 	     IndexInFinalStateParticleArray++;
	     
 	     npip++;
 	     loc_pipl_v4.push_back(pip);
 	   }

	   if(p->par()->getPid()==321){
	     SetLorentzVector(kpl,kpls[nkpl]);

 	     genie_pdgf[IndexInFinalStateParticleArray] = 321;
 	     genie_Ef[IndexInFinalStateParticleArray] = kpl.E();
 	     genie_pxf[IndexInFinalStateParticleArray] = kpl.Px();
 	     genie_pyf[IndexInFinalStateParticleArray] = kpl.Py();
 	     genie_pzf[IndexInFinalStateParticleArray] = kpl.Pz();
 	     genie_pf[IndexInFinalStateParticleArray] = kpl.Rho();
 	     genie_cthf[IndexInFinalStateParticleArray] = kpl.CosTheta();
 	     IndexInFinalStateParticleArray++;
	     

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
	     
 	     genie_pdgf[IndexInFinalStateParticleArray] = -211;
 	     genie_Ef[IndexInFinalStateParticleArray] = PiMinusE;
 	     genie_pxf[IndexInFinalStateParticleArray] = pim.Vect().X();
 	     genie_pyf[IndexInFinalStateParticleArray] = pim.Vect().Y();
 	     genie_pzf[IndexInFinalStateParticleArray] = pim.Vect().Z();
 	     genie_pf[IndexInFinalStateParticleArray] = PiMinusP;
 	     genie_cthf[IndexInFinalStateParticleArray] = pim.Vect().CosTheta();
 	     IndexInFinalStateParticleArray++;


 	     npim++;
 	     loc_pimi_v4.push_back(pim);
 	   }

	   if(p->par()->getPid()==-321){
	     SetLorentzVector(kmi,kmis[nkmi]);

 	     genie_pdgf[IndexInFinalStateParticleArray] = -321;
 	     genie_Ef[IndexInFinalStateParticleArray] = kmi.E();
 	     genie_pxf[IndexInFinalStateParticleArray] = kmi.Px();
 	     genie_pyf[IndexInFinalStateParticleArray] = kmi.Py();
 	     genie_pzf[IndexInFinalStateParticleArray] = kmi.Pz();
 	     genie_pf[IndexInFinalStateParticleArray] = kmi.Rho();
 	     genie_cthf[IndexInFinalStateParticleArray] = kmi.CosTheta();
 	     IndexInFinalStateParticleArray++;
	     

	     nkmi++;
	     loc_kmi_v4.push_back(kpl);
	   }
	 }
	 
	 
	 //what about neutrals???
	 
       }


    //---------------------------------------------
    //----------END OF HADRON SELECTION------------
    //---------------------------------------------

       genie_nfp = nprot;
       genie_nfn = 0;
       genie_nfpip = npip;
       genie_nfpim = npim;

       genie_nfkp = nkpl;
       genie_nfkm = nkmi;

       genie_nfpi0 = 0;
       //genie_nfpi0 = ec_num_n;

       //genie_nf = genie_nfp + genie_nfn + genie_nfpip + genie_nfpim + genie_nfpi0;
       genie_nf = genie_nfp + genie_nfn + genie_nfpip + genie_nfpim + genie_nfkp + genie_nfkm + genie_nfpi0;
       
       genie_iev = NEventsTotal;
       NEventsTotal++;


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
void RGMtoGenie(){

  // Data files to process
  //TString inFile("/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v0/dst/train/skim4/skim4_005*.hipo");
  //TString inFile("/volatile/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v0/dst/train/dst/train/skim4/*.hipo");
  TString inFile("/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v0/dst/train/skim4/skim4_005051.hipo");
  // TString inputFile2("/volatile/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train/inc/*.hipo");

  TString outFile("/home/sfegan/CTOF_rcdbTest5051.root");
  
  const std::string databaseF("/home/sfegan/rcdb.root");
  
  cout << "calling function" <<endl;
  
  RGMtoGenie(inFile, outFile, databaseF); //call the analysis function with this filename 

}

