#include "CaloCalibration/SourceCalib/inc/MakeAnalysisTree.hh"
#include "CaloCalibration/SourceCalib/inc/SourceFitter.hh"
#include <chrono>
using namespace std::chrono;

using namespace CaloSourceCalib;

TString filepath = "/pnfs/mu2e/tape/usr-nts/nts/sophie/SourceCalibSimAna/ADCCuts/root/56/e6/nts.sophie.SourceCalibSimAna.ADCCuts.0.root";

//TString filepath = "/exp/mu2e/app/users/hjafree/SourceFitDir/mu2e_caloSimu_crysEdep_770.root";
/*function to extract the TTree from the SourceCalibAna output*/
TH1F* get_data_histogram(int cryNum){
std::cout<< "line12" <<std::endl;
    TFile *f =  new TFile(filepath);
std::cout<< "line15" <<std::endl;
    TString crystalNumber = to_string(cryNum);
std::cout<< "line17" <<std::endl;
    TH1F* hist = (TH1F*)f->Get("CaloSourceCalibDigiAna/crystals_ADC/hspec_" + crystalNumber); 
std::cout<< "line19" <<std::endl;
    return hist;
std::cout<< "line21" <<std::endl; // NOT printing
}
//std::cout<< "line19" <<std::endl;
//TString truthfilepath = "/pnfs/mu2e/tape/usr-nts/nts/hjafree/SourceCalibAna/MDC2020ae/root/0a/1f/nts.hjafree.SourceCalibAna.MDC2020ae.0.root"; 

TString truthfilepath = "/exp/mu2e/app/users/hjafree/SourceFitDir/10M.root"; 
/*function to extract the TTree from the SourceCalibAna output*/
TTree* get_data_tree(){
    TFile *f =  new TFile(truthfilepath);
    TTree *t = (TTree*)f->Get("CaloExample/Calo");
    return t;
  std::cout<< "line28" <<std::endl;//NOT printing
}

//new function void-- loop over ttree (use above func to get ttree) use t2hf func to get histogram-- bin width is 674 bins 
/* bin for a given crystal */
/*void AnalyzeTruthCrystal(int crystalNo){

    TString cryNum;
    cryNum = to_string(crystalNo);

    // create the output file.
    TString trueouptName = "mu2e_caloSimu_crysEdep_" + cryNum + ".root";
    TFile *trueouptFile = new TFile(trueouptName, "RECREATE");
    //create output tree
    Float_t trueEdep;//crysEdep, ratio, ntrig, stim, time, tErg;
    TTree *outTree = new TTree("trueEdep","trueEdep");
    outTree->Branch("trueEdep", &trueEdep, "trueEdep/F"); //True Energy Deposited
//ToDO change things to make sure correct branch names
    // make histogram
    TH1F *true_spec = new TH1F("simEdep", "simEdep", 300, 0.0, 7.5);
    // input tree
    TTree *inTree = get_data_tree();
    // extract branches - these are arrays as there will be multiple crystal hits per event
    int  nCry;//number of crystal hits in event
    int cryId[10];
    float cryTime[10];
    float cryPosX[10], cryPosY[10], cryPosZ[10];
    float cryEdep[10];
    float simEdep[10];
    inTree -> SetBranchAddress("nCry", &nCry);//ncalhitHit
    inTree -> SetBranchAddress("cryId", &cryId);//same
    inTree -> SetBranchAddress("simEdep", &simEdep);//same
    inTree -> SetBranchAddress("cryEdep", &cryEdep);//calhitRecoEdep
    inTree -> SetBranchAddress("cryTime", &cryTime);//calhitRecoTime
    inTree -> SetBranchAddress("cryPosX", &cryPosX);//calhitRecoPosX
    inTree -> SetBranchAddress("cryPosY", &cryPosY);//calhitRecoPosY
    inTree -> SetBranchAddress("cryPosZ", &cryPosZ);//calhitRecoPosZ
   
    //unsigned int nEvt = (int)inTree -> GetEntries();
    unsigned int nEvtCrys = 0; //events in this crystal
    for(unsigned int iEvt=0; iEvt<1e7; iEvt++)//1e7
    {
      // extract single entry in the Tree
      inTree -> GetEntry(iEvt);
      
      int idExist = 0;
      std::vector<float> sameCryTime;
      std::vector<float> sameCryEdep;
      
      for(int icry=0; icry<nCry; icry++)
      {
        
        if(cryId[icry] == crystalNo)
        {
          idExist += 1;
          sameCryTime.push_back(cryTime[icry]);
          sameCryEdep.push_back(simEdep[icry]);
        }
      }

      if(idExist == 0) continue;
      sort(sameCryTime.begin(), sameCryTime.end());
      float edepTarget = 0.0;
      float edepOthers   = 0.0;

      std::vector<float> deltaTime;

      for(int icry=0; icry<nCry; icry++)
      {
        if(cryId[icry] == crystalNo)
        {
          edepTarget = simEdep[icry];
          deltaTime.push_back(cryTime[icry]);

        }
        else
        {
          edepOthers += simEdep[icry];
          deltaTime.push_back(cryTime[icry]);
        }
      }
      sort(deltaTime.begin(), deltaTime.end());
      float difTime = deltaTime.back() - deltaTime.front();
      
      if((edepTarget / (edepTarget + edepOthers)) >= 0.8 && difTime < 4)//0.8, 4
      {
        trueEdep =edepTarget;
        true_spec->Fill(edepTarget);
        nEvtCrys+=1;
      }
      outTree->Fill();
      
    }
    std::cout<<" Events analyzed in this crystal : "<<nEvtCrys<<std::endl;
  trueouptFile -> Write();
  trueouptFile -> Close();

}*/


void AnalyzeCrystal(int crystalNo){
    
    TString cryNum;
    cryNum = to_string(crystalNo);

    // create the output file.
    TString ouptName = "mu2e_caloSimu_crysEdep_" + cryNum + ".root";
    TFile *ouptFile_mev = new TFile(ouptName, "RECREATE");
    //create output tree
    Float_t crysEdep;//, ratio, ntrig, stim, time, tErg;
    TTree *outTree = new TTree("crysTree","crysTree");
    outTree->Branch("crysEdep", &crysEdep, "crysEdep/F"); //Energy Deposited
    
    // make histogram
    TH1F *h_spec = new TH1F("h_spec", "h_spec", 300, 0.0, 7.5);
    // input tree
    TTree *inTree = get_data_tree();

    // extract branches - these are arrays as there will be multiple crystal hits per event
    int  nCry;//number of crystal hits in event
    int cryId[10];
    float cryTime[10];
    float cryPosX[10], cryPosY[10], cryPosZ[10];
    float cryEdep[10];
    
    inTree -> SetBranchAddress("nCrystals", &nCry);//ncalhitHit
    inTree -> SetBranchAddress("cryId", &cryId);//same
    inTree -> SetBranchAddress("cryEdep", &cryEdep);//calhitRecoEdep
    inTree -> SetBranchAddress("cryTime", &cryTime);//calhitRecoTime
    inTree -> SetBranchAddress("cryPosX", &cryPosX);//calhitRecoPosX
    inTree -> SetBranchAddress("cryPosY", &cryPosY);//calhitRecoPosY
    inTree -> SetBranchAddress("cryPosZ", &cryPosZ);//calhitRecoPosZ

    unsigned int nEvt = (int)inTree -> GetEntries();
    unsigned int nEvtCrys = 0; //events in this crystal
    for(unsigned int iEvt=0; iEvt<nEvt; iEvt++)
    {
      // extract single entry in the Tree
      inTree -> GetEntry(iEvt);
      
      int idExist = 0;
      std::vector<float> sameCryTime;
      std::vector<float> sameCryEdep;
      
      for(int icry=0; icry<nCry; icry++)
      {
        
        if(cryId[icry] == crystalNo)
        {
          idExist += 1;
          sameCryTime.push_back(cryTime[icry]);
          sameCryEdep.push_back(cryEdep[icry]);
        }
      }

      if(idExist == 0) continue;
      sort(sameCryTime.begin(), sameCryTime.end());
      float edepTarget = 0.0;
      float edepOthers   = 0.0;

      std::vector<float> deltaTime;

      for(int icry=0; icry<nCry; icry++)
      {
        if(cryId[icry] == crystalNo)
        {
          edepTarget = cryEdep[icry];
          deltaTime.push_back(cryTime[icry]);

        }
        else
        {
          edepOthers += cryEdep[icry];
          deltaTime.push_back(cryTime[icry]);
        }
      }
      sort(deltaTime.begin(), deltaTime.end());
      float difTime = deltaTime.back() - deltaTime.front();
      
      if((edepTarget / (edepTarget + edepOthers)) >= 0.8 && difTime < 4)
      {
        crysEdep =edepTarget;
        h_spec->Fill(edepTarget);
        nEvtCrys+=1;
      }
      outTree->Fill();
      
    }
    std::cout<<" Events analyzed in this crystal : "<<nEvtCrys<<std::endl;
  ouptFile_mev -> Write();
  ouptFile_mev -> Close();
   std::cout<< "line219" <<std::endl;
}

  /*function calls RooFit and fits a given crystal file*/
void RunRooFit(int crystalNo) {
  SourceFitter *fit = new SourceFitter();
  fit->MCFitCrystal(crystalNo,"chi2_MC");
std::cout<< "line226" <<std::endl;//NOT printing
}
/*main function allows a loop over all crystals or a choice of a single crystal*/
int main(int argc, char* argv[]){
  std::cout<<"========== Welcome to the Mu2e Source Calibration Analysis =========="<<std::endl;
  int anacrys_start = 770; //starting crystal//675
  int anacrys_end = 800; //final crystal//680
  std::cout<<"choosing crystal "<<std::endl;
  if(strcmp( argv[1], "chooseCrystal") == 0 ){
    cout<<"crystal to be analyzed (int) : "<<endl;
    cin>>anacrys_start;
    anacrys_end = anacrys_start+1;
std::cout<< "line238" <<std::endl;//NOT printing
  }
  int debug = 0;
  if (debug == 0) {
 
		TFile *ouptFile = new TFile("paraFile.root", "RECREATE");
		Float_t fpeak, fsigma, chiSq, fstpeak, fstsigma, scdpeak,scdsigma,fcbalphaparam,fcbndegparam,Aparam,Bparam,Cparam,fullResparam,fstResparam,scdResparam,comCnstparam,
		combetaparam,frFullparam,frFrstparam,frScndparam,crystalNoparam,frBKGparam,
		calibconst;
		TTree *covar = new TTree("covar","Covariance Plot");
		covar->Branch("Peak", &fpeak,"fpeak/F");
		covar->Branch("Width", &fsigma,"fsigma/F");
		covar->Branch("ChiSq", &chiSq,"chiSq/F");
		covar->Branch("1stPeak", &fstpeak,"fstpeak/F");
		covar->Branch("1stWidth", &fstsigma,"fstsigma/F");
		covar->Branch("2ndPeak", &scdpeak,"scdpeak/F");
		covar->Branch("2ndWidth", &scdsigma,"scdsigma/F");
		covar->Branch("Alpha", &fcbalphaparam,"fcbalphaparam/F");
		covar->Branch("Ndeg", &fcbndegparam,"fcbndegparam/F");
		covar->Branch("A", &Aparam,"Aparam/F");
		covar->Branch("B", &Bparam,"Bparam/F");
		covar->Branch("C", &Cparam,"Cparam/F");
		covar->Branch("fullRes", &fullResparam,"fullResparam/F");
		covar->Branch("fstRes", &fstResparam,"fstResparam/F");
		covar->Branch("scdRes", &scdResparam,"scdResparam/F");
		covar->Branch("comCnst", &comCnstparam,"comCnstparam/F");
		covar->Branch("combeta", &combetaparam,"combetaparam/F");
		covar->Branch("frFull", &frFullparam,"frFullparam/F");
		covar->Branch("frFrst", &frFrstparam,"frFrstparam/F");
		covar->Branch("frScnd", &frScndparam,"frScndparam/F");
		covar->Branch("frBKG", &frBKGparam,"frBKGparam/F");
		covar->Branch("crystalNo", &crystalNoparam,"crystalNoparam/F");
		covar->Branch("1/calibconst", &calibconst,"calibconst/F");
		
	 for(int cryNum=770; cryNum<771; cryNum++){//674 --> 1348 
		TH1F* h = get_data_histogram(cryNum);
		std::cout<<"start value"<<anacrys_start<<std::endl;
		std::cout<<"end value"<<anacrys_end<<std::endl;
		std::cout<<"Running crystal binning ....."<<std::endl;
		//MakeCrystalListOutputs(anacrys_start,anacrys_end);
		//MakeCrystalBinsOutputs(anacrys_start,anacrys_end);
		SourceFitter *fit = new SourceFitter();
		fit->FitCrystal(h,"chi2", cryNum, covar, fpeak,fsigma,chiSq, fstpeak, fstsigma, scdpeak,scdsigma,fcbalphaparam,fcbndegparam,Aparam,Bparam,Cparam,
		fullResparam,fstResparam,scdResparam,comCnstparam,combetaparam,
		frFullparam,frFrstparam,frScndparam,crystalNoparam,frBKGparam,calibconst);
	 };
	 	ouptFile -> Write();
		ouptFile -> Close();
  std::cout<< "line286" <<std::endl;
		std::cout<<"Finished processing ..."<<std::endl;
		};
	if (debug == 1){		
	 for(int crystalNo = 770; crystalNo < 771; crystalNo++){
    //AnalyzeTruthCrystal(crystalNo);
    AnalyzeCrystal(crystalNo);    
    RunRooFit(crystalNo);  
std::cout<< "line294" <<std::endl;//NOT printing
};

			};
/*	if (debug == 2){
			
		TFile *ouptFile = new TFile("paraFile.root", "RECREATE");
		Float_t fpeak, fsigma, chiSq,fcbalphaparam,fcbndegparam,comCnstparam,
		combetaparam,crystalNoparam;
		TTree *covar = new TTree("covar","Covariance Plot");
		covar->Branch("Peak", &fpeak,"fpeak/F");
		covar->Branch("Width", &fsigma,"fsigma/F");
		
		
		
		covar->Branch("Alpha", &fcbalphaparam,"fcbalphaparam/F");
		covar->Branch("Ndeg", &fcbndegparam,"fcbndegparam/F");
		covar->Branch("comCnst", &comCnstparam,"comCnstparam/F");
		covar->Branch("combeta", &combetaparam,"combetaparam/F");
		covar->Branch("crystalNo", &crystalNoparam,"crystalNoparam/F");
		
		for(int cryNum=770; cryNum<771; cryNum++){//674 --> 1348 
		TH1F* h = get_data_histogram(cryNum);
		std::cout<<"start value"<<anacrys_start<<std::endl;
		std::cout<<"end value"<<anacrys_end<<std::endl;
		std::cout<<"Running crystal binning ....."<<std::endl;
		//MakeCrystalListOutputs(anacrys_start,anacrys_end);
		//MakeCrystalBinsOutputs(anacrys_start,anacrys_end);
		SourceFitter *fit = new SourceFitter();
		fit->CBFitCrystal(h,"chi2", cryNum, covar, fpeak,fsigma,chiSq,fcbalphaparam,fcbndegparam,
		comCnstparam,combetaparam,crystalNoparam);//add fractions
	 };
	 	ouptFile -> Write();
		ouptFile -> Close();

		std::cout<<"Finished processing ..."<<std::endl;
		
}*/
  std::cout<< "line332" <<std::endl;
  return 0;
}
