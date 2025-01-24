#include "CaloCalibration/SourceCalib/inc/SourceFitter.hh"

using namespace TMath;
using namespace RooFit;
using namespace CaloSourceCalib;

void SourceFitter::FitCrystal(TH1F* h_spec, TString opt, int crystalNo,  TTree *covar, Float_t &fpeak, Float_t &fsigma, Float_t &chiSq, Float_t &fstpeak,Float_t &fstsigma, Float_t &scdpeak,Float_t &scdsigma,Float_t &fcbalphaparam,Float_t &fcbndegparam,Float_t &Aparam,Float_t &Bparam, Float_t &Cparam, Float_t &fullResparam, Float_t &fstResparam,Float_t &scdResparam,Float_t &comCnstparam, Float_t &combetaparam, Float_t &frFullparam, Float_t &frFrstparam,Float_t &frScndparam,Float_t &crystalNoparam,Float_t &frBKGparam, Float_t &calibconst){
    
  // set stlye optionsr
  gStyle -> SetOptFit(1111);
  gStyle -> SetOptStat(0);
  gStyle -> SetPadBottomMargin(0.125);
  gStyle -> SetPadTopMargin(0.075);
  gStyle -> SetPadRightMargin(0.05);
  gStyle -> SetPadLeftMargin(0.15);
  gStyle -> SetTitleOffset(1.0, "x");
  gStyle -> SetTitleOffset(1.75, "y");

  TCanvas *can = new TCanvas("can", "", 100, 100, 600, 600);
  can -> Draw();
  TString cryNum = to_string(crystalNo);
  TString oName = "mu2e_simu_fitSpec_"+ opt+"_" + cryNum + ".root";
  TString title = "Energy Spectrum of Crystal " + cryNum;

  TFile *ouptFile = new TFile("mu2e_caloSourceFit_" + cryNum + ".root", "RECREATE");

  // obtain the two initial values for the parameters
  double par1 = 2500;
  double par2 = 50;
  //double ADC_conv = 0.0625;

  //parameters
  RooRealVar crysADC("crysADC", "ADC[counts]", 48, 120);
  RooRealVar ergElec("ergElec", "electron energy in ADC", 8.176);
	RooRealVar fcbalpha("fcbalpha", "fcbalpha", 0.5, 0.1, 5);
  RooRealVar fcbndeg("fcbndeg", "fcbndeg", 10, 1, 40.);
  RooRealVar A("A", "coeff of E", 0.02,0.001,3);//0.16, 0.08, 16);0.006
  RooRealVar B("B", "const B",0.00405,0.0001,2);//8, 0.08, 16);0.00405
  RooRealVar C("C", "const C", 0.0027,0.0001,3);//Electronic noise in MeV 2.4  0.0027
  
  //Full peak:
  RooRealVar fullPeak("fullPeak", "full peak", 96.8, 80, 108);
  RooFormulaVar fullRes("fullRes","full peak resolution","0.98650796*sqrt(pow(A/pow(fullPeak/1000,0.25),2)+pow(B/0.0625,2)+pow(C/(fullPeak/1000),2))", RooArgSet(A,B,C,fullPeak));    
  RooFormulaVar fullWidth("fullWidth", "width of the full peak", "fullPeak*fullRes",RooArgSet(fullPeak, fullRes));
  
  // 1st escape peak:
  RooFormulaVar fstEsPeak("fstEsPeak", "first escape peak", "fullPeak - ergElec", RooArgSet(fullPeak,ergElec));
  RooFormulaVar fstRes("fstRes","first peak resolution","0.98650796*sqrt(pow(A/pow(fstEsPeak /1000,0.25),2)+pow(B/0.0625,2)+pow(C/(fstEsPeak/1000),2))", RooArgSet(A,B,C,fstEsPeak));
  RooFormulaVar fstWidth("fstWidth", "width of first escape peak","fstEsPeak*fstRes",RooArgSet(fstEsPeak,fstRes));

  // 2nd escape peak:
  RooFormulaVar scdEsPeak("scdEsPeak", "second escape peak", "fullPeak - 2*ergElec", RooArgSet(fullPeak,ergElec));
  RooFormulaVar scdRes("scdRes","second peak resolution","0.98650796*sqrt(pow(A/pow(scdEsPeak/1000,0.25),2)+pow(B/0.0625,2)+pow(C/(scdEsPeak/1000),2))", RooArgSet(A,B,C,scdEsPeak));
  RooFormulaVar scdWidth("scdWidth", "width of second escape peak","scdEsPeak*scdRes",RooArgSet(scdEsPeak,scdRes));

  // construct CBs:
  RooCBShape fullErg("fullErg", "full energy peak", crysADC, fullPeak, fullWidth, fcbalpha, fcbndeg);
  RooCBShape firsErg("firsErg", "first escape peak", crysADC, fstEsPeak, fstWidth, fcbalpha, fcbndeg);
  RooCBShape secdErg("secdErg", "second escape peak", crysADC, scdEsPeak, scdWidth, fcbalpha, fcbndeg);
  
  // logistic background:
  RooRealVar comCnst("comCnst", "comCnst", par1, 160, 5120);
  RooRealVar combeta("combeta", "combeta", par2, 4, 480);
  RooGenericPdf comPdf("comPdf", "logistic", "1.0/(1.0+exp((crysADC-comCnst)/combeta))",
             RooArgSet(crysADC, comCnst, combeta));

  // Fraction of events in each peak
  RooRealVar frFull("frFull", "Fraction of full peak",0.3,0.25, 1);
  RooRealVar frFrst("frFrst", "Fraction of first escape peak", 0.5,0.1, 1);
  RooRealVar frScnd("frScnd", "Fraction of second escape peak", 0.1, 0.1, 1);

  Float_t dpeak;
  //preparing RooPlot
  RooPlot *chFrame = crysADC.frame(Title(title));
  RooDataHist chSpec("crysADC","crysADC", crysADC, h_spec);
  // combined fit function
  RooAddPdf fitFun("fitFun", "firsErg + (secdErg + (fullErg + comPdf))", RooArgList(firsErg, secdErg, fullErg, comPdf), RooArgList(frFrst, frScnd, frFull) );

  if(opt == "chi2"){ //binned chi2 fit
    RooFitResult *fitRes = fitFun.chi2FitTo(chSpec, Range(48,115.2), Strategy(1),Hesse(kTRUE), PrintLevel(1),Save());
    fitRes->Print("v");
   }  
  if(opt == "nll"){ //binned nll fit

    //RooAddPdf fitFun("fitFun", "firsErg + (secdErg + (fullErg + comPdf))", RooArgList(firsErg, secdErg, fullErg, comPdf), RooArgList(frFrst, frScnd, frFull), kTRUE);
    RooAbsReal* nll = fitFun.createNLL(chSpec, Range(48,115.2)); 
    RooMinimizer m(*nll);
    m.migrad();
    m.hesse();
    RooFitResult *fitRes = m.save();
    fitRes->Print("v");
  }
  // plot components
  chSpec.plotOn(chFrame, MarkerColor(kBlack), LineColor(kBlack), MarkerSize(0.5), Name("chSpec"));
  fitFun.plotOn(chFrame, LineColor(kRed), LineStyle(1), Name("fullFit"));
  chiSq = chFrame->chiSquare(11);
  fitFun.plotOn(chFrame, Components(fullErg), LineColor(kOrange), LineStyle(5));
  fitFun.plotOn(chFrame, Components(firsErg), LineColor(kViolet), LineStyle(5));
  fitFun.plotOn(chFrame, Components(secdErg), LineColor(kCyan), LineStyle(5));
  fitFun.plotOn(chFrame, Components(comPdf), LineColor(kBlue), LineStyle(5));
	

  fpeak = fullPeak.getVal();
  dpeak = fullPeak.getError();
  fsigma = fullWidth.getVal();
  fstpeak = fstEsPeak.getVal();
  fstsigma = fstWidth.getVal();
  scdpeak = scdEsPeak.getVal();
  scdsigma = scdWidth.getVal();
  fcbalphaparam = fcbalpha.getVal();
  fcbndegparam = fcbndeg.getVal();
  Aparam = A.getVal();
  Bparam = B.getVal();
  Cparam = C.getVal();
  fullResparam = fullRes.getVal();
  fstResparam = fstRes.getVal();
  scdResparam = scdRes.getVal();
  comCnstparam = comCnst.getVal();
  combetaparam = combeta.getVal();
  frFullparam = frFull.getVal();
  frFrstparam = frFrst.getVal();
  frScndparam = frScnd.getVal();
  frBKGparam= 1-(frFullparam+frFrstparam+frScndparam);
  std::cout<< "fraction of background" <<frBKGparam <<std::endl;
  crystalNoparam = crystalNo;
  calibconst = 6.13/fpeak;
  std::cout<< "calib const" <<calibconst <<std::endl;

  //make pretty plots
  TPaveLabel *pchi2 = new TPaveLabel(0.20, 0.70, 0.35, 0.80, Form("#chi^{2}/ndf = %4.2f", chiSq), "brNDC");
  pchi2 -> SetFillStyle(0);
  pchi2 -> SetBorderSize(0);
  pchi2 -> SetTextSize(0.25);
  pchi2 -> SetTextColor(kBlack);
  pchi2 -> SetFillColor(kWhite);
  chFrame -> addObject(pchi2);
  TPaveLabel *fpk = new TPaveLabel(0.20, 0.85, 0.35, 0.95, Form("peak = %.2f#pm%.2f", fpeak, dpeak), "brNDC");
  fpk -> SetFillStyle(0);
  fpk -> SetBorderSize(0);
  fpk -> SetTextSize(0.25);
  fpk -> SetTextColor(kBlack);
  fpk -> SetFillColor(kWhite);
  chFrame -> addObject(fpk);
  TPaveLabel *fsg = new TPaveLabel(0.20, 0.775, 0.35, 0.875, Form("sigma = %4.2f", fsigma), "brNDC");//dsigma
  fsg -> SetFillStyle(0);
  fsg -> SetBorderSize(0);
  fsg -> SetTextSize(0.25);
  fsg -> SetTextColor(kBlack);
  fsg -> SetFillColor(kWhite);
  chFrame -> addObject(fsg);
  std::cout << "chi2: " << chiSq << "; Probability: " << Prob(chiSq, 151) << std::endl;

  chFrame -> SetYTitle("Events per 25 keV");
  chFrame -> Draw();
  can -> SaveAs(oName);

  covar->Fill();
  ouptFile -> Write();
  ouptFile -> Close();
}
//MC source fitter-- will be deleted before push
void SourceFitter::MCFitCrystal(int crystalNo, TString opt){//Float_t &fpeak_MC,Float_t &fsigma_MC,Float_t &chiSq_MC

    // set stlye options
    gStyle -> SetOptFit(1111);
    gStyle -> SetOptStat(0);
    gStyle -> SetPadBottomMargin(0.125);
    gStyle -> SetPadTopMargin(0.075);
    gStyle -> SetPadRightMargin(0.05);
    gStyle -> SetPadLeftMargin(0.15);
    gStyle -> SetTitleOffset(1.0, "x");
    gStyle -> SetTitleOffset(1.75, "y");


    //for(unsigned int crystalNo = startCry; crystalNo < endCry; crystalNo++){

    TCanvas *can_MC = new TCanvas("can_MC", "", 100, 100, 600, 600);
    can_MC -> Draw();
    TString cryNum_MC = to_string(crystalNo);
    TString fName_MC = "mu2e_caloSimu_crysEdep_"+cryNum_MC+".root" ;
    TString oName_MC = "mu2e_simu_MCfitSpec_"+ opt+"_" + cryNum_MC + ".root";
    TString title_MC = "Energy Spectrum of Crystal " + cryNum_MC;
    TFile *inFile_MC = TFile::Open(fName_MC);

    TFile *ouptFile_MC = new TFile("mu2e_MCcaloSourceFit_" + cryNum_MC + ".root", "RECREATE");

    
    // obtain the two initial values for the parameters
    double par1_MC = 100.;
    double par2_MC = 10.;

    //parameters
    RooRealVar crysEdep_MC("crysEdep_MC", "E_{reco}_MC [MeV]", 3.0, 7.0);
    //RooRealVar ergElec_MC("ergElec_MC", "electron energy_MC", 0.511);
    RooRealVar fcbalpha_MC("fcbalpha_MC", "fcbalpha_MC", 2.5, 0.05, 20.0);
    RooRealVar fcbndeg_MC("fcbndeg_MC", "fcbndeg_MC", 10., 0.25, 80.);
    RooRealVar A_MC("A_MC", "coeff of E_MC", 0.0003,0.00001,10);
    RooRealVar B_MC("B_MC", "const B_MC",0.00253125,0.00001,10);
    RooRealVar C_MC("C_MC", "const C_MC", 0.00027,0.00001,10);
    

    //Full peak:
    RooRealVar fullPeak_MC("fullPeak_MC", "full peak_MC", 6.05, 5.0, 6.75);
    RooFormulaVar fullRes_MC("fullRes_MC","full peak resolution_MC","sqrt(pow(A_MC/pow(6.13/1000,0.25),2)+pow(B_MC/0.0625,2)+pow(C_MC/(6.13/1000),2))", RooArgSet(A_MC,B_MC,C_MC));    
    RooFormulaVar fullWidth_MC("fullWidth_MC", "width of the full peak_MC", "fullPeak_MC*fullRes_MC",RooArgSet(fullPeak_MC, fullRes_MC));

    // 1st escape peak:
    RooFormulaVar fstEsPeak_MC("fstEsPeak_MC", "first escape peak_MC", "fullPeak_MC - 0.511", RooArgSet(fullPeak_MC));
    RooFormulaVar fstRes_MC("fstRes_MC","first peak resolution_MC","sqrt(pow(A_MC/pow(5.619/1000,0.25),2)+pow(B_MC/0.0625,2)+pow(C_MC/(5.619/1000),2))", RooArgSet(A_MC,B_MC,C_MC));
    RooFormulaVar fstWidth_MC("fstWidth_MC", "width of first escape peak_MC","fstEsPeak_MC*fstRes_MC",RooArgSet(fstEsPeak_MC,fstRes_MC));

    // 2nd escape peak:
    RooFormulaVar scdEsPeak_MC("scdEsPeak_MC", "second escape peak_MC", "fullPeak_MC - 2*0.511", RooArgSet(fullPeak_MC));
    RooFormulaVar scdRes_MC("scdRes_MC","second peak resolution_MC","sqrt(pow(A_MC/pow(5.1075/1000,0.25),2)+pow(B_MC/0.0625,2)+pow(C_MC/(5.1075/1000),2))", RooArgSet(A_MC,B_MC,C_MC));
    RooFormulaVar scdWidth_MC("scdWidth_MC", "width of second escape peak_MC","scdEsPeak_MC*scdRes_MC",RooArgSet(scdEsPeak_MC,scdRes_MC));

    // construct CBs:
    RooCBShape fullErg_MC("fullErg_MC", "full energy peak_MC", crysEdep_MC, fullPeak_MC, fullWidth_MC, fcbalpha_MC, fcbndeg_MC);
    RooCBShape firsErg_MC("firsErg_MC", "first escape peak_MC", crysEdep_MC, fstEsPeak_MC, fstWidth_MC, fcbalpha_MC, fcbndeg_MC);
    RooCBShape secdErg_MC("secdErg_MC", "second escape peak_MC", crysEdep_MC, scdEsPeak_MC, scdWidth_MC, fcbalpha_MC, fcbndeg_MC);

    // logistic background:
    RooRealVar comCnst_MC("comCnst_MC", "comCnst_MC", par1_MC, 10, 320.);
    RooRealVar combeta_MC("combeta_MC", "combeta_MC", par2_MC, 0.25, 30);
    RooGenericPdf comPdf_MC("comPdf_MC", "logistic_MC", "1.0/(1.0+exp((crysEdep_MC-comCnst_MC)/combeta_MC))",
               RooArgSet(crysEdep_MC, comCnst_MC, combeta_MC));

    // Fraction of events in each peak
    RooRealVar frFull_MC("frFull_MC", "Fraction of full peak_MC", 0.5, 0.25, 0.7);
    RooRealVar frFrst_MC("frFrst_MC", "Fraction of first escape peak_MC", 0.3, 0.1, 0.5);
    RooRealVar frScnd_MC("frScnd_MC", "Fraction of second escape peak_MC", 0.2, 0.1, 0.5);
    RooAddition constrSum_MC("add","add", RooArgList(frFull_MC,frFrst_MC,frScnd_MC));
    
    //RooRealVar frBKG("frBKG", "Fraction of BKG peak", 0.05, 0.1, 0.0);

    // the parameter is 11 after setting width of escape peaks opened
    Float_t fsigma_MC, fpeak_MC,dpeak_MC;//fsigma, dsigma, fpeak,
    Float_t chiSq_MC = 0;
    //preparing RooPlot
    RooPlot *chFrame_MC = crysEdep_MC.frame(Title(title_MC));

   if(opt == "chi2_MC"){ //binned chi2 fit
      TH1F *h_spec = (TH1F*)inFile_MC->Get("h_spec");//simEdep
      RooDataHist chSpec_MC("chSpec_MC", "chSpec", crysEdep_MC, h_spec);
      // combined fit function
      RooAddPdf fitFun_MC("fitFun", "firsErg_MC + (secdErg_MC + (fullErg_MC + comPdf_MC))", RooArgList(firsErg_MC, secdErg_MC, fullErg_MC, comPdf_MC), RooArgList(frFrst_MC, frScnd_MC, frFull_MC) );

      //chi2 fit
      //RooFitResult *fitRes = fitFun.chi2FitTo(chSpec, Range(3.0, 7.2), Strategy(3), PrintLevel(1), Hesse(kTRUE), Extended(),Save());
       RooChi2Var chi2_MC("chi2_MC", "chi2_MC", fitFun_MC, chSpec_MC, Range(48, 115.2), Extended());
       RooMinimizer minimizer_MC(chi2_MC);
			 minimizer_MC.minimize("Minuit2", "Migrad");
			 minimizer_MC.minos();
			 //minimizer.hesse();
			 RooFitResult* fitRes_MC = minimizer_MC.save();
			 fitRes_MC->Print("v");
      // plot components
      chSpec_MC.plotOn(chFrame_MC, MarkerColor(kBlack), LineColor(kBlack), MarkerSize(0.5), Name("chSpec"));
      fitFun_MC.plotOn(chFrame_MC, LineColor(kRed), LineStyle(1), Name("fullFit"));
      chiSq_MC = chFrame_MC->chiSquare(11);
      fitFun_MC.plotOn(chFrame_MC, Components(fullErg_MC), LineColor(kOrange), LineStyle(5));
      fitFun_MC.plotOn(chFrame_MC, Components(firsErg_MC), LineColor(kViolet), LineStyle(5));
      fitFun_MC.plotOn(chFrame_MC, Components(secdErg_MC), LineColor(kCyan), LineStyle(5));
      fitFun_MC.plotOn(chFrame_MC, Components(comPdf_MC), LineColor(kBlue), LineStyle(5));
    } 
    
    //get fit parameters
    fpeak_MC = fullPeak_MC.getVal();
    dpeak_MC = fullPeak_MC.getError();
    fsigma_MC = fullWidth_MC.getVal();

    //make pretty plots
    TPaveLabel *pchi2_MC = new TPaveLabel(0.20, 0.70, 0.35, 0.80, Form("#chi^{2}/ndf = %4.2f", chiSq_MC), "brNDC");
    pchi2_MC -> SetFillStyle(0);
    pchi2_MC -> SetBorderSize(0);
    pchi2_MC -> SetTextSize(0.25);
    pchi2_MC -> SetTextColor(kBlack);
    pchi2_MC -> SetFillColor(kWhite);
    chFrame_MC -> addObject(pchi2_MC);
    TPaveLabel *fpk_MC = new TPaveLabel(0.20, 0.85, 0.35, 0.95, Form("peak = %.2f#pm%.2f", fpeak_MC, dpeak_MC), "brNDC");
    fpk_MC -> SetFillStyle(0);
    fpk_MC -> SetBorderSize(0);
    fpk_MC -> SetTextSize(0.25);
    fpk_MC -> SetTextColor(kBlack);
    fpk_MC -> SetFillColor(kWhite);
    chFrame_MC -> addObject(fpk_MC);
    TPaveLabel *fsg_MC = new TPaveLabel(0.20, 0.775, 0.35, 0.875, Form("sigma = %4.2f", fsigma_MC), "brNDC");
    fsg_MC -> SetFillStyle(0);
    fsg_MC -> SetBorderSize(0);
    fsg_MC -> SetTextSize(0.25);
    fsg_MC -> SetTextColor(kBlack);
    fsg_MC -> SetFillColor(kWhite);
    chFrame_MC -> addObject(fsg_MC);
    std::cout << "chi2_MC: " << chiSq_MC << "; Probability: " << Prob(chiSq_MC, 151) << std::endl;

    chFrame_MC -> SetYTitle("Events per 25 keV");
    chFrame_MC -> Draw();
    can_MC -> SaveAs(oName_MC);

    // the file must be closed
    inFile_MC -> Close();

    ouptFile_MC -> Write();
    ouptFile_MC -> Close();
}
