#ifndef _SourceFitter_hh
#define _SourceFitter_hh

#include <fstream>
#include <iostream>
#include "TSystem.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"

#include "TF1.h"
#include "TMath.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMath.h"
#include <Riostream.h>

#include "TPaveStats.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPaveLabel.h"
#include "TAttFill.h"
// add roofit header files
#include "RooHist.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooKeysPdf.h"
#include "RooHistPdf.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooCBShape.h"
#include "RooGenericPdf.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
#include "RooUniform.h"
#include "RooMinimizer.h"
#include "RooAbsReal.h"
#include "RooFormulaVar.h"
#include "RooAddition.h"
#include "RooChi2Var.h" //included to add minimiser alternations
using namespace std;
using namespace TMath;
using namespace RooFit;

namespace CaloSourceCalib{
  class SourceFitter  {
      public:
        explicit SourceFitter(){};
        explicit SourceFitter(const SourceFitter &){};
        SourceFitter& operator = (const SourceFitter &);
        virtual ~SourceFitter() = default;
        #ifndef __CINT__
        void FitCrystal(TH1F* histogram, TString opt, int cryNum, TTree* covar, Float_t &fpeak, Float_t &fsigma, Float_t &chiSq, Float_t &fstpeak, Float_t &fstsigma, Float_t &scdpeak,Float_t &scdsigma,Float_t &fcbalphaparam,Float_t &fcbndegparam,Float_t &Aparam,Float_t &Bparam, Float_t &Cparam, Float_t &fullResparam, Float_t &fstResparam,Float_t &scdResparam,Float_t &comCnstparam, Float_t &combetaparam, Float_t &frFullparam, Float_t &frFrstparam,Float_t &frScndparam,Float_t &crystalNoparam,Float_t &frBKGparam,Float_t &calibconst);
        void MCFitCrystal(int crystalNo, TString opt);
        #endif
        ClassDef (SourceFitter,1);
    };
}
#endif
