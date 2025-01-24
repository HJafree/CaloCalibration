#ifndef _SourcePlotter_hh
#define _SourcePlotter_hh

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
using namespace std;
using namespace TMath;
using namespace RooFit;

namespace CaloSourceCalib{
  class SourcePlotter  {
      public:
        explicit SourcePlotter(){};
        explicit SourcePlotter(const SourcePlotter &){};
        SourcePlotter& operator = (const SourcePlotter &);
        virtual ~SourcePlotter() = default;
        #ifndef __CINT__
        void ScatterPlot();
        #endif
        ClassDef (SourcePlotter,1);
    };
}
#endif
