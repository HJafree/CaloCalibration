#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <iostream>
#include <vector>

TString filepath = "/exp/mu2e/app/users/hjafree/SourceFitDir/paraFile.root";

/*TTree* get_data(){
    TFile *f =  new TFile(filepath);
    TTree *t = (TTree*)f->Get("covar");
    return t;
}*/
TTree* get_data() {
    TFile* f = new TFile(filepath, "READ");
    if (f->IsZombie()) {
        std::cerr << "Error opening file: " << filepath << std::endl;
        delete f;
        return nullptr;
    }
    TTree* t = (TTree*)f->Get("covar");
    if (!t) {
        std::cerr << "Error: Tree 'covar' not found in file." << std::endl;
        delete f;
        return nullptr;
    }
    return t;
  std::cout<< "line30" <<std::endl;
}
void ScatterPlot(const char* filepath) {
  std::cout<< "line33" <<std::endl;
    // Open the ROOT file
    /*TFile file(filepath, "READ");
    TTree* tree = (TTree*)file.Get("covar"); 
    double crystalNo, Peak;
    tree->SetBranchAddress("crystalNo", &crystalNo);
    tree->SetBranchAddress("Peak", &Peak);*/

    /*TString ouptName = "param_distribution_plots.root";
    TFile *ouptFile = new TFile(ouptName, "RECREATE");
    Float_t paramdist;
    TTree *outTree = new TTree("paramTree","paramTree");
    outTree->Branch("paramdist", &paramdist, "paramdist/F"); 

    // input tree
    TTree *inTree = get_data_tree();
    int  crystalNo;
    float Peak;
    inTree -> SetBranchAddress("cryId", &crystalNo);
    inTree -> SetBranchAddress("peak", &Peak);
    unsigned int nEvt = (int)inTree -> GetEntries();
    unsigned int nEvtCrys = 0; //events in this crystal
    for(unsigned int iEvt=0; iEvt<nEvt; iEvt++)
    {
      // extract single entry in the Tree
      inTree -> GetEntry(iEvt);
      std::vector<double> crystalNos, Peaks;
      
      for(int icry=0; icry<crystalNo; icry++)
      {
        	crystalNo.push_back(crystalNo);
        	Peak.push_back(Peak);
        
      }

      outTree->Fill();
      
    }
    std::cout<<" Events analyzed in this crystal : "<<nEvtCrys<<std::endl;
  ouptFile-> Write();
  ouptFile -> Close();
}*/
    TString ouptName = "param_distribution_plots.root";
    TFile *ouptFile = new TFile(ouptName, "RECREATE");
    
    Float_t paramdist;
    TTree *outTree = new TTree("paramTree","paramTree");
    outTree->Branch("paramdist", &paramdist, "paramdist/F");
     
    TTree *inTree = get_data();
    //if (!inTree) return;

    int  crystalNo;
    float Peak;
    inTree -> SetBranchAddress("crystalNo", &crystalNo);
    inTree -> SetBranchAddress("Peak", &Peak);
    // Vectors to hold data 
    std::vector<double> crystalNos, Peaks;   
    // Fill the vectors with data from the TTree
    Long64_t nEntries = inTree->GetEntries();
    for (Long64_t entry = 0; entry < nEntries; ++entry) {

        inTree->GetEntry(entry);
        crystalNos.push_back(crystalNo);
        Peaks.push_back(Peak);
    
     paramdist = Peak;
     outTree->Fill();
}
    // Create a TGraph for the scatter plot
    TGraph graph(crystalNos.size(), crystalNos.data(), Peaks.data());

    // Customize the plot
    graph.SetTitle("Cry Number vs Full Peak;Crystal Number;Full Peak");
    graph.SetMarkerStyle(20);
    graph.SetMarkerSize(0.8);
    graph.SetMarkerColor(kBlue);

    // Create a canvas to draw the plot
    TCanvas canvas("canvas", "Scatter Plot", 800, 600);
    graph.Draw("AP"); // A for axis, P for points
  	//canvas -> Draw();
    // Save the canvas
    canvas.SaveAs("scatter_plot.eps");
 		ouptFile -> Write();
  	ouptFile -> Close();
  	//delete ouptFile;
  	//delete inTree->GetCurrentFile();

    // Close the file
   // filepath.Close();
}

