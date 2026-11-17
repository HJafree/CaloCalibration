// SourceCalibDigiAna: fills the source-calibration histograms from the crystals selected by
// SourceCalibSelector. All cuts live in the selector; this module only fills. One entry per
// selected crystal per time group (one CaloHit):
//   crystals_ADC/cry_i         crystal energy = average over its SiPMs that have a peak (ADC)
//   sipm_ADC/sipm_i            summed amplitude of that SiPM's peaks (ADC)
//   crystals_edep_truth/cry_i  energyDepG4 of every CaloShowerSim in that crystal in the event
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"

#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloRecoDigi.hh"
#include "Offline/MCDataProducts/inc/CaloShowerSim.hh"
#include "TH1F.h"
#include <iostream>
#include <map>
#include <string>

constexpr int ncrystals = 1348;
constexpr int nsipms = 2696;

//std::vector<int> badcrys = {774, 775,776,792, 793,794,795,796,798,814,815,816,818,819,820,822,851,865,881,941, 965,969,1006,1029,1039,1071,1097,1144,1183,1209,1243,1275,1281,1316,1336};
namespace mu2e {

  class SourceCalibDigiAna : public art::EDAnalyzer
  {
     public:
        struct Config
        {
           using Name    = fhicl::Name;
           using Comment = fhicl::Comment;
           fhicl::Atom<art::InputTag>                          caloHitCollection        { Name("caloHitCollection"),        Comment("CaloHits of the selected crystals (SourceCalibSelector)") };
           fhicl::Atom<art::InputTag>                          caloShowerSimCollection  { Name("caloShowerSimCollection"),  Comment("Calo Sim module label") };
           fhicl::Atom<int>                                    diagLevel                { Name("diagLevel"),                Comment("Diagnosis level") };
           fhicl::Atom<double>                                 endBin                   { Name("endBin"),                   Comment("finalBin in ADC"),120 };
        };

        explicit SourceCalibDigiAna(const art::EDAnalyzer::Table<Config>& config) :
           EDAnalyzer{config},
           caloHitsToken_      {consumes<CaloHitCollection>(config().caloHitCollection())},
           caloShowerSimToken_ {consumes<CaloShowerSimCollection>(config().caloShowerSimCollection())},
           diagLevel_          (config().diagLevel()),
           endBin_             (config().endBin())
        {}
        virtual ~SourceCalibDigiAna() {}
        virtual void beginJob();
        virtual void analyze(const art::Event& e);

     private:
        const  art::ProductToken<CaloHitCollection>       caloHitsToken_;
        const  art::ProductToken<CaloShowerSimCollection> caloShowerSimToken_;
        int                                               diagLevel_;
        double                                            endBin_;
        TH1F* list_of_crys_hists[ncrystals];
        TH1F* list_of_sipm_hists[nsipms];
        TH1F* list_of_crys_hists_truth[ncrystals];
  };


  void SourceCalibDigiAna::analyze(const art::Event& event)
  {
      if (diagLevel_ > 0) std::cout<<"[SourceCalibDigiAna::analyze] begin"<<std::endl;
      const auto& caloHits = *event.getValidHandle(caloHitsToken_);
      const auto& caloSims = *event.getValidHandle(caloShowerSimToken_);

      for (unsigned int h = 0; h < caloHits.size(); ++h)
      {
          const auto& hit = caloHits[h];
          int id       = hit.crystalID();
          int id_sipm1 = id*2;
          int id_sipm2 = id*2+1;

          std::map<int, double> total_energy_in_sipm;             // summed amplitude per SiPM
          for (const auto& recoDigi : hit.recoCaloDigis()) total_energy_in_sipm[recoDigi->SiPMID()] += recoDigi->energyDep();

          std::cout << event.event() << " " << h << " " << id << " " << id_sipm1 << " " << id_sipm2 << "\n";
          list_of_crys_hists[id]->Fill(hit.energyDep());
          if(total_energy_in_sipm[id_sipm1]!=0) list_of_sipm_hists[id_sipm1]->Fill(total_energy_in_sipm[id_sipm1]);
          if(total_energy_in_sipm[id_sipm2]!=0) list_of_sipm_hists[id_sipm2]->Fill(total_energy_in_sipm[id_sipm2]);
          for (const auto& caloSim : caloSims)
          {
              if (caloSim.crystalID() == id) list_of_crys_hists_truth[id]->Fill(caloSim.energyDepG4());
          }
      }
      if (diagLevel_ > 0) std::cout<<"[SourceCalibDigiAna::analyze] end"<<std::endl;
  }

  void SourceCalibDigiAna::beginJob(){
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory crydir = tfs->mkdir( "crystals_ADC" );
      art::TFileDirectory crytruedir = tfs->mkdir( "crystals_edep_truth" );
      art::TFileDirectory sipmdir = tfs->mkdir( "sipm_ADC" );
      for(int i = 0; i < ncrystals ; i++){
        TString histname = "cry_"+std::to_string(i);
        list_of_crys_hists[i] = crydir.make<TH1F>( histname , histname, 120, 0.0, endBin_);
        list_of_crys_hists_truth[i] = crytruedir.make<TH1F>( histname , histname, 150, 0.0, 10.0); //assume MeV
      }
      for(int i = 0; i < nsipms ; i++){
        TString histname = "sipm_"+std::to_string(i);
        list_of_sipm_hists[i] = sipmdir.make<TH1F>( histname , histname, 120, 0.0, endBin_);
      }
      std::cout << "event hit crystal sipm1 sipm2\n";
  }
}
DEFINE_ART_MODULE(mu2e::SourceCalibDigiAna)
