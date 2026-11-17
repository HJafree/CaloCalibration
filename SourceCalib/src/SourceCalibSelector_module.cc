// SourceCalibSelector: source-calibration selection as a filter (waveform chi2, time grouping,
// isolation against the surrounding crystals). For every crystal that passes in a time group it
// writes one CaloHit (energy = average over the crystal's SiPMs that have a peak, in ADC) pointing to
// that crystal's CaloRecoDigis (one per waveform peak, amplitude in ADC). It also writes the selected
// CaloDigis for EventNtuple. The event passes if at least one crystal was selected.
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/RecoDataProducts/inc/CaloRecoDigi.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/CaloReco/inc/CaloWaveformProcessor.hh"
#include "Offline/CaloReco/inc/CaloTemplateWFProcessor.hh"
#include "Offline/CaloReco/inc/CaloRawWFProcessor.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace mu2e {

  class SourceCalibSelector : public art::EDFilter
  {
     public:
        enum processorStrategy {NoChoice, RawExtract, Template};

        struct Config
        {
           using Name    = fhicl::Name;
           using Comment = fhicl::Comment;
           fhicl::Table<mu2e::CaloRawWFProcessor::Config>      proc_raw_conf      { Name("RawProcessor"),       Comment("Raw processor config") };
           fhicl::Table<mu2e::CaloTemplateWFProcessor::Config> proc_templ_conf    { Name("TemplateProcessor"),  Comment("Log normal fit processor config") };
           fhicl::Atom<art::InputTag>                          caloDigiCollection { Name("caloDigiCollection"), Comment("Calo Digi module label") };
           fhicl::Atom<std::string>                            processorStrategy  { Name("processorStrategy"),  Comment("Digi reco processor name") };
           fhicl::Atom<double>                                 maxChi2Cut         { Name("maxChi2Cut"),         Comment("Max chi2/ndf of the waveform fit") };
           fhicl::Atom<double>                                 ratioCut           { Name("ratioCut"),           Comment("Min E_crystal/(E_crystal + E_surrounding crystals) in the time group") };
           fhicl::Atom<double>                                 timeCut            { Name("timeCut"),            Comment("Time-group window (ns), measured from the group's first digi") };
           fhicl::Atom<int>                                    diagLevel          { Name("diagLevel"),          Comment("Diagnosis level"), 0 };
        };

        explicit SourceCalibSelector(const art::EDFilter::Table<Config>& config);
        bool beginRun(art::Run& run) override;
        bool filter(art::Event& event) override;
        void endJob() override;

     private:
        art::ProductToken<CaloDigiCollection>   caloDigisToken_;
        std::string                             processorStrategy_;
        double                                  maxChi2Cut_;
        double                                  ratioCut_;
        double                                  timeCut_;
        int                                     diagLevel_;
        std::unique_ptr<CaloWaveformProcessor>  waveformProcessor_;
        unsigned long                           nEvents_       = 0;
        unsigned long                           nEventsPassed_ = 0;
        unsigned long                           nCrystalsIn_   = 0;
        unsigned long                           nCrystalsOut_  = 0;
        unsigned long                           nDigisIn_      = 0;
        unsigned long                           nDigisOut_     = 0;
        unsigned long                           nDigisNoPeak_  = 0;   // no fitted peak (too small to fit)
        unsigned long                           nDigisBadFit_  = 0;   // fitted, but chi2/ndf > maxChi2Cut
        unsigned long                           nDigisFitted_  = 0;   // fitted with a good chi2
  };

  SourceCalibSelector::SourceCalibSelector(const art::EDFilter::Table<Config>& config) :
     EDFilter{config},
     caloDigisToken_    {consumes<CaloDigiCollection>(config().caloDigiCollection())},
     processorStrategy_ (config().processorStrategy()),
     maxChi2Cut_        (config().maxChi2Cut()),
     ratioCut_          (config().ratioCut()),
     timeCut_           (config().timeCut()),
     diagLevel_         (config().diagLevel())
  {
      std::map<std::string, processorStrategy> spmap;
      spmap["RawExtract"]  = RawExtract;
      spmap["TemplateFit"] = Template;
      switch (spmap[processorStrategy_])
      {
          case RawExtract:
              waveformProcessor_ = std::make_unique<CaloRawWFProcessor>(config().proc_raw_conf());
              break;
          case Template:
              waveformProcessor_ = std::make_unique<CaloTemplateWFProcessor>(config().proc_templ_conf());
              break;
          default:
              throw cet::exception("CATEGORY") << "Unrecognized processor in SourceCalibSelector module";
      }
      produces<CaloDigiCollection>();
      produces<CaloRecoDigiCollection>();
      produces<CaloHitCollection>();
  }

  bool SourceCalibSelector::beginRun(art::Run&)
  {
      waveformProcessor_->initialize();
      return true;
  }

  bool SourceCalibSelector::filter(art::Event& event)
  {
      const auto caloDigisH = event.getValidHandle(caloDigisToken_);
      const auto& caloDigis = *caloDigisH;

      // 1. fit every digi, keep peaks with a good fit
      struct DigiHit { double time; int crystalID; int sipmID; double eDep; double eDepErr;
                       double peakTime; double peakTimeErr; double chi2; int ndf; bool pileUp; size_t digiIdx; };
      std::vector<DigiHit> hits;
      std::vector<double> x, y;

      for (size_t idx = 0; idx < caloDigis.size(); ++idx)
      {
          const auto& caloDigi = caloDigis[idx];
          double t0 = caloDigi.t0()*CaloConst::_digitizationPeriod;   // CaloDigi t0 is in digitizer clock ticks
          const std::vector<int>& waveform = caloDigi.waveform();

          x.clear(); y.clear();
          for (unsigned int i = 0; i < waveform.size(); ++i)
          {
              x.push_back(t0 + (i+0.5)*CaloConst::_digitizationPeriod);
              y.push_back(waveform.at(i));
          }

          waveformProcessor_->reset();
          waveformProcessor_->extract(x, y);

          if      (waveformProcessor_->nPeaks() == 0)                                          ++nDigisNoPeak_;
          else if (waveformProcessor_->chi2()/float(waveformProcessor_->ndf()) > maxChi2Cut_) ++nDigisBadFit_;
          else                                                                                ++nDigisFitted_;

          for (int i = 0; i < waveformProcessor_->nPeaks(); ++i)
          {
              double chi2 = waveformProcessor_->chi2();
              int    ndf  = waveformProcessor_->ndf();
              if (chi2/float(ndf) > maxChi2Cut_) continue;
              hits.push_back(DigiHit{t0, caloDigi.SiPMID()/2, caloDigi.SiPMID(),
                                     waveformProcessor_->amplitude(i), waveformProcessor_->amplitudeErr(i),
                                     waveformProcessor_->time(i), waveformProcessor_->timeErr(i),
                                     chi2, ndf, waveformProcessor_->isPileUp(i), idx});
          }
      }

      // 2. time groups: sort by digi t0, start a new group when a hit is more than timeCut_ after the group's first hit
      std::sort(hits.begin(), hits.end(), [](const DigiHit& a, const DigiHit& b){ return a.time < b.time; });

      std::vector<std::vector<DigiHit>> groups;
      for (const auto& hit : hits)
      {
          if (groups.empty() || hit.time - groups.back().front().time > timeCut_) groups.push_back({hit});
          else                                                                    groups.back().push_back(hit);
      }

      // 3. isolation: this crystal / (this crystal + its surrounding crystals) >= ratioCut_
      auto selectedDigis = std::make_unique<CaloDigiCollection>();
      auto recoDigis     = std::make_unique<CaloRecoDigiCollection>();
      auto caloHits      = std::make_unique<CaloHitCollection>();
      const art::ProductID         recoDigisPID    = event.getProductID<CaloRecoDigiCollection>();
      const art::EDProductGetter*  recoDigisGetter = event.productGetter(recoDigisPID);

      const Calorimeter& cal = *GeomHandle<Calorimeter>();
      std::set<size_t> keep;
      for (unsigned int g = 0; g < groups.size(); ++g)
      {
          const auto& group = groups[g];
          std::map<int, double> eSiPM;                             // summed amplitude per SiPM in this time group
          for (const auto& hit : group) eSiPM[hit.sipmID] += hit.eDep;

          std::map<int, double> eCrystal;                          // crystal energy: average over its SiPMs that have a peak
          std::map<int, int>    nSiPM;
          for (const auto& [sipmID, e] : eSiPM) { eCrystal[sipmID/2] += e; ++nSiPM[sipmID/2]; }
          for (auto& [cryID, e] : eCrystal) e /= nSiPM[cryID];

          for (const auto& [cryID, eCry] : eCrystal)
          {
              ++nCrystalsIn_;
              double eNeigh = 0.0;                                 // summed energy of the up-to-8 surrounding crystals
              for (int n : cal.crystal(cryID).neighbors())
              {
                  auto it = eCrystal.find(n);                      // only neighbours that have a digi in this group
                  if (it != eCrystal.end()) eNeigh += it->second;
              }
              bool passes = eCry > 0 && eCry/(eCry + eNeigh) >= ratioCut_;
              if (diagLevel_ > 0)
                  std::cout << "[SourceCalibSelector] event " << event.event() << " group " << g
                            << " crystal " << cryID << " eCry " << eCry << " eNeigh " << eNeigh
                            << " frac " << eCry/(eCry + eNeigh) << " pass " << passes << "\n";
              if (!passes) continue;

              // this crystal's peaks as CaloRecoDigis (ADC), and one CaloHit pointing to them
              std::vector<art::Ptr<CaloRecoDigi>> hitRecoDigis;
              double timeSum(0), timeErr2(0), eDepErr2(0);
              for (const auto& hit : group)
              {
                  if (hit.crystalID != cryID) continue;
                  keep.insert(hit.digiIdx);
                  recoDigis->emplace_back(art::Ptr<CaloDigi>(caloDigisH, hit.digiIdx), hit.eDep, hit.eDepErr,
                                          hit.peakTime, hit.peakTimeErr, hit.chi2, hit.ndf, hit.pileUp);
                  hitRecoDigis.emplace_back(recoDigisPID, recoDigis->size()-1, recoDigisGetter);
                  timeSum  += hit.peakTime;
                  timeErr2 += hit.peakTimeErr*hit.peakTimeErr;
                  eDepErr2 += hit.eDepErr*hit.eDepErr;
              }
              int    nRO   = nSiPM[cryID];
              double nPeak = hitRecoDigis.size();
              caloHits->emplace_back(cryID, nRO, timeSum/nPeak, std::sqrt(timeErr2)/nPeak,
                                     eCry, std::sqrt(eDepErr2)/nRO, hitRecoDigis);
          }
      }

      // 4. selected digis for EventNtuple, in their original order
      for (size_t idx : keep) selectedDigis->push_back(caloDigis[idx]);

      bool pass = !caloHits->empty();
      ++nEvents_;
      if (pass) ++nEventsPassed_;
      nCrystalsOut_ += caloHits->size();
      nDigisIn_     += caloDigis.size();
      nDigisOut_    += selectedDigis->size();

      event.put(std::move(selectedDigis));
      event.put(std::move(recoDigis));
      event.put(std::move(caloHits));
      return pass;
  }

  void SourceCalibSelector::endJob()
  {
      std::cout << "[SourceCalibSelector] " << nEventsPassed_ << " of " << nEvents_ << " events passed; selected "
                << nCrystalsOut_ << " of " << nCrystalsIn_ << " crystal entries (crystal per time group); kept "
                << nDigisOut_ << " of " << nDigisIn_ << " CaloDigis" << std::endl;

      // cut flow, in the order the cuts are applied: what is left after each cut, what that cut removed,
      // removed as % of the line above, and what is left as % of the first line
      struct Step { std::string name; unsigned long n; };
      const std::vector<std::pair<std::string, std::vector<Step>>> flows{
          {"CaloDigis (one per SiPM pulse window)", {
              {"in",                                                      nDigisIn_},
              {"have a fitted peak (3 samples >= minPeakAmplitude)",      nDigisIn_ - nDigisNoPeak_},
              {"good fit (chi2/ndf <= maxChi2Cut)",                       nDigisFitted_},
              {"their crystal passes isolation (kept)",                   nDigisOut_}}},
          {"crystal entries (one crystal in one time group)", {
              {"with a fitted peak",                                      nCrystalsIn_},
              {"pass isolation (selected, = SourceAna entries)",          nCrystalsOut_}}},
          {"events", {
              {"in",                                                      nEvents_},
              {"with at least one selected crystal",                      nEventsPassed_}}}};

      printf("[SourceCalibSelector] cutflow  (maxChi2Cut %g, ratioCut %g, timeCut %g ns)\n", maxChi2Cut_, ratioCut_, timeCut_);
      for (const auto& [title, steps] : flows)
      {
          printf("[SourceCalibSelector] cutflow  %s\n", title.c_str());
          printf("[SourceCalibSelector] cutflow    %-52s %12s %12s %9s %9s\n", "step", "left", "removed", "%removed", "%left");
          for (size_t i = 0; i < steps.size(); ++i)
          {
              unsigned long prev    = (i == 0) ? steps[0].n : steps[i-1].n;
              unsigned long removed = prev - steps[i].n;
              printf("[SourceCalibSelector] cutflow    %-52s %12lu %12lu %8.1f%% %8.1f%%\n", steps[i].name.c_str(),
                     steps[i].n, removed, prev ? 100.0*removed/prev : 0.0, steps[0].n ? 100.0*steps[i].n/steps[0].n : 0.0);
          }
      }
  }
}

DEFINE_ART_MODULE(mu2e::SourceCalibSelector)
