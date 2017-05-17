#include "NNLOPSWeights/NNLOPSWeights.h"

#include <iostream>
#include <memory>
#include <type_traits>

#include <TLorentzVector.h>

#include <HGamAnalysisFramework/HgammaIncludes.h>

namespace detail {
  template <typename S, typename T>
  inline void cat_impl(S& s, const T& t) { s << t; }
  template <typename S, typename T, typename... TT>
  inline void cat_impl(S& s, const T& t, const TT&... tt) {
    s << t;
    cat_impl(s,tt...);
  }
}
template <typename... TT>
inline std::string cat(const TT&... tt) {
  std::stringstream ss;
  detail::cat_impl(ss,tt...);
  return ss.str();
}

// this is needed to distribute the algorithm to the workers
ClassImp(NNLOPSWeights)

using std::cout;
using std::endl;

class hist_bundle {
  static constexpr size_t bundle_size = 1;
  static const std::array<std::string,bundle_size> wnames;
  static std::vector<const hist_bundle*> all;

  std::array<TH1D*,bundle_size> hh;

public:
  static std::array<double,bundle_size> weights;

  template <size_t N>
  hist_bundle(const char* name, const double(&bins)[N]) {
    all.push_back(this);
    for (size_t i=0; i<bundle_size; ++i)
      hh[i] = new TH1D( cat(name,':',wnames[i]).c_str(), "", N-1, bins );
  }
  void operator()(double x) {
    for (unsigned i=0, n=hh.size(); i<n; ++i)
      hh[i]->Fill(x,weights[i]);
  }
  inline TH1D* operator[](size_t i) const noexcept { return hh[i]; }

  template <typename WK>
  static void add_all(WK* wk) {
    for (const hist_bundle* hb : all)
      for (TH1D* h : hb->hh)
        wk->addOutput(h);
  }
};
decltype(hist_bundle::wnames) hist_bundle::wnames {
  "nominal"
};
decltype(hist_bundle::weights) hist_bundle::weights;
decltype(hist_bundle::all) hist_bundle::all;

constexpr double b_pT_yy[] = { 0,20,30,45,60,80,120,170,220,350 };

struct NNLOPSWeights::impl {
  hist_bundle h_pT_yy;

  impl()
  : h_pT_yy("pT_yy",b_pT_yy)
  {
    // const double b_pT_yy[] = { 0,20,30,45,60,80,120,170,220,350 };
    // h_pT_yy = new TH1D("pT_yy","",std::extent<decltype(b_pT_yy)>::value-1,b_pT_yy);

  }
};

NNLOPSWeights::NNLOPSWeights(const char *name): HgammaAnalysis(name) { }
NNLOPSWeights::~NNLOPSWeights() { }

EL::StatusCode NNLOPSWeights::createOutput() {
  p.reset(new impl);
  // p = new impl;
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode NNLOPSWeights::initialize() {
  const EL::StatusCode code = HgammaAnalysis::initialize();
  hist_bundle::add_all(wk());
  return code;
}

EL::StatusCode NNLOPSWeights::execute() {
  HgammaAnalysis::execute();

  auto* event = eventHandler();
  auto* truth = truthHandler();

  // select photons -------------------------------------------------
  auto all_photons = truth->getPhotons();
  auto photons = truth->applyPhotonSelection(all_photons);
  // ----------------------------------------------------------------

  if (truth->passFiducial(&photons)) {
    // select other objects -----------------------------------------
    auto all_electrons = truth->getElectrons();
    auto all_muons = truth->getMuons();
    auto all_jets = truth->getJets();
    auto electrons = truth->applyElectronSelection(all_electrons);
    auto muons = truth->applyMuonSelection(all_muons);
    auto jets = truth->applyJetSelection(all_jets);

    truth->removeOverlap(photons,jets,electrons,muons);
    // --------------------------------------------------------------

    // Obtain weights -----------------------------------------------
    xAOD::HiggsWeights w = event->higgsWeights();
    // --------------------------------------------------------------

    // Compute variables --------------------------------------------
    const TLorentzVector yy = photons[0]->p4() + photons[1]->p4();
    // --------------------------------------------------------------

    hist_bundle::weights[0] = w.nominal;

    // Fill histograms ----------------------------------------------
    // p->h_pT_yy->Fill(yy.Pt(),w.nominal);
    p->h_pT_yy(yy.Pt()*1e-3);
    // --------------------------------------------------------------
  }

  return EL::StatusCode::SUCCESS;
}

