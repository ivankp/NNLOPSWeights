#include "NNLOPSWeights/NNLOPSWeights.h"

#include <type_traits>

#include <TParameter.h>

#include <HGamAnalysisFramework/HgammaIncludes.h>
#include <HGamAnalysisFramework/HGamVariables.h> // for var::

template <typename T>
using decay_t = typename std::decay<T>::type;

// this is needed to distribute the algorithm to the workers
ClassImp(NNLOPSWeights)

using std::cout;
using std::endl;

struct NNLOPSWeights::impl {
  TTree *tree;
  TParameter<double> xs_br_fe;

#ifdef VAR
#error "Macro name VAR already in use"
#endif
#define VAR(NAME) decay_t<decltype(var::NAME.truth())> NAME;
  VAR(N_j_30)
  VAR(pT_yy)
  VAR(pT_j1_30)
  VAR(m_jj_30)
  VAR(Dphi_j_j_30)
  VAR(Dphi_j_j_30_signed)
#undef VAR

  xAOD::HiggsWeights hw;

  impl(): tree(new TTree("tree","")), xs_br_fe("crossSectionBRfilterEff",0) {
#define VAR(NAME) tree->Branch(#NAME,&NAME);
    VAR(N_j_30)
    VAR(pT_yy)
    VAR(pT_j1_30)
    VAR(m_jj_30)
    VAR(Dphi_j_j_30)
    VAR(Dphi_j_j_30_signed)
#undef VAR
    tree->Branch("w_nominal",&hw.nominal);
    tree->Branch("w_pdf4lhc_unc",&hw.pdf4lhc_unc);
    tree->Branch("w_nnpdf30_unc",&hw.nnpdf30_unc);
    tree->Branch("w_qcd",&hw.qcd);
    tree->Branch("w_qcd_nnlops",&hw.qcd_nnlops);
  }
};

NNLOPSWeights::NNLOPSWeights(const char *name): HgammaAnalysis(name) { }
NNLOPSWeights::~NNLOPSWeights() { }

EL::StatusCode NNLOPSWeights::createOutput() {
  p.reset(new impl);
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode NNLOPSWeights::initialize() {
  const EL::StatusCode status = HgammaAnalysis::initialize();
  wk()->addOutput(p->tree);
  p->xs_br_fe.SetVal( HgammaAnalysis::crossSectionBRfilterEff() );
  wk()->addOutput(&p->xs_br_fe);
  return status;
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

    HG::VarHandler::getInstance()->setTruthContainers(
      &all_photons, &electrons, &muons, &jets);
    // --------------------------------------------------------------

    // Compute variables --------------------------------------------
    // NOTE: var is a namespace
#define VAR(NAME) p->NAME = var::NAME.truth();
    VAR(N_j_30)
    VAR(pT_yy)
    VAR(pT_j1_30)
    VAR(m_jj_30)
    VAR(Dphi_j_j_30)
    VAR(Dphi_j_j_30_signed)
#undef VAR
    // --------------------------------------------------------------

    // Obtain weights -----------------------------------------------
    p->hw = event->higgsWeights();
    // --------------------------------------------------------------

    p->tree->Fill();
  }

  return EL::StatusCode::SUCCESS;
}

