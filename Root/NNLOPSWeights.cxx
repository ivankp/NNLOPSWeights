#include "NNLOPSWeights/NNLOPSWeights.h"

#include <iostream>
#include <memory>
#include <type_traits>

#include <TLorentzVector.h>

// #include <EventLoop/Worker.h>
// #include <xAODTruth/TruthEventContainer.h>
#include <HGamAnalysisFramework/HgammaIncludes.h>
// #include <TruthWeightTools/HiggsWeightTool.h>
// #include <HepMC/GenEvent.h>

// this is needed to distribute the algorithm to the workers
ClassImp(NNLOPSWeights)

using std::cout;
using std::endl;

struct NNLOPSWeights::impl {

  // Tools
  // xAODtoHepMCTool *m_hepmc_tool;
  // HiggsTruthCategoryTool *m_htxs_tool;
  // xAOD::HiggsWeightTool *m_higgsMCtool;
  // xAOD::HiggsWeightTool *weightTool;

  // Histograms
  TH1 *h_pT_yy;

  impl()
  // : m_hepmc_tool(new xAODtoHepMCTool("xAODtoHepMCTool")),
  //   m_htxs_tool(new HiggsTruthCategoryTool("HiggsTruthCategoryTool")),
  //   m_higgsMCtool(new xAOD::HiggsWeightTool("HiggsWeightTool"))
  // : weightTool(new xAOD::HiggsWeightTool("HiggsWeightTool"))
  {
    // weightTool->setProperty("ForceNNLOPS",true);
    // weightTool->initialize();

    const double b_pT_yy[] = { 0,20,30,45,60,80,120,170,220,350 };
    h_pT_yy = new TH1D("pT_yy","",std::extent<decltype(b_pT_yy)>::value-1,b_pT_yy);
  }
  ~impl() {
    // delete m_higgsMCtool;
    // delete weightTool;
    // delete h_pT_yy;
  }
};

NNLOPSWeights::NNLOPSWeights(const char *name): HgammaAnalysis(name) { }
NNLOPSWeights::~NNLOPSWeights() {
  // if (p) {
  //   delete p;
  //   p = nullptr;
  // }
}

EL::StatusCode NNLOPSWeights::createOutput() {
  // p = new impl;
  p.reset(new impl);
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode NNLOPSWeights::initialize() {
  const EL::StatusCode code = HgammaAnalysis::initialize();
  return code;
}

EL::StatusCode NNLOPSWeights::execute() {
  HgammaAnalysis::execute();

  // std::vector<std::string> weightNames = p->weightTool->getWeightNames();
  // for (const auto& name : weightNames)
  //   cout << name << endl;

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
    // xAOD::HiggsWeights higgsWeights = truth->weights();

    // const xAOD::TruthEventContainer *truthEvents = nullptr;
    // const std::vector<float>& weights = truthEvents->at(0)->weights();

    // const xAOD::TruthEventContainer *truthEvents = nullptr;
    // wk()->xaodEvent()->retrieve(truthEvents,"TruthEvents");

    // std::vector<HepMC::GenEvent> hepmcEvents =
    //   p->m_hepmc_tool->getHepMCEvents( truthEvents, eventInfo );
    // HTXS::HiggsClassification htxs =
    //   p->m_htxs_tool->getHiggsTruthCategoryObject( hepmcEvents[0], m_prodMode );
    // xAOD::HiggsWeights higgsWeights =
    //   p->m_higgsMCtool->getHiggsWeights( htxs.jets30.size(), htxs.higgs.Pt() );

    // const xAOD::HiggsWeights w = p->weightTool->getHiggsWeights(
    //   );
    // --------------------------------------------------------------

    // Compute variables --------------------------------------------
    const TLorentzVector yy = photons[0]->p4() + photons[1]->p4();
    // --------------------------------------------------------------

    // Fill histograms ----------------------------------------------
    p->h_pT_yy->Fill(yy.Pt(),w.nominal);
    // --------------------------------------------------------------
  }

  return EL::StatusCode::SUCCESS;
}

