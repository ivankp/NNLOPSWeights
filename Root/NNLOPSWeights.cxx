#include "NNLOPSWeights/NNLOPSWeights.h"
#include <HGamAnalysisFramework/HgammaIncludes.h>

#include <iostream>
#include <memory>
#include <TLorentzVector.h>

// this is needed to distribute the algorithm to the workers
ClassImp(NNLOPSWeights)

struct NNLOPSWeights::impl {
  TH1 *h_pT_yy;
  impl() {
    const double b_pT_yy[] = { 0,20,30,45,60,80,120,170,220,350 };
    h_pT_yy = new TH1D("pT_yy","",9,b_pT_yy);
  }
  ~impl() { delete h_pT_yy; }
};

NNLOPSWeights::NNLOPSWeights(const char *name): HgammaAnalysis(name) { }
NNLOPSWeights::~NNLOPSWeights() {
  delete p;
  p = nullptr;
}

EL::StatusCode NNLOPSWeights::createOutput() {
  p = new impl;
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode NNLOPSWeights::initialize() {
  const EL::StatusCode code = HgammaAnalysis::initialize();
  return code;
}

EL::StatusCode NNLOPSWeights::execute() {
  HgammaAnalysis::execute();
  auto* event = eventHandler();
  auto* truth = truthHandler();

  auto photons1 = truth->getPhotons();
  const auto photons = truth->applyPhotonSelection(photons1);

  // Must have at least 2 truth photons
  if (photons.size() < 2) return EL::StatusCode::SUCCESS;

  // Truth diphoton
  const TLorentzVector yy = photons[0]->p4() + photons[1]->p4();
  const double m_yy = yy.M();
  if (m_yy<105. || 160.<m_yy) return EL::StatusCode::SUCCESS;

  auto electrons1 = truth->getElectrons();
  auto muons1 = truth->getMuons();
  auto jets1 = truth->getJets();
  const auto electrons = truth->applyElectronSelection(electrons1);
  const auto muons = truth->applyMuonSelection(muons1);
  const auto jets = truth->applyJetSelection(jets1);

  if (truth->passFiducial(&photons,&electrons,&muons,&jets)) {
    xAOD::HiggsWeights higgsWeights = event->higgsWeights();

    const double nominal = higgsWeights.nominal;

    p->h_pT_yy->Fill(yy.Pt(),nominal);
  }

  return EL::StatusCode::SUCCESS;
}

