#ifndef NNLOPSWeights_H
#define NNLOPSWeights_H

#include "HGamAnalysisFramework/HgammaAnalysis.h"

class NNLOPSWeights : public HgammaAnalysis {
private:
  class impl;
  std::shared_ptr<impl> p; //!

public:
  NNLOPSWeights() { }
  NNLOPSWeights(const char* name);
  virtual ~NNLOPSWeights();

  // these are the functions inherited from HgammaAnalysis
  virtual EL::StatusCode createOutput();
  virtual EL::StatusCode execute();
  virtual EL::StatusCode initialize();

  // this is needed to distribute the algorithm to the workers
  ClassDef(NNLOPSWeights, 1);
};

#endif
