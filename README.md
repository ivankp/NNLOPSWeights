# Setup

```
setupATLAS
rcSetup Base,2.4.28
rc checkout_pkg atlasoff/PhysicsAnalysis/HiggsPhys/Run2/HGamma/xAOD/HGamAnalysisFramework/tags/HGamAnalysisFramework-00-02-77-12
./HGamAnalysisFramework/scripts/setupRelease
rc find_packages
rc compile
```

# Run

```
runNNLOPSWeights data/test.cfg
```
