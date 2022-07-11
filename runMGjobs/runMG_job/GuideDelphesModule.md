# Guide of adding modules on Delphes 
It is easy to modify and add modules to Delphes, this guide shows how to
1) add jets modules with different algorithm and parameters,
2) add B-tagging algorithm,
3) add Tau-tagging algorithm.
## Jet algorithm
Reconstruction Jet algorithms could be easily apply using the FastJetFinder module with following syntax:
```
module FastJetFinder <ModuleNameForJetsAlgo> {
    
    set InputArray EFlowMerger/eflow
    set OuputArray <jetName>
    
    set JetAlgorithm <a number>
    # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt, 7 anti-kt with winner-take-all axis (for N-subjettiness), 8 N-jettiness, 9 Valencia
    
    Set ParameterR <cone size>
    set JetPTMin <jet's minimum transverse momentum cut>

}
```

Generator level jet algorithms is very similar, you only need to change the InputArray from "EFlowMerger/eflow" to "NeutrinoFilter/filteredParticles".
## B-tagging algorithm
B-tagging module is dependent on the Jet Flavor Association modules. To implementing Jet Flavor Association module, adding module to "MuonColliderDet_JetFlavorAssociation.tcl" in the card/MuonCollider directory with following syntax:
```
module JetFlavorAssociation <ModuleNameForThis> {
    
    set PartonInputArray Delphes/partons
    set ParticleInputArray Delphes/allParticles
    set ParticleLHEFInputArray Delphes/allParticlesLHEF
    set JetInputArray <ModuleNameForJetsAlgo>/<jetName>

    set DeltaR 0.5 
    set PartonPTMin 1.0
    set PartonEtaMax 2.5
    # these parameters are for Muon Collider

}
```
Then, we can add the B-tagging algorithms for different Working Points to "MuonColliderDet_BTagging.tcl" in the same directory with following syntax:
```
module BTagging <ModuleNameForBTaggingWP50/70/90> {
    
    set JetInputArray <ModuleNameForJetsAlgo>/<jetName>
    setBitNumber <0forWP50,1forWP70,2forWP90>
    source MuonCollider/MuonColliderDet_BTag_<50/70/90>.tcl

}
```
Notice that here we did not apply Jet Smearing algorithm, if needed, change JetInputArray to "\<JetMomentumSmearingModuleName\>/\<jetName\>".
## Tau-tagging algorithm
Similar to B-tagging, Tau-tagging also requires Jet Flavor Association. Then we can add the module to "MuonColliderDet_TauTagging.tcl" with following syntax:
```
module Tau-Tagging <ModuleNameForTauTagging> {

    set ParticleInputArray Delphes/allParticles
    set PartonInputArray Delphes/partons
    set JetInputArray <ModuleNameForJetsAlgo>/<jetName>
    set DeltaR 0.5
    set TauPTMin 1.0
    set TauEtaMax 2.5
    add EfficiencyFormula {0} {0.02} #consider set this as zero, as 2% false rate is too high
    add EfficiencyFormula {11} {0.001}
    add EfficiencyFormula {15} {
        (pt < 10) * (0.0) +
	(pt >=10) * (0.80)
    }
}
```
## Execution Path of Delphes  
New modules need to be add into the execution path of Delphes to be able to run:
```
set ExecutionPath {

    <Modules already exist in your main Delphes card>
    <Your new modules in correct order>
    <Modules already exist in your main Delphes card>

}
```
Information you would like to be save into the tree also need to be added, take jet as an example:
```
module TreeWriter TreeWriter {

    ...
    add Branch <JetModuleName>/<OutputArrayName> <NameYouWantJetToBeShowInTree> Jet
    ...

}
```
