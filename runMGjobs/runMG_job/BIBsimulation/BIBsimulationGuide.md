# Guide of BIB simulation in Delphes
This is two version of BIB simulation in Delphes, this guide shows how to use them and how they work.
## Direct BIB simultion
One of the most direct way to simulate BIB in Delphes is to directly generate randomized BIB particles according to given distribution in a root file and add them at the very begining. In this way you need to have very large memory space for large number of BIB particle.
The module have following syntax:
```
module BIBModule BIBModule {
    set InputArray Delphes/allParticles
    set StableInputArray Delphes/stableParticles
	    
    set OutArray allParticles
    set StableOutArray stableParticles

    set NumParticles 10000 

    set FileName /Delphes/cards/MuonCollider/histograms_MCPar_muComb_files1To8_Allevts_new.root
    set xHistName x
    set PositionHistName z_r
    set MomentumHistName px_py_pz
    set ThetaHistName theta 
    set PdgIDHistName pdgid 
    set PhiHistName phi 
}
```
An exmaple script is delphes_card_MuonColliderDet_HHstudy_BIB.tcl.
## BIB simulation in proximity
In order to avoid very large memory useage, we could only considered BIB particles in proximity of the original particles. We have design a way such that the BIB particles is adding in calorimeter in two step: ECal and HCal. The modules have following syntax:
```
module BIBECal BIBECal {
    set PhotonsInputArray ECal/eflowPhotons
    set TracksInputArray ECal/eflowTracks

    set PhotonsOutputArray eflowPhotons
    set TracksOutputArray eflowTracks

    set NumParticles 10000000    

    set FileName /Delphes/cards/MuonCollider/histograms_MCPar_muComb_files1To8_Allevts_new.root
    set xHistName x
    set PositionHistName z_r
    set MomentumHistName px_py_pz
    set PdgEnergyHistName pdgid_energy 
    set PhotonsDeltaR 0.05 
    set TracksDeltaR 0.07
    set Bz 4.0
}
#############
#   HCAL
#############
module BIBNeutralHadrons BIBNeutralHadrons {
    set InputArray BIBHCal/eflowNeutralHadrons

    set OutputArray eflowNeutralHadrons

    set NumParticles 10000000 

    set FileName /Delphes/cards/MuonCollider/histograms_MCPar_muComb_files1To8_Allevts_new.root
    set xHistName x
    set PositionHistName z_r
    set MomentumHistName px_py_pz
    set PdgEnergyHistName pdgid_energy 
    set DeltaR 0.70

}
```
An exmaple script is delphes_card_MuonColliderDetCaloBIB.tcl. 
## Debugging and other Issues
The scripts are written specifically for Muon Collider simulation. Please feel free to contact Kenny Jia through email: hjia38@wisc.edu, hjia625@stanford.edu, or haoyi.jia@cern.ch. I would be happy to help with the debugging process or questions!
