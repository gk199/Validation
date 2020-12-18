#import ROOT
from ROOT import *
from DataFormats.FWLite import Events, Handle
import numpy as np
#import tifffile
import math
import intersect

def etaVal(ieta):
  if (ieta <= -24):
    etavl = .1695*ieta + 1.9931
  elif (ieta <= -1):
    etavl = .0875*ieta + .0489
  elif (ieta < 24):
    etavl = .0875*ieta - .0489
  else:
    etavl = .1695*ieta - 1.9931
  return etavl

def phiVal(iphi):
  phiBins=72.
  phivl=float(iphi)*(2.*math.pi/phiBins);
  if (iphi > 36):
    phivl -= 2.*math.pi;
  return phivl

gROOT.SetBatch(True)

#f = TFile.Open("SinglePiPt100RECO.root")
#f = TFile.Open("SinglePiPt50ieta23RECO.root")
#f = TFile.Open("SinglePiPt1-5ieta22-24RECO.root")
#f = TFile.Open("SinglePiPt5-15ieta22-24RECO.root")
#f = TFile.Open("MultiPartPt50ieta22-24RECO.root")
#f = TFile.Open("SingleK0Pt50ieta23RECO.root")
#f = TFile.Open("SingleN0Pt50ieta23RECO.root")
#f = TFile.Open("output_new_tau0_step0_1000.root")
#f = TFile.Open("output_new_step0_1000.root")
#f = TFile.Open("output_new_step0_100.root")
#f = TFile.Open("output_decay_tau0_step0_1000.root")
#f = TFile.Open("output_tau0_step0_1000.root")
#f = TFile.Open("output_step0_1000.root")
#f = TFile.Open("output_step0_100.root")
#f = TFile.Open("output_step0.root")
#f = TFile.Open("step1tdc_prompt.root")
#f = TFile.Open("/eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/mh2000_mx975_pl10000_ev1000/ppTohToSS1SS2_SS1Tobb_SS2Toveve_1_withISR_step1_TDC.root")
#f = TFile.Open("/eos/cms/store/user/lowang/LLP_htobbbb/step1/MH-125_MFF-50_CTau-1000mm_step1.root")
#f = TFile.Open("/eos/cms/store/user/lowang/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_step1.root")
#f = TFile.Open("/eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/PionGun/SinglePion211_E10_PU_00_step1.root")
#f = TFile.Open("/eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/PionGun/SinglePion211_E10_PU_step1.root")
f = TFile.Open("/eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/PionGun/SinglePion211_E10_PU_00_eta1phi0_timeslew-false_step1_CaloSamples_10events.root") #SinglePion211_E10_PU_00_eta1phi0_timeslew-false_step1.root")
#f = TFile.Open("step1tdc.root")

events = f.Events

eclust = TH1D('eclust','Cluster Energy',50,0.0,400.0)
tclust = TH1D('tclust','Cluster Time',50,-1.0,10.0)
lclust = TH1D('lclust','Cluster Layer',10,-0.5,9.5)
dclust = TH1D('dclust','Cluster Depth',10,-0.5,9.5)
eclust2 = TH1D('eclust2','HBHE Cluster Energy',50,0.0,400.0)
tclust2 = TH1D('tclust2','HBHE Cluster Time',50,-1.0,10.0)
lclust2 = TH1D('lclust2','HBHE Cluster Layer',10,-0.5,9.5)
dclust2 = TH1D('dclust2','HBHE Cluster Depth',10,-0.5,9.5)
ehit = TH1D('ehit','Hit Energy',50,0.0,10.0)
thit = TH1D('thit','Hit Time',50,-1.0,10.0)
dhit = TH1D('dhit','Hit Depth',10,-0.5,9.5)
phit = TH1D('phit','Hit Phi',74,-0.5,73.5)
nhit = TH1D('nhit','Hit Eta',31,-0.5,30.5)
ehit2 = TH1D('ehit2','HBHE Hit Energy',50,0.0,10.0)
thit2 = TH1D('thit2','HBHE Hit Time',50,-1.0,10.0)
lhit2 = TH1D('lhit2','HBHE Hit Layer',10,-0.5,9.5)
dhit2 = TH1D('dhit2','HBHE Hit Depth',10,-0.5,9.5)
#ehit3 = TH1D('ehit3','HBHE SimHit Energy',50,0.0,0.1)
ehit3 = TH1D('ehit3','HBHE SimHit Energy',50,0.0,10.0)
ehit4 = TH1D('ehit4','HB SimHit CaloRegion Energy',100,0.0,250.0)
thit3 = TH1D('thit3','HBHE SimHit Time',502,-1.0,501.0)
dhit3 = TH1D('dhit3','HBHE SimHit Depth',20,-0.5,19.5)
phit3 = TH1D('phit3','HBHE SimHit Phi',74,-0.5,73.5)
nhit3 = TH1D('nhit3','HBHE SimHit Eta',31,-0.5,30.5)
thit4 = TH1D('thit4','HBHE SimHit Time',502,-1.0,21.0)
vhit4 = TH1D('vhit4','HBHE SimHit Vertex Distance',121,0.5,60.5)
vhit5 = TH1D('vhit5','HBHE SimHit Vertex Distance',121,0.5,60.5)
chit1 = TH1D('chit1','Cosine(theta) of Gen Momentum and Displaced Vertex Direction',100,-1.0,1.0)
chit2d = TH2D('chit2d','delayed time versus 1/sin(theta) ',100,1,5,502,-1.0,51.0)
chit3 = TH1D('chit3','Cosine(theta) of Gen Momentum and Displaced Vertex Direction',100,-1.0,1.0)
qhit1 = TH1D('qhit1','Path Length Difference Time Delay',24,-1.0,5.0)
qhit2d = TH2D('qhit2d','Delayed Time versis Path Length Difference Time Delay',24,-1.0,5.0,50,-1.0,24.0)
qhit3 = TH1D('qhit3','Path Length Difference Time Delay',24,-1.0,5.0)
qhit4d = TH2D('qhit4d','Delayed Time versis Path Length Difference Time Delay',24,-1.0,5.0,50,-1.0,24.0)
qhit5d = TH2D('qhit5d','Delayed Time versis CaloRegion Energy',25,0.0,250.0,50,-1.0,24.0)
qhit6d = TH2D('qhit6d','CaloRegion Delayed Time versis Path Length Difference Time Delay',24,-1.0,5.0,50,-1.0,24.0)
thit5 = TH1D('thit5','HB SimHit Delayed Time (Delta R<0.5)',62,-1.0,31.0)
thit6 = TH1D('thit6','HB SimHit Delayed Time',62,-1.0,31.0)
thit7 = TH1D('thit7','HB SimHit CaloRegion Time Delay',100,-1.0,24.0)
thit8 = TH1D('thit8','HB SimHit Time',502,-1.0,21.0)
thit9 = TH1D('thit9','HB SimHit Time',502,-1.0,21.0)
mhit1 = TH1D('mhit1','Multiplicity of 1ns delayed cells above 3 GeV',140,-0.5,139.5)
mhit3 = TH1D('mhit3','Multiplicity of 2ns delayed cells above 3 GeV',140,-0.5,139.5)
mhit5 = TH1D('mhit5','Multiplicity of 1ns delayed cells above 1000 fC',140,-0.5,139.5)
mhit7 = TH1D('mhit7','Multiplicity of 2ns delayed cells above 1000 fC',140,-0.5,139.5)
adchit1 = TH1D('adchit1','HB ADC',256,-0.5,255.5)
tdchit1 = TH1D('tdchit1','HB TDC',64,-0.5,63.5)
tdchit2d = TH2D('tdchit2d','HB TDC versus TOF',46,3.5,49.5,40,-1.5,18.5)
tdchit3 = TH1D('tdchit3','HB TDC',64,-0.5,63.5)
tdchit4 = TH1D('tdchit4','HB TDC(ns)-TOF (Depth 1), 1k fC',100,-25.0,25.0)
tdchit5 = TH1D('tdchit5','HB TDC(ns)-TOF (Depth 2), 1k fC',100,-25.0,25.0)
tdchit6 = TH1D('tdchit6','HB TDC(ns)-TOF (Depth 3), 1k fC',100,-25.0,25.0)
tdchit7 = TH1D('tdchit7','HB TDC(ns)-TOF (Depth 4), 1k fC',100,-25.0,25.0)
tdchit4HE = TH1D('tdchit4HE','HE TDC(ns)-TOF (Depth 1), 1k fC',100,-25.0,25.0)
tdchit5HE = TH1D('tdchit5HE','HE TDC(ns)-TOF (Depth 2), 1k fC',100,-25.0,25.0)
tdchit6HE = TH1D('tdchit6HE','HE TDC(ns)-TOF (Depth 3), 1k fC',100,-25.0,25.0)
tdchit7HE = TH1D('tdchit7HE','HE TDC(ns)-TOF (Depth 4), 1k fC',100,-25.0,25.0)
tdchit8HE = TH1D('tdchit8HE','HE TDC(ns)-TOF (Depth 5), 1k fC',100,-25.0,25.0)
tdchit9HE = TH1D('tdchit9HE','HE TDC(ns)-TOF (Depth 6+), 1k fC',100,-25.0,25.0)
tdchit4_3kfC = TH1D('tdchit4_3kfC','HB TDC(ns)-TOF (Depth 1), 3k fC',100,-25.0,25.0)
tdchit5_3kfC = TH1D('tdchit5_3kfC','HB TDC(ns)-TOF (Depth 2), 3k fC',100,-25.0,25.0)
tdchit6_3kfC = TH1D('tdchit6_3kfC','HB TDC(ns)-TOF (Depth 3), 3k fC',100,-25.0,25.0)
tdchit7_3kfC = TH1D('tdchit7_3kfC','HB TDC(ns)-TOF (Depth 4), 3k fC',100,-25.0,25.0)
tdchit4HE_3kfC = TH1D('tdchit4HE_3kfC','HE TDC(ns)-TOF (Depth 1), 3k fC',100,-25.0,25.0)
tdchit5HE_3kfC = TH1D('tdchit5HE_3kfC','HE TDC(ns)-TOF (Depth 2), 3k fC',100,-25.0,25.0)
tdchit6HE_3kfC = TH1D('tdchit6HE_3kfC','HE TDC(ns)-TOF (Depth 3), 3k fC',100,-25.0,25.0)
tdchit7HE_3kfC = TH1D('tdchit7HE_3kfC','HE TDC(ns)-TOF (Depth 4), 3k fC',100,-25.0,25.0)
tdchit8HE_3kfC = TH1D('tdchit8HE_3kfC','HE TDC(ns)-TOF (Depth 5), 3k fC',100,-25.0,25.0)
tdchit9HE_3kfC = TH1D('tdchit9HE_3kfC','HE TDC(ns)-TOF (Depth 6+), 3k fC',100,-25.0,25.0)
fhit1 = TH1D('fhit1','HB fC',50,-50.0,1500.0)
vt = TH1D('vt','v t',50,-1.0,1.0)
vx = TH1D('vx','v x',50,-0.01,0.01)
vy = TH1D('vy','v y',50,-0.01,0.01)
vz = TH1D('vz','v z',50,-15.0,15.0)


#floatROOTMathCartesian3DROOTMathDefaultCoordinateSystemTagROOTMathPositionVector3D_genParticles_xyz0_HLT
events.Draw("float_genParticles_t0_HLT.obj>>vt", "","goff")
events.Draw("floatROOTMathCartesian3DROOTMathDefaultCoordinateSystemTagROOTMathPositionVector3D_genParticles_xyz0_HLT.obj.x()>>vx","","goff")
events.Draw("floatROOTMathCartesian3DROOTMathDefaultCoordinateSystemTagROOTMathPositionVector3D_genParticles_xyz0_HLT.obj.y()>>vy","","goff")
events.Draw("floatROOTMathCartesian3DROOTMathDefaultCoordinateSystemTagROOTMathPositionVector3D_genParticles_xyz0_HLT.obj.z()>>vz","","goff")
#events.Draw("float_genParticles_t0_SIM.obj>>vt", "","goff")
#events.Draw("floatROOTMathCartesian3DROOTMathDefaultCoordinateSystemTagROOTMathPositionVector3D_genParticles_xyz0_SIM.obj.x()>>vx","","goff")
#events.Draw("floatROOTMathCartesian3DROOTMathDefaultCoordinateSystemTagROOTMathPositionVector3D_genParticles_xyz0_SIM.obj.y()>>vy","","goff")
#events.Draw("floatROOTMathCartesian3DROOTMathDefaultCoordinateSystemTagROOTMathPositionVector3D_genParticles_xyz0_SIM.obj.z()>>vz","","goff")
#events.Draw("recoPFClusters_particleFlowClusterHCAL__RECO.obj.energy()>>eclust","","goff")
#events.Draw("recoPFClusters_particleFlowClusterHCAL__RECO.obj.time()>>tclust","","goff")
#events.Draw("recoPFClusters_particleFlowClusterHCAL__RECO.obj.layer()>>lclust","","goff")
#events.Draw("recoPFClusters_particleFlowClusterHCAL__RECO.obj.depth()>>dclust","","goff")
#events.Draw("recoPFClusters_particleFlowClusterHBHE__RECO.obj.energy()>>eclust2","","goff")
#events.Draw("recoPFClusters_particleFlowClusterHBHE__RECO.obj.time()>>tclust2","","goff")
#events.Draw("recoPFClusters_particleFlowClusterHBHE__RECO.obj.layer()>>lclust2","","goff")
#events.Draw("recoPFClusters_particleFlowClusterHBHE__RECO.obj.depth()>>dclust2","","goff")
#events.Draw("HBHERecHitsSorted_reducedHcalRecHits_hbhereco_RECO.obj.obj.energy()>>ehit","","goff")
#events.Draw("HBHERecHitsSorted_reducedHcalRecHits_hbhereco_RECO.obj.obj.time()>>thit","","goff")
#events.Draw("(((HBHERecHitsSorted_reducedHcalRecHits_hbhereco_RECO.obj.obj.detid().rawId())>>20)&0xF)>>dhit","","goff")
#events.Draw("(((HBHERecHitsSorted_reducedHcalRecHits_hbhereco_RECO.obj.obj.detid().rawId()))&0x3FF)>>phit","","goff")
#events.Draw("(((HBHERecHitsSorted_reducedHcalRecHits_hbhereco_RECO.obj.obj.detid().rawId())>>10)&0x1FF)>>nhit","","goff")
#events.Draw("recoPFRecHits_particleFlowRecHitHBHE__RECO.obj.energy()>>ehit2","","goff")
#events.Draw("recoPFRecHits_particleFlowRecHitHBHE__RECO.obj.time()>>thit2","","goff")
#events.Draw("recoPFRecHits_particleFlowRecHitHBHE__RECO.obj.layer()>>lhit2","","goff")
#events.Draw("recoPFRecHits_particleFlowRecHitHBHE__RECO.obj.depth()>>dhit2","","goff")
#events.Draw("PCaloHits_g4SimHits_HcalHits_SIM.obj.energy()*178.0>>ehit3","","goff")
events.Draw("PCaloHits_g4SimHits_HcalHits_SIM.obj.time()>>thit3","","goff")
#events.Draw("PCaloHits_g4SimHits_HcalHits_SIM.obj.time()>>thit4","","goff")
events.Draw("(((PCaloHits_g4SimHits_HcalHits_SIM.obj.id())>>21)&0x1F)>>dhit3","","goff")
events.Draw("(((PCaloHits_g4SimHits_HcalHits_SIM.obj.id()))&0x3FF)>>phit3","","goff")
events.Draw("(((PCaloHits_g4SimHits_HcalHits_SIM.obj.id())>>10)&0x1FF)>>nhit3","","goff")

#events = Events("SinglePiPt100RECO.root")
#events = Events("SinglePiPt50ieta23RECO.root")
#events = Events("SinglePiPt1-5ieta22-24RECO.root")
#events = Events("SinglePiPt5-15ieta22-24RECO.root")
#events = Events("output_new_tau0_step0_1000.root")
#events = Events("output_new_step0_1000.root")
#events = Events("output_new_step0_100.root")
#events = Events("output_decay_tau0_step0_1000.root")
#events = Events("output_tau0_step0_1000.root")
#events = Events("output_step0_1000.root")
#events = Events("output_step0_100.root")
#events = Events("output_step0.root")
#events = Events("step1tdc_prompt.root")
#events = Events("/eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/mh2000_mx975_pl10000_ev1000/ppTohToSS1SS2_SS1Tobb_SS2Toveve_1_withISR_step1_TDC.root")
#events = Events("/eos/cms/store/user/lowang/mh1000_pl1000_step1.root")
#events = Events("/eos/cms/store/user/lowang/LLP_htobbbb/step1/MH-125_MFF-50_CTau-1000mm_step1.root")
#events = Events("/eos/cms/store/user/lowang/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_step1.root")
events = Events("/eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/PionGun/SinglePion211_E10_PU_00_eta1phi0_timeslew-false_step1_CaloSamples_10events.root") #SinglePion211_E10_PU_00_eta1phi0_timeslew-false_step1.root")
#events = Events("/eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/PionGun/SinglePion211_E10_PU_00_step1.root")
#events = Events("/eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/PionGun/SinglePion211_E10_PU_step1.root")
#events = Events("step1tdc.root")

handleSim  = Handle ("std::vector<PCaloHit>")
labelSim = ("g4SimHits", "HcalHits", "SIM")

handleGen  = Handle ("std::vector<reco::GenParticle>")
#labelGen = ("genParticles", "", "SIM")
labelGen = ("genParticles", "", "HLT")

handlexyz0 = Handle("ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag>")
labelxyz0 = ("genParticles","xyz0")

handlet0 = Handle("float")
labelt0 = ("genParticles","t0")

#QIE11DataFrameHcalDataFrameContainer_simHcalUnsuppressedDigis_HBHEQIE11DigiCollection_HLT
handleQIE11Digi = Handle ("HcalDataFrameContainer<QIE11DataFrame>")
labelQIE11Digi = ("simHcalUnsuppressedDigis", "HBHEQIE11DigiCollection", "HLT")

#handleReco  = Handle ("edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >")
#handleReco  = Handle ("std::vector<HBHERecHitsSorted>")
#labelReco = ("hcalRecHits", "hbhereco", "RECO")
#handleReco  = Handle ("std::vector<HBHERecHitsSorted>")
#labelReco = ("reducedHcalRecHits", "hbhereco", "RECO")
#labelReco = ("hbhereco", "", "RECO")
#handleReco  = Handle ("std::vector<recoPFRecHits>")
#labelReco = ("particleFlowRecHitHBHE", "", "RECO")

#mdepthimg = np.zeros((6,5,36))
mdepthimg = np.zeros((1000,6,5,36))
nevents = 0
sevents = 0
radiusHB = 179.00
# depth layer radius and distances given here: https://cmssdt.cern.ch/lxr/source/SimG4CMS/Calo/src/HcalTestAnalysis.cc
rLay = np.zeros(19)
rLay[0] = 183.60

rLay[1] = 190.20
rLay[2] = 196.20
rLay[3] = 202.20
rLay[4] = 208.20

rLay[5] = 214.20
rLay[6] = 220.20
rLay[7] = 226.20
rLay[8] = 232.20
rLay[9] = 238.20

rLay[10] = 244.80
rLay[11] = 251.40
rLay[12] = 258.00
rLay[13] = 264.60
rLay[14] = 271.20
rLay[15] = 277.60
rLay[16] = 286.25

rLay[17] = 384.70 # HO
rLay[18] = 405.20 # HO
rDepth = np.zeros(4)
rDepth[0] = 183.60
rDepth[1] = 190.20
rDepth[2] = 214.20
rDepth[3] = 244.80

zDepth = np.zeros(7)
zDepth[0] = 403.4
zDepth[1] = 403.2
zDepth[2] = 429.7
zDepth[3] = 455.8
zDepth[4] = 481.9
zDepth[5] = 516.7
zDepth[6] = 516.7

tofadj = np.zeros(4)
#tofadj[0] = -3.5
#tofadj[1] = -3.5
#tofadj[2] = -4.5
#tofadj[3] = -5.5

tofadjHE = np.zeros(6)
#tofadjHE[0] = -10.5
#tofadjHE[1] = -10.5
#tofadjHE[2] = -11.5
#tofadjHE[3] = -12
#tofadjHE[4] = -13
#tofadjHE[5] = -14


for event in events:
    nevents+=1
#    if(nevents>1):
#       break
    print "\n PROCESSING EVENT ", nevents

    event.getByLabel(labelGen, handleGen)
    event.getByLabel (labelxyz0, handlexyz0)
    event.getByLabel (labelt0, handlet0)

    t0 = handlet0.product()
    v0 = handlexyz0.product()

    print " vx", v0.x(), " vy", v0.y(), " vz", v0.z(), " time",t0[0]

    genParticles = handleGen.product()

    llp=False
    maxpt=0.
    xPos=0.
    yPos=0.
    zPos=0.
    igen=0
    intereta=0.
    interphi=0.
    for gen in genParticles:
        igen+=1
        if (igen <= -1):
           print igen, gen.pdgId(), gen.numberOfMothers(), gen.numberOfDaughters(), gen.pt(), gen.px(), gen.py(), gen.pz(), gen.pt(), gen.p(), gen.energy(), gen.mass(), gen.eta(), gen.phi(), gen.vx(), gen.vy(), gen.vz()
           if (gen.numberOfMothers() > 0):
              print igen, " mothers ", gen.mother(0).pdgId()
           if (gen.numberOfDaughters() > 1):
              print igen, " daughters ", gen.daughter(0).pdgId(), gen.daughter(1).pdgId()
        decaylength = math.sqrt((gen.vx()-v0.x())**2+(gen.vy()-v0.y())**2+(gen.vz()-v0.z())**2)
        if (gen.numberOfDaughters() > 1):
            decaylength = math.sqrt((gen.daughter(0).vx()-v0.x())**2+(gen.daughter(0).vy()-v0.y())**2+(gen.daughter(0).vz()-v0.z())**2)
        if (gen.numberOfMothers() > 0):
            if (gen.mother(0).pdgId() == 35 or gen.mother(0).pdgId() == 36 or gen.pdgId() == 6000113):

                print " Particle ", gen.pdgId(), " from ", gen.mother(0).pdgId(), " decaylength = ", decaylength
                if (decaylength>1.):
                    print " Long-lived particle ***********"
                if (decaylength>1. and gen.pt() > 10.):
                    llp=True
#                    print(" gen vertex ", gen.vx(), gen.vertex().X())
                    momdecaydotprod = (gen.vx()-v0.x())*gen.px()+(gen.vy()-v0.y())*gen.py()+(gen.vz()-v0.z())*gen.pz()
                    cosmd = momdecaydotprod/(decaylength*gen.p())
                    if (gen.pt() > maxpt):
                        intereta, interphi = intersect.intersect(gen.vx()-v0.x(),gen.vy()-v0.y(),gen.vz()-v0.z(),gen.px(),gen.py(),gen.pz())
                        xPos = gen.vx()
                        yPos = gen.vy()
                        zPos = gen.vz()
                        maxpt = gen.pt()
                    vhit5.Fill(decaylength)
                    chit1.Fill(cosmd)
                    if (cosmd < 0.9):
                        chit3.Fill(cosmd)
                vhit4.Fill(decaylengths)
#                print gen.status(), gen.pdgId(), gen.pt(), gen.px(), gen.py(), gen.pz(), gen.energy(), gen.mass(), gen.eta(), gen.phi(), gen.vx(), gen.vy(), gen.vz()

    event.getByLabel(labelSim, handleSim)
    sim = handleSim.product()

# void HcalTestNumbering::unpackHcalIndex(
#   const uint32_t& idx, int& det, int& z, int& depth, int& eta, int& phi, int& lay) {
#   det = (idx >> 28) & 15;
#   depth = (idx >> 26) & 3;
#   depth += 1;
#  lay = (idx >> 21) & 31;
#  lay += 1;
#  z = (idx >> 20) & 1;
#  eta = (idx >> 10) & 1023;
#  phi = (idx & 1023);
# }
#   double theta = 2.0 * atan(exp(-eta));
#   double dist = 0.;
#   if (det == static_cast<int>(HcalBarrel)) {
#     const double rLay[19] = {1836.0,
#                              1902.0,
#                              1962.0,
#                              2022.0,
#                              2082.0,
#                              2142.0,
#                              2202.0,
#                              2262.0,
#                              2322.0,
#                              2382.0,
#                              2448.0,
#                              2514.0,
#                              2580.0,
#                              2646.0,
#                              2712.0,
#                              2776.0,
#                              2862.5,
#                              3847.0,
#                              4052.0};
#     if (layer > 0 && layer < 20)
#       dist += rLay[layer - 1] * mm / sin(theta);

    print sim
    print "NUMBER OF SIM HITS", len(sim)

    mult1ns3gev=0
    mult2ns3gev=0
    simenergy = 0.0
    ewdt = np.zeros((18,4,2))
    ew = np.zeros((18,4,2))
    for simhit in sim:
        simenergy+=simhit.energy()*117.0
        ehit3.Fill(simhit.energy()*117.0)
#HE        simenergy+=simhit.energy()*178.0
#HE        ehit3.Fill(simhit.energy()*178.0)
        thit4.Fill(simhit.time())
        lay = ((simhit.id())>>21)&0x1F
        iphi = (simhit.id())&0x3FF
        ieta = ((simhit.id())>>10)&0x1FF
        zside = (simhit.id()>>20)&0x1
        if (ieta < 15):
            phi = phiVal(iphi)
#            if (zside==0):
#              ieta = -ieta 
            eta = etaVal(ieta)
            theta = 2.*math.atan(math.exp(-eta))
            distance = rLay[lay] / math.sin(theta)
            tof = 0
#            tof = 1e9*distance/2.99793e10  # cm/sec
#            delayedtime = simhit.time() - t0[0] - tof
            delayedtime = simhit.time() - tof
            chit2d.Fill(delayedtime,1./math.sin(theta))
            ewdt[(iphi-1)/4,(ieta-1)/4,zside]+=delayedtime*simhit.energy()*117.0
#HE            ewdt[(iphi-1)/4,(ieta-1)/4,zside]+=delayedtime*simhit.energy()*178.0
            ew[(iphi-1)/4,(ieta-1)/4,zside]+=simhit.energy()*117.0
#HE            ew[(iphi-1)/4,(ieta-1)/4,zside]+=simhit.energy()*178.0
            if (simhit.energy()*117.0>3.0):
#HE            if (simhit.energy()*178.0>5.0):
                thit8.Fill(delayedtime)
                thit9.Fill(delayedtime)
                if (delayedtime>1.):
                    mult1ns3gev+=1
                if (delayedtime>2.):
                    mult2ns3gev+=1
            if (llp):
                if (intersect.deltaR(eta,phi,intereta,interphi)<0.5):
                    intertheta = 2.*math.atan(math.exp(-intereta))
                    interpath1 = math.sqrt(xPos**2+yPos**2+zPos**2)
                    interpath2 = math.sqrt((radiusHB*math.cos(interphi)-xPos)**2+(radiusHB*math.sin(interphi)-yPos)**2+(radiusHB/min((1e9,math.tan(intertheta)))-zPos)**2)
                    directpath = radiusHB / math.sin(theta)
                    pathdifftime = 1e9*(interpath1 + interpath2 - directpath)/2.99793e10 # cm/sec
                    qhit1.Fill(pathdifftime)
                    qhit2d.Fill(pathdifftime,delayedtime)
                    if (simhit.energy()*117.0>3.0):
#HE                    if (simhit.energy()*178.0>5.0):
                        qhit3.Fill(pathdifftime)
                        qhit4d.Fill(pathdifftime,delayedtime)
                    thit5.Fill(delayedtime)
#                    print " sim (layer,iphi,ieta,zside) ", lay, iphi, ieta, zside, " d, t, tof, ns ", distance, theta, tof, delayedtime
                else:
                    thit6.Fill(delayedtime)
#timemap[lay].Fill(iphi,ieta,delayedtime)
    mhit1.Fill(mult1ns3gev)
    mhit3.Fill(mult2ns3gev)
    for iphi in range(18):
        for ieta in range(4):
            for zside in range(2):
                phi = phiVal(4*iphi+2)
                eta = etaVal(4*ieta+2)
                theta = 2.*math.atan(math.exp(-eta))
                bigenergy = ew[iphi,ieta,zside]
                directpath = radiusHB / math.sin(theta)
                interpath1 = math.sqrt(xPos**2+yPos**2+zPos**2)
                ehit4.Fill(bigenergy)
                if (bigenergy>10.0):
                    bigdelay = ewdt[iphi,ieta,zside]/bigenergy
                    qhit5d.Fill(bigenergy,bigdelay)
                    if (bigenergy>40.0):
                        thit7.Fill(bigdelay)
                    if (llp):
                        if (intersect.deltaR(eta,phi,intereta,interphi)<0.5):
                            interpath2 = math.sqrt((radiusHB*math.cos(interphi)-xPos)**2+(radiusHB*math.sin(interphi)-yPos)**2+(radiusHB/min((1e9,math.tan(intertheta)))-zPos)**2)
                            pathdifftime = 1e9*(interpath1 + interpath2 - directpath)/2.99793e10 # cm/sec
                            qhit6d.Fill(pathdifftime,delayedtime)

    event.getByLabel(labelQIE11Digi,handleQIE11Digi)
    qie11digi = handleQIE11Digi.product()

    print qie11digi
    print "NUMBER OF QIE11 Digi", len(qie11digi)

#       QIE11DataFrame qiedf = static_cast<QIE11DataFrame>(qie11dc[j]);
#       DetId detid = qiedf.detid();

    mult1ns1000fC=0
    mult2ns1000fC=0
    for i in range(len(qie11digi)):
#        print qie11digi[i]
#        print qie11digi[i].id(), qie11digi[i].size()
        id = qie11digi[i].id()
        zside = (id >> 19) & 0x1
        depth = (id >> 20) & 0xF
        iphi = id & 0x3FF
        ieta = (id >> 10) & 0x1FF
#        print ' qie11 digi', zside, depth, iphi, ieta
        if (ieta < 15):
          phi = phiVal(iphi)
          eta = etaVal(ieta)
#          if (zside==0):
#            ieta = -ieta 
          theta = 2.*math.atan(math.exp(-eta))
          distance = rDepth[depth-1] / math.sin(theta)
          tof = 0
#          tof = 1e9*distance/2.99793e10  # cm/sec
          tof += tofadj[depth-1]
          qie11fCtot = 0.0
          qie11tdc = 0.0
          for j in range(qie11digi[i].size()):
#              adc = qie11digi[i][j].adc()
#              tdc = qie11digi[i][j].tdc()
#              capid = qie11digi[i][j].capid()
#              soi = qie11digi[i][j].soi()
              digi = qie11digi[i][j]
              adc = digi & 0xFF
#              tdc = (digi >> 8) & 0x3
#              capid = (digi >> 10) & 0x3
              tdc = (digi >> 8) & 0x3F
              capid = (digi >> 12) & 0x3
              soi = (digi >> 14) & 0x1
              if (soi == 1):
                  qie11tdc = tdc
              adchit1.Fill(adc)
              tdchit1.Fill(tdc)
              tdchit2d.Fill(tdc,tof)
              if (j > 3 and j < 8):
                  if ( adc < 16 ):
                      qie11fC = ((adc-4))*3.1
                  elif ( adc < 36 ):
                      qie11fC = ((15-4) + (adc-15)*2)*3.1
                  elif ( adc < 57 ):
                      qie11fC = ((15-4) + (35-15)*2 + (adc-35)*4)*3.1
                  else:
                      qie11fC = ((15-4) + (35-15)*2 + (56-35)*4 + (adc-56)*8)*3.1
                  qie11fCtot+=qie11fC
          fhit1.Fill(qie11fCtot)
          if (qie11fCtot > 1000.0):
              tdchit3.Fill(qie11tdc)
              delayedtime = qie11tdc*0.5-tof
              if (depth == 1):
                  tdchit4.Fill(delayedtime)
              elif (depth == 2):
                  tdchit5.Fill(delayedtime)
              elif (depth == 3):
                  tdchit6.Fill(delayedtime)
              else:
                  tdchit7.Fill(delayedtime)
              if (delayedtime>1.):
                  mult1ns1000fC+=1
              if (delayedtime>2.):
                  mult2ns1000fC+=1
          if (qie11fCtot > 3000.0):
              delayedtime = qie11tdc*0.5-tof
              if (depth == 1):
                  tdchit4_3kfC.Fill(delayedtime)
              elif (depth == 2):
                  tdchit5_3kfC.Fill(delayedtime)
              elif (depth == 3):
                  tdchit6_3kfC.Fill(delayedtime)
              else:
                  tdchit7_3kfC.Fill(delayedtime)

        if (ieta > 15): # HE, to see if need tdc_adjust to shift prompt peak to 0
          phi = phiVal(iphi)
          eta = etaVal(ieta)
          theta = 2.*math.atan(math.exp(-eta))
          distance = zDepth[depth-1] / math.cos(theta)
          tof = 0
#          tof = 1e9*distance/2.99793e10  # cm/sec
          if (depth <= 5): tof += tofadjHE[depth-1]
          else: tof += tofadjHE[5]
          qie11fCtot = 0.0
          qie11tdc = 0.0
          for j in range(qie11digi[i].size()):
              digi = qie11digi[i][j]
              adc = digi & 0xFF
              tdc = (digi >> 8) & 0x3F
              capid = (digi >> 12) & 0x3
              soi = (digi >> 14) & 0x1
              if (soi == 1):
                  qie11tdc = tdc
              if (j > 3 and j < 8):
                  if ( adc < 16 ):
                      qie11fC = ((adc-4))*3.1
                  elif ( adc < 36 ):
                      qie11fC = ((15-4) + (adc-15)*2)*3.1
                  elif ( adc < 57 ):
                      qie11fC = ((15-4) + (35-15)*2 + (adc-35)*4)*3.1
                  else:
                      qie11fC = ((15-4) + (35-15)*2 + (56-35)*4 + (adc-56)*8)*3.1
                  qie11fCtot+=qie11fC
          if (qie11fCtot > 1000.0):
              delayedtime = qie11tdc*0.5-tof
              if (depth == 1):
                  tdchit4HE.Fill(delayedtime)
              elif (depth == 2):
                  tdchit5HE.Fill(delayedtime)
              elif (depth == 3):
                  tdchit6HE.Fill(delayedtime)
              elif (depth == 4):
                  tdchit7HE.Fill(delayedtime)
              elif (depth == 5):
                  tdchit8HE.Fill(delayedtime)
              else:
                  tdchit9HE.Fill(delayedtime)
#              if (delayedtime>1.):
#                  mult1ns1000fC+=1
#              if (delayedtime>2.):
#                  mult2ns1000fC+=1
          if (qie11fCtot > 3000.0):
              delayedtime = qie11tdc*0.5-tof
              if (depth == 1):
                  tdchit4HE_3kfC.Fill(delayedtime)
              elif (depth == 2):
                  tdchit5HE_3kfC.Fill(delayedtime)
              elif (depth == 3):
                  tdchit6HE_3kfC.Fill(delayedtime)
              elif (depth == 4):
                  tdchit7HE_3kfC.Fill(delayedtime)
              elif (depth == 5):
                  tdchit8HE_3kfC.Fill(delayedtime)
              else:
                  tdchit9HE_3kfC.Fill(delayedtime)

    mhit5.Fill(mult1ns1000fC)
    mhit7.Fill(mult2ns1000fC)
            

# QIE11
# pedestal ~4.5 (between codes 4 and 5)
# 0-15 3.1 fC/bin
# 16-35 6.2 fC/bin
# 36-56 12.4 fC/bin
# 57-255 24.8 fC/bin
# if ( adc < 16 ):
#   qie11fC = (adc-4)*3.1
# elif ( adc < 36 ):
#   qie11fC = (15-4)*3.1 + (adc-15)*6.2
# elif ( adc < 57 ):
#   qie11fC = (15-4)*3.1 + (35-15)*6.2 + (adc-35)*12.4
# else:
#   qie11fC = (15-4)*3.1 + (35-15)*6.2 + (56-35)*12.4 + (adc-56)*24.8
#              print " digi ", j, adc, tdc, capid, soi

# MASK_ADC = 0xFF;
# MASK_TDC_HE = 0x3F;
# MASK_TDC_HB = 0x3;
# OFFSET_TDC = 8;  // 8 bits
# MASK_SOI = 0x4000;
# MASK_LE_HB = 0x2000;
# MASK_CAPID = 0x3;
# MASK_CAPID_INV_HB = 0xF3FF;
# MASK_CAPID_KEEP_HB = 0x0C00;
# OFFSET_CAPID_HE = 8;
# OFFSET_CAPID_HB = 10;

# kHcalPhiMask1 = 0x7F;
# kHcalPhiMask2 = 0x3FF;
# kHcalEtaOffset1 = 7;
# kHcalEtaOffset2 = 10;
# kHcalEtaMask1 = 0x3F;
# kHcalEtaMask2 = 0x1FF;
# kHcalZsideMask1 = 0x2000;
# kHcalZsideMask2 = 0x80000;
# kHcalDepthOffset1 = 14;
# kHcalDepthOffset2 = 20;
# kHcalDepthMask1 = 0x1F;
# kHcalDepthMask2 = 0xF;
# kHcalDepthSet1 = 0x1C000;
# kHcalDepthSet2 = 0xF00000;
# kHcalIdFormat2 = 0x1000000;
# kHcalIdMask = 0xFE000000;

# std::ostream& operator<<(std::ostream& s, const QIE11DataFrame& digi) {
#   if (digi.detid().det() == DetId::Hcal) {
#     s << "DetID=" << HcalGenericDetId(digi.detid()) << " flavor=" << digi.flavor();
#   } else {
#     s << "DetId(" << digi.detid().rawId() << ")";
#   }
#   s << " " << digi.samples() << " samples";
#   if (digi.linkError())
#     s << " LinkError ";
#   if (digi.capidError())
#     s << " CapIdError ";
#   if (digi.zsMarkAndPass())
#     s << " M&P ";
#   s << std::endl;
#   for (int i = 0; i < digi.samples(); i++) {
#     QIE11DataFrame::Sample sam = digi[i];
#     s << "  ADC=" << sam.adc() << " TDC=" << sam.tdc() << " CAPID=" << sam.capid();
#     if (sam.soi())
#       s << " SOI ";
#     s << std::endl;
#   }

#    event.getByLabel(labelReco, handleReco)
#    event.getByLabel("HcalRecHit", "hbhereco", "RECO", handleReco)
#    reco = handleReco.product()

#    print reco
#    print "NUMBER OF RECO HITS", len(reco)
#    if (len(reco)>0):
#        sevents+=1

    recoenergy = 0.0
#    for hit in reco:
#        recoenergy+=hit.energy()
#        ehit.Fill(hit.energy())
#        thit.Fill(hit.time())
#        dhit.Fill(double(((hit.detid().rawId())>>20)&0xF))
#        phit.Fill(double(((hit.detid().rawId()))&0x3FF))
#        nhit.Fill(double(((hit.detid().rawId())>>10)&0x1FF))
#        if (nevents==7): 
#        idepth = ((hit.detid().rawId())>>20)&0xF
#        iphi = (((hit.detid().rawId()))&0x3FF-1)/2+1
#        ieta = ((hit.detid().rawId())>>10)&0x1FF
#        if (hit.energy() > 0.0):
#            print ' idepth, iphi, ieta, energy',idepth,iphi,ieta,hit.energy()
#        if ((iphi > 0) and (iphi <= 36) and (ieta > 20) and (ieta <= 25) and (idepth > 0) and (idepth <= 6)):
#            mdepthimg[sevents-1,idepth-1,ieta-21,iphi-1] = hit.energy()
#            if (idepth == 1):
#                mdepthimg[sevents-1,ieta-21,iphi-1] = hit.energy()

    print "sim/reco energy=",simenergy,recoenergy

#    if (nevents==7): 
#        print(mdepthimg.shape)
#        tifffile.imwrite('mdepth.tif',mdepthimg)
#        mdepthimg2 = tifffile.imread('mdepth.tif').astype(np.float32)
#        print(mdepthimg2.shape)
#        print(mdepthimg2[0,0,0],mdepthimg2[0,2,18])
#        mdepthimg3 = tifffile.imread('mdepth.tif')
#        print(mdepthimg3.shape)
#        mdepthimg3 = (mdepthimg3 > 0).astype(np.float32)
#        print(mdepthimg3[0,0,0],mdepthimg3[0,2,18])
#        break
#tifffile.imwrite('mdepth.tif',mdepthimg)

gStyle.SetOptStat(1111100);

c3 = TCanvas()
vt.Draw("EHIST")
vt.GetXaxis().SetTitle("t")
vt.GetYaxis().SetTitle("# of v")
c3.SaveAs("vt.png")

c4 = TCanvas()
vx.Draw("EHIST")
vx.GetXaxis().SetTitle("x")
vx.GetYaxis().SetTitle("# of v")
c4.SaveAs("vx.png")

c5 = TCanvas()
vy.Draw("EHIST")
vy.GetXaxis().SetTitle("y")
vy.GetYaxis().SetTitle("# of v")
c5.SaveAs("vy.png")

c6 = TCanvas()
vz.Draw("EHIST")
vz.GetXaxis().SetTitle("z")
vz.GetYaxis().SetTitle("# of v")
c6.SaveAs("vz.png")

#c3 = TCanvas()
#c3.SetLogy(True)
#eclust2.Draw("EHIST")
#eclust2.GetXaxis().SetTitle("energy")
#eclust2.GetYaxis().SetTitle("# of entries")
#c3.SaveAs("eclust2.png")

#c4 = TCanvas()
#tclust2.Draw("EHIST")
#tclust2.GetXaxis().SetTitle("time")
#tclust2.GetYaxis().SetTitle("# of entries")
#c4.SaveAs("tclust2.png")

#c5 = TCanvas()
#lclust2.Draw("EHIST")
#lclust2.GetXaxis().SetTitle("layer")
#lclust2.GetYaxis().SetTitle("# of entries")
#c5.SaveAs("lclust2.png")

#c6 = TCanvas()
#dclust2.Draw("EHIST")
#dclust2.GetXaxis().SetTitle("depth")
#dclust2.GetYaxis().SetTitle("# of entries")
#c6.SaveAs("dclust2.png")

#c7 = TCanvas()
#c7.SetLogy(True)
#eclust.Draw("EHIST")
#eclust.GetXaxis().SetTitle("energy")
#eclust.GetYaxis().SetTitle("# of entries")
#c7.SaveAs("eclust.png")

#c8 = TCanvas()
#tclust.Draw("EHIST")
#tclust.GetXaxis().SetTitle("time")
#tclust.GetYaxis().SetTitle("# of entries")
#c8.SaveAs("tclust.png")

#c9 = TCanvas()
#lclust.Draw("EHIST")
#lclust.GetXaxis().SetTitle("layer")
#lclust.GetYaxis().SetTitle("# of entries")
#c9.SaveAs("lclust.png")

#c10 = TCanvas()
#dclust.Draw("EHIST")
#dclust.GetXaxis().SetTitle("depth")
#dclust.GetYaxis().SetTitle("# of entries")
#c10.SaveAs("dclust.png")

#c11 = TCanvas()
#ehit.Draw("EHIST")
#ehit.GetXaxis().SetTitle("energy")
#ehit.GetYaxis().SetTitle("# of entries")
#c11.SaveAs("ehit.png")

#c12 = TCanvas()
#thit.Draw("EHIST")
#thit.GetXaxis().SetTitle("time")
#thit.GetYaxis().SetTitle("# of entries")
#c12.SaveAs("thit.png")

#c13 = TCanvas()
#c13.SetLogy(True)
#thit.Draw("EHIST")
#thit.GetXaxis().SetTitle("time")
#thit.GetYaxis().SetTitle("# of entries")
#c13.SaveAs("thitlogy.png")

#c14 = TCanvas()
#dhit.Draw("EHIST")
#dhit.GetXaxis().SetTitle("depth")
#dhit.GetYaxis().SetTitle("# of entries")
#c14.SaveAs("dhit.png")

#c1 = TCanvas()
#phit.Draw("EHIST")
#phit.GetXaxis().SetTitle("iphi")
#phit.GetYaxis().SetTitle("# of entries")
#c1.SaveAs("phit.png")

#c2 = TCanvas()
#nhit.Draw("EHIST")
#nhit.GetXaxis().SetTitle("ieta")
#nhit.GetYaxis().SetTitle("# of entries")
#c2.SaveAs("nhit.png")

#c15 = TCanvas()
#ehit2.Draw("EHIST")
#ehit2.GetXaxis().SetTitle("energy")
#ehit2.GetYaxis().SetTitle("# of entries")
#c15.SaveAs("ehit2.png")

#c16 = TCanvas()
#thit2.Draw("EHIST")
#thit2.GetXaxis().SetTitle("time")
#thit2.GetYaxis().SetTitle("# of entries")
#c16.SaveAs("thit2.png")

#c17 = TCanvas()
#lhit2.Draw("EHIST")
#lhit2.GetXaxis().SetTitle("layer")
#lhit2.GetYaxis().SetTitle("# of entries")
#c17.SaveAs("lhit2.png")

#c18 = TCanvas()
#dhit2.Draw("EHIST")
#dhit2.GetXaxis().SetTitle("depth")
#dhit2.GetYaxis().SetTitle("# of entries")
#c18.SaveAs("dhit2.png")

c19 = TCanvas()
ehit3.Draw("EHIST")
ehit3.GetXaxis().SetTitle("energy")
ehit3.GetYaxis().SetTitle("# of entries")
c19.SaveAs("ehit3.png")

c24 = TCanvas()
c24.SetLogy(True)
ehit3.Draw("EHIST")
ehit3.GetXaxis().SetTitle("energy")
ehit3.GetYaxis().SetTitle("# of entries")
c24.SaveAs("ehit3logy.png")

c20 = TCanvas()
thit3.Draw("EHIST")
thit3.GetXaxis().SetTitle("time")
thit3.GetYaxis().SetTitle("# of entries")
c20.SaveAs("thit3.png")

c21 = TCanvas()
dhit3.Draw("EHIST")
dhit3.GetXaxis().SetTitle("depth")
dhit3.GetYaxis().SetTitle("# of entries")
c21.SaveAs("dhit3.png")

c22 = TCanvas()
phit3.Draw("EHIST")
phit3.GetXaxis().SetTitle("phi")
phit3.GetYaxis().SetTitle("# of entries")
c22.SaveAs("phit3.png")

c23 = TCanvas()
nhit3.Draw("EHIST")
nhit3.GetXaxis().SetTitle("eta")
nhit3.GetYaxis().SetTitle("# of entries")
c23.SaveAs("nhit3.png")

c24 = TCanvas()
thit4.Draw("EHIST")
thit4.GetXaxis().SetTitle("time")
thit4.GetYaxis().SetTitle("# of entries")
c24.SaveAs("thit4.png")

c25 = TCanvas()
vhit4.Draw("EHIST")
vhit4.GetXaxis().SetTitle("distance")
vhit4.GetYaxis().SetTitle("# of entries")
c25.SaveAs("vhit4.png")

c26 = TCanvas()
vhit5.Draw("EHIST")
vhit5.GetXaxis().SetTitle("distance")
vhit5.GetYaxis().SetTitle("# of entries")
c26.SaveAs("vhit5.png")

c27 = TCanvas()
chit1.Draw("EHIST")
chit1.GetXaxis().SetTitle("cos(theta)")
chit1.GetYaxis().SetTitle("# of entries")
c27.SaveAs("chit1.png")

c28 = TCanvas()
c28.SetLogy(True)
thit5.Draw("EHIST")
thit5.GetXaxis().SetTitle("delayed time")
thit5.GetYaxis().SetTitle("# of entries")
c28.SaveAs("thit5.png")

c29 = TCanvas()
c29.SetLogy(True)
thit6.Draw("EHIST")
thit6.GetXaxis().SetTitle("delayed time")
thit6.GetYaxis().SetTitle("# of entries")
c29.SaveAs("thit6.png")

c30 = TCanvas()
#c30.SetLogy(True)
thit7.Draw("EHIST")
thit7.GetXaxis().SetTitle("delayed time")
thit7.GetYaxis().SetTitle("# of entries")
c30.SaveAs("thit7.png")

c31 = TCanvas()
thit8.Draw("EHIST")
thit8.GetXaxis().SetTitle("time")
thit8.GetYaxis().SetTitle("# of entries")
c31.SaveAs("thit8.png")

c32 = TCanvas()
chit2d.Draw("BOX")
chit2d.GetXaxis().SetTitle("1/sin(theta)")
chit2d.GetYaxis().SetTitle("delayed time")
c32.SaveAs("chit2d.png")

c33 = TCanvas()
c33.SetLogy(True)
thit9.Draw("EHIST")
thit9.GetXaxis().SetTitle("time")
thit9.GetYaxis().SetTitle("# of entries")
c33.SaveAs("thit9.png")

c34 = TCanvas()
qhit1.Draw("EHIST")
qhit1.GetXaxis().SetTitle("delayed time")
qhit1.GetYaxis().SetTitle("# of entries")
c34.SaveAs("qhit1.png")

c35 = TCanvas()
qhit2d.Draw("BOX")
qhit2d.GetXaxis().SetTitle("path length difference time delay")
qhit2d.GetYaxis().SetTitle("delayed time")
c35.SaveAs("qhit2d.png")

c36 = TCanvas()
chit3.Draw("EHIST")
chit3.GetXaxis().SetTitle("cos(theta)")
chit3.GetYaxis().SetTitle("# of entries")
c36.SaveAs("chit3.png")

c37 = TCanvas()
qhit3.Draw("EHIST")
qhit3.GetXaxis().SetTitle("delayed time")
qhit3.GetYaxis().SetTitle("# of entries")
c37.SaveAs("qhit3.png")

c38 = TCanvas()
qhit4d.Draw("BOX")
qhit4d.GetXaxis().SetTitle("path length difference time delay")
qhit4d.GetYaxis().SetTitle("delayed time")
c38.SaveAs("qhit4d.png")

c39 = TCanvas()
c39.SetLogy(True)
ehit4.Draw("EHIST")
ehit4.GetXaxis().SetTitle("energy")
ehit4.GetYaxis().SetTitle("# of entries")
c39.SaveAs("ehit4.png")

c40 = TCanvas()
qhit5d.Draw("BOX")
qhit5d.GetXaxis().SetTitle("CaloRegion energy")
qhit5d.GetYaxis().SetTitle("delayed time")
c40.SaveAs("qhit5d.png")

c41 = TCanvas()
qhit6d.Draw("BOX")
qhit6d.GetXaxis().SetTitle("path length difference time delay")
qhit6d.GetYaxis().SetTitle("CaloRegion delayed time")
c41.SaveAs("qhit6d.png")

c42 = TCanvas()
mhit1.Draw("EHIST")
mhit1.GetXaxis().SetTitle("multiplicity (1ns,3GeV)")
mhit1.GetYaxis().SetTitle("# of entries")
c42.SaveAs("mhit1.png")

c43 = TCanvas()
c43.SetLogy(True)
mhit1.Draw("EHIST")
mhit1.GetXaxis().SetTitle("multiplicity (1ns,3GeV)")
mhit1.GetYaxis().SetTitle("# of entries")
c43.SaveAs("mhit2.png")

c44 = TCanvas()
mhit3.Draw("EHIST")
mhit3.GetXaxis().SetTitle("multiplicity (2ns,3GeV)")
mhit3.GetYaxis().SetTitle("# of entries")
c44.SaveAs("mhit3.png")

c45 = TCanvas()
c45.SetLogy(True)
mhit3.Draw("EHIST")
mhit3.GetXaxis().SetTitle("multiplicity (2ns,3GeV)")
mhit3.GetYaxis().SetTitle("# of entries")
c45.SaveAs("mhit4.png")

c46 = TCanvas()
c46.SetLogy(True)
adchit1.Draw("EHIST")
adchit1.GetXaxis().SetTitle("ADC")
adchit1.GetYaxis().SetTitle("# of entries")
c46.SaveAs("adchit1.png")

c47 = TCanvas()
c47.SetLogy(True)
tdchit1.Draw("EHIST")
tdchit1.GetXaxis().SetTitle("TDC")
tdchit1.GetYaxis().SetTitle("# of entries")
c47.SaveAs("tdchit1.png")

c48 = TCanvas()
tdchit2d.Draw("BOX")
tdchit2d.GetXaxis().SetTitle("TDC versus TOF")
tdchit2d.GetYaxis().SetTitle("# of entries")
c48.SaveAs("tdchit2d.png")

c49 = TCanvas()
c49.SetLogy(True)
fhit1.Draw("EHIST")
fhit1.GetXaxis().SetTitle("fC")
fhit1.GetYaxis().SetTitle("# of entries")
c49.SaveAs("fhit1.png")

c50 = TCanvas()
tdchit3.Draw("EHIST")
tdchit3.GetXaxis().SetTitle("TDC")
tdchit3.GetYaxis().SetTitle("# of entries")
c50.SaveAs("tdchit3.png")

c51 = TCanvas()
tdchit4.Draw("EHIST")
tdchit4.GetXaxis().SetTitle("TDC(ns)")
tdchit4.GetYaxis().SetTitle("# of entries")
c51.SaveAs("tdchit4_1kfC.png")

c52 = TCanvas()
tdchit5.Draw("EHIST")
tdchit5.GetXaxis().SetTitle("TDC(ns)")
tdchit5.GetYaxis().SetTitle("# of entries")
c52.SaveAs("tdchit5_1kfC.png")

c53 = TCanvas()
tdchit6.Draw("EHIST")
tdchit6.GetXaxis().SetTitle("TDC(ns)")
tdchit6.GetYaxis().SetTitle("# of entries")
c53.SaveAs("tdchit6_1kfC.png")

c54 = TCanvas()
tdchit7.Draw("EHIST")
tdchit7.GetXaxis().SetTitle("TDC(ns)")
tdchit7.GetYaxis().SetTitle("# of entries")
c54.SaveAs("tdchit7_1kfC.png")

c51HE = TCanvas()
tdchit4HE.Draw("EHIST")
tdchit4HE.GetXaxis().SetTitle("TDC(ns)")
tdchit4HE.GetYaxis().SetTitle("# of entries")
c51HE.SaveAs("tdchit4HE_1kfC.png")

c52HE = TCanvas()
tdchit5HE.Draw("EHIST")
tdchit5HE.GetXaxis().SetTitle("TDC(ns)")
tdchit5HE.GetYaxis().SetTitle("# of entries")
c52HE.SaveAs("tdchit5HE_1kfC.png")

c53HE = TCanvas()
tdchit6HE.Draw("EHIST")
tdchit6HE.GetXaxis().SetTitle("TDC(ns)")
tdchit6HE.GetYaxis().SetTitle("# of entries")
c53HE.SaveAs("tdchit6HE_1kfC.png")

c54HE = TCanvas()
tdchit7HE.Draw("EHIST")
tdchit7HE.GetXaxis().SetTitle("TDC(ns)")
tdchit7HE.GetYaxis().SetTitle("# of entries")
c54HE.SaveAs("tdchit7HE_1kfC.png")

c55HE = TCanvas()
tdchit8HE.Draw("EHIST")
tdchit8HE.GetXaxis().SetTitle("TDC(ns)")
tdchit8HE.GetYaxis().SetTitle("# of entries")
c55HE.SaveAs("tdchit8HE_1kfC.png")

c56HE = TCanvas()
tdchit9HE.Draw("EHIST")
tdchit9HE.GetXaxis().SetTitle("TDC(ns)")
tdchit9HE.GetYaxis().SetTitle("# of entries")
c56HE.SaveAs("tdchit9HE_1kfC.png")

c51_3kfC = TCanvas()
tdchit4_3kfC.Draw("EHIST")
tdchit4_3kfC.GetXaxis().SetTitle("TDC(ns)")
tdchit4_3kfC.GetYaxis().SetTitle("# of entries")
c51_3kfC.SaveAs("tdchit4_3kfC.png")

c52_3kfC = TCanvas()
tdchit5_3kfC.Draw("EHIST")
tdchit5_3kfC.GetXaxis().SetTitle("TDC(ns)")
tdchit5_3kfC.GetYaxis().SetTitle("# of entries")
c52_3kfC.SaveAs("tdchit5_3kfC.png")

c53_3kfC = TCanvas()
tdchit6_3kfC.Draw("EHIST")
tdchit6_3kfC.GetXaxis().SetTitle("TDC(ns)")
tdchit6_3kfC.GetYaxis().SetTitle("# of entries")
c53_3kfC.SaveAs("tdchit6_3kfC.png")

c54_3kfC = TCanvas()
tdchit7_3kfC.Draw("EHIST")
tdchit7_3kfC.GetXaxis().SetTitle("TDC(ns)")
tdchit7_3kfC.GetYaxis().SetTitle("# of entries")
c54_3kfC.SaveAs("tdchit7_3kfC.png")

c51HE_3kfC = TCanvas()
tdchit4HE_3kfC.Draw("EHIST")
tdchit4HE_3kfC.GetXaxis().SetTitle("TDC(ns)")
tdchit4HE_3kfC.GetYaxis().SetTitle("# of entries")
c51HE_3kfC.SaveAs("tdchit4HE_3kfC.png")

c52HE_3kfC = TCanvas()
tdchit5HE_3kfC.Draw("EHIST")
tdchit5HE_3kfC.GetXaxis().SetTitle("TDC(ns)")
tdchit5HE_3kfC.GetYaxis().SetTitle("# of entries")
c52HE_3kfC.SaveAs("tdchit5HE_3kfC.png")

c53HE_3kfC = TCanvas()
tdchit6HE_3kfC.Draw("EHIST")
tdchit6HE_3kfC.GetXaxis().SetTitle("TDC(ns)")
tdchit6HE_3kfC.GetYaxis().SetTitle("# of entries")
c53HE_3kfC.SaveAs("tdchit6HE_3kfC.png")

c54HE_3kfC = TCanvas()
tdchit7HE_3kfC.Draw("EHIST")
tdchit7HE_3kfC.GetXaxis().SetTitle("TDC(ns)")
tdchit7HE_3kfC.GetYaxis().SetTitle("# of entries")
c54HE_3kfC.SaveAs("tdchit7HE_3kfC.png")

c55HE_3kfC = TCanvas()
tdchit8HE_3kfC.Draw("EHIST")
tdchit8HE_3kfC.GetXaxis().SetTitle("TDC(ns)")
tdchit8HE_3kfC.GetYaxis().SetTitle("# of entries")
c55HE_3kfC.SaveAs("tdchit8HE_3kfC.png")

c56HE_3kfC = TCanvas()
tdchit9HE_3kfC.Draw("EHIST")
tdchit9HE_3kfC.GetXaxis().SetTitle("TDC(ns)")
tdchit9HE_3kfC.GetYaxis().SetTitle("# of entries")
c56HE_3kfC.SaveAs("tdchit9HE_3kfC.png")







c55 = TCanvas()
mhit5.Draw("EHIST")
mhit5.GetXaxis().SetTitle("multiplicity (1ns,1000fC)")
mhit5.GetYaxis().SetTitle("# of entries")
c55.SaveAs("mhit5.png")

c56 = TCanvas()
c56.SetLogy(True)
mhit5.Draw("EHIST")
mhit5.GetXaxis().SetTitle("multiplicity (1ns,1000fC)")
mhit5.GetYaxis().SetTitle("# of entries")
c56.SaveAs("mhit6.png")

c57 = TCanvas()
mhit7.Draw("EHIST")
mhit7.GetXaxis().SetTitle("multiplicity (2ns,1000fC)")
mhit7.GetYaxis().SetTitle("# of entries")
c57.SaveAs("mhit7.png")

c58 = TCanvas()
c58.SetLogy(True)
mhit7.Draw("EHIST")
mhit7.GetXaxis().SetTitle("multiplicity (2ns,1000fC)")
mhit7.GetYaxis().SetTitle("# of entries")
c58.SaveAs("mhit8.png")

#raw_input("press enter to continue")
