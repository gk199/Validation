import ROOT
from array import array
import numpy as np
import matplotlib.pyplot as plt
from root_numpy import hist2array, root2array, tree2array
from root_numpy import testdata

nvariables = 4
ndata = 2000

# read in signal and background files
background_filename = ROOT.TFile("rates_new_cond_QCD_4L1Jets.root")
background_tree = (background_filename.Get("MultForTMVA")) # Get() returns a ROOT.TTree, which is correct
#background_hist = (background_filename.Get("dt3GeV3nsHBJetMult_emu")) # Get() returns a ROOT.TH1F (correct)
background_data = tree2array(background_tree) # access histogram from the bkg ROOT file
signal_filename = ROOT.TFile("rates_new_cond_pl1000_4L1Jets.root") 
signal_tree = (signal_filename.Get("MultForTMVA"))
#signal_hist = (signal_filename.Get("dt3GeV3nsHBJetMult_emu")) # Get() returns a ROOT.TH1F (correct)   
signal_data = tree2array(signal_tree) # access histogram from the signal ROOT file  

MethodNames = ["Fisher", "SVM", "BDT"] #, "KNN"]
methodOptions = {"BDT":"VarTransform=G", "Fisher":"VarTransform=G", "SVM":"C=1.0:Gamma=0.005:Tol=0.001:VarTransform=None"} #, "KNN":"VarTransform=G"}
outputfile = "output.root"

#ROOT Stuff
datatree = ROOT.TTree("datatree", "datatree")
xa = []
for i in range(nvariables):
    xa.append(array('f', [0.0]))
    datatree.Branch('x'+str(i), xa[i], 'x'+str(i) + '/F')
trutha = array('b', [0])
datatree.Branch('truth', trutha, 'truth/B')

for x in signal_data:
    for j in range(nvariables):
        xa[j][0] = x[j]
    trutha[0] = 1
    datatree.Fill()
    
for x in background_data:
    for j in range(nvariables):
        xa[j][0] = x[j]
    trutha[0] = 0
    datatree.Fill()
    
f_out = ROOT.TFile(outputfile, "RECREATE")
ROOT.TMVA.Tools.Instance()
factory = ROOT.TMVA.Factory("TMVAClass", f_out, "AnalysisType=Classification")

dataloader = ROOT.TMVA.DataLoader("dataset")
for i in range(nvariables):
    dataloader.AddVariable('x'+str(i), "F")
dataloader.AddSignalTree(datatree)
dataloader.AddBackgroundTree(datatree)
SigCut = ROOT.TCut("truth>0.5")
BkgCut = ROOT.TCut("truth<=0.5")
dataloader.PrepareTrainingAndTestTree(SigCut, BkgCut, "SplitMode=Random:NormMode=NumEvents:!V")

for method in MethodNames:
    Type = "ROOT.TMVA.Types.k" + method
    factory.BookMethod(dataloader, eval(Type), method, "")
    factory.BookMethod(dataloader, eval(Type), method+"Gauss", methodOptions[method])
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
f_out.Close()
