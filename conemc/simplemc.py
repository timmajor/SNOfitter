import numpy as np
import ROOT

#generate photons just off a cone to see what the distribution looks like.
conetheta=np.pi/2
width = .01

ROOT.gROOT.Reset()
c1=ROOT.TCanvas("c1")

coshist=ROOT.TH1F("coshist","coshist",10000,-1,1) 
for i in range(100000):
    offx=np.random.normal(0,width)
    offy=np.random.normal(0,width)
    offphi=np.random.rand()*2*np.pi
    #start with a vector along the z axis and rotate by off-angles
    y = offr * np.cos(offphi) 

    coshist.Fill(offr)
coshist.Draw()
c1.SetLogy()
c1.Update()

raw_input()
