import numpy as np
import ROOT

#generate photons just off a cone to see what the distribution looks like.
conetheta=.78
#width = .6

ROOT.gROOT.Reset()
c1=ROOT.TCanvas("c1")
widths=[]
conethetas=[]
normparam=[]
sqrtparam=[]
expabsparam=[]
expdiffparam=[]
gausparam=[]
coneparam=[]
normparame=[]
sqrtparame=[]
expabsparame=[]
expdiffparame=[]
gausparame=[]
coneparame=[]

for width in np.arange(.6,1,.4):
    coshist=ROOT.TH1F("coshist","coshist",1000,-1,1) 
    for i in range(100000):
        startphi=np.random.rand()*2*np.pi
        offtheta=np.random.normal(0,width)
        offphi=np.random.rand()*2*np.pi
        #start with a vector along the z axis and rotate by off-angles
        x = np.sin(offtheta)*np.cos(offphi)
        y = np.sin(offtheta)*np.sin(offphi)
        z = np.cos(offtheta)
    
        #now rotate it to the cone angle
        midx = x * np.cos(conetheta) - z * np.sin(conetheta)
        midy = y
        midz = x * np.sin(conetheta) + z * np.cos(conetheta)
    
        #now rotate it to the right startphi
        endx = midx * np.cos(startphi) + midy * np.sin(startphi)
        endy =-midx * np.sin(startphi) + midy * np.cos(startphi)
        endz = midz
    
        #and now histogram it
        coshist.Fill(endz)
    fit=ROOT.TF1("fit","[0]*exp([1]*pow(abs(x-[4]),[6])+[2]*abs(x-[4])+[3]*(x-[4])*(x-[4])+[5]*x)",-1,1)
    fit.SetParameters(700,-4,-4,0.5/width*width,np.cos(conetheta),0,.5)
    fit.FixParameter(0, (13.24*width*width - 27.52*width +18.49)*(370*conetheta+171)/441)
    fit.FixParameter(1, 3.07*width*width - 3.52*width -2.21)
    fit.FixParameter(2, -10.8*width*width + 17.2*width -5.5)
    fit.FixParameter(3, 2.83*width -3.21)
    fit.FixParameter(4, np.cos(conetheta))
    fit.FixParameter(5, (-2.6*width*width + 4.2*width -1)*(-1.05*conetheta +1.46)/.694)
    fit.FixParameter(6, .37)
    fit.SetNpx(10000)
    coshist.Draw()
    coshist.Scale(.01)
    coshist.Fit(fit)
#    c1.SetLogy()
    c1.Update()
    c1.SaveAs("coneAngle"+str(conetheta)+"blur"+str(width)+".png")
    widths.append(width)
    conethetas.append(conetheta)
    normparam.append(fit.GetParameter(0))
    sqrtparam.append(fit.GetParameter(1))
    expabsparam.append(fit.GetParameter(2))
    expdiffparam.append(fit.GetParameter(5))
    gausparam.append(fit.GetParameter(3))
    coneparam.append(fit.GetParameter(4))
    normparame.append(fit.GetParError(0))
    sqrtparame.append(fit.GetParError(1))
    expabsparame.append(fit.GetParError(2))
    expdiffparame.append(fit.GetParError(5))
    gausparame.append(fit.GetParError(3))
    coneparame.append(fit.GetParError(4))

print "widths:", widths
print "conethetas:", conethetas
print "norms:", normparam
print "sqrts:", sqrtparam
print "expabs:", expabsparam
print "expdiff:", expdiffparam
print "gaus:", gausparam
print "cone:", coneparam
print "normse:", normparame
print "sqrtse:", sqrtparame
print "expabse:", expabsparame
print "expdiffe:", expdiffparame
print "gause:", gausparame
print "conee:", coneparame
    
raw_input()
