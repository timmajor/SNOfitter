import numpy as np
import ROOT

STEP = 1 #mm

#Generate a photon in a random direction
def getPhoton():
    costhetai=1-2*np.random.rand()
    zstep = STEP * costhetai
    rhostep = STEP * np.sqrt(1-costhetai*costhetai)
    return zstep, rhostep

#Find where it meets an interface
def whereami(z, rho):
    r=np.sqrt(z*z+rho*rho)
    if r<6000:
        region=1
    elif r<6050:
        region=2
    elif r<8400:
        region=3
    else:
        region=4
    return region

def whatindex(region):
    indices=[0,1.33,1.5,1.33,1.33]
    return indices[region]

def walk(z, rho, zstep, rhostep):
    newz=z+zstep
    newrho=rho+rhostep
    return newz, newrho

def nextinterface(z, rho, zstep, rhostep):
    newz, newrho=walk(z, rho, zstep, rhostep)
    region=whereami(newz,newrho)
    while (region is whereami(newz, newrho)):
        newz, newrho = walk(newz, newrho, zstep, rhostep)
    return newz, newrho

#Where does it go?
def newdir(z, rho, zstep, rhostep, n1, n2, nrefl):
    #components
#    #print(np.sqrt(zstep*zstep+rhostep*rhostep))
    r=np.sqrt(z*z+rho*rho)
    rhat_z=z/r
    rhat_rho=rho/r
    radial1=rhat_z*zstep+rhat_rho*rhostep
    radial1_z=rhat_z*radial1
    radial1_rho=rhat_rho*radial1
    tangential1_z=zstep-radial1_z
    tangential1_rho=rhostep-radial1_rho
    #refract
    tangential2_z=n1/n2*tangential1_z
    tangential2_rho=n1/n2*tangential1_rho
    if (STEP*STEP>tangential2_z*tangential2_z+tangential2_rho*tangential2_rho):
#        print "refract"
        radial2=np.sqrt(STEP*STEP-tangential2_z*tangential2_z-tangential2_rho*tangential2_rho)
        radial2_z=radial2*rhat_z*np.sign(radial1)
        radial2_rho=radial2*rhat_rho*np.sign(radial1)
    #reflect
    else:
#        print "reflect"
        tangential2_z=tangential1_z
        tangential2_rho=tangential1_rho
        radial2_z=-radial1_z
        radial2_rho=-radial1_rho
        nrefl+=1
    newzstep=tangential2_z+radial2_z
    newrhostep=tangential2_rho+radial2_rho
#    print(np.sqrt(newzstep*newzstep+newrhostep*newrhostep))
    return newzstep, newrhostep, nrefl

def getfinalphoton():
    startzstep, startrhostep = getPhoton()
#    startzstep=0
#    startrhostep=1
#    print startzstep, startrhostep
    nrefl=0
    startz=6001+np.random.rand()*48
    startrho=0
    startregion=whereami(startz, startrho)
    thisz, thisrho = nextinterface(startz, startrho, startzstep, startrhostep)
#    print thisz/6000., thisrho/6000
    thisregion=whereami(thisz, thisrho)
    thiszstep=startzstep
    thisrhostep=startrhostep
    while thisregion is not 4:
        thiszstep, thisrhostep, nrefl = newdir(thisz, thisrho, thiszstep, thisrhostep, whatindex(whereami(thisz-thiszstep,thisrho-thisrhostep)), whatindex(whereami(thisz, thisrho)),nrefl)
#        print thiszstep, thisrhostep
        thisz, thisrho = nextinterface(thisz, thisrho, thiszstep, thisrhostep)
        startregion=thisregion
        thisregion=whereami(thisz, thisrho)
        if nrefl>0:
            return 0, 0
    costhetai = startzstep/STEP
    costhetaf = (thisz-startz)/np.sqrt((thisz-startz)*(thisz-startz)+thisrho*thisrho)
#    print costhetai, costhetaf
    return costhetai, costhetaf

c1=ROOT.TCanvas("c1","c1",1200,400)
c1.Divide(3,1)
histo=ROOT.TH2D("histo","histo",100,-1,1,100,-1,1)

#getfinalphoton()

#thiszstep=0
#thisrhostep=1
#thisz, thisrho = nextinterface(6025, 0, thiszstep, thisrhostep)
#startregion=whereami(6025,0)
#thisregion=whereami(thisz, thisrho)
#nrefl=0
#thiszstep, thisrhostep, nrefl = newdir(thisz, thisrho, thiszstep, thisrhostep, whatindex(startregion), whatindex(thisregion),nrefl)
#print thisz, thisrho
#print np.arccos(thiszstep)
#thisz, thisrho = nextinterface(thisz, thisrho, thiszstep, thisrhostep)
#thiszstep, thisrhostep nrefl = newdir(thisz, thisrho, thiszstep, thisrhostep, whatindex(whereami(thisz-thiszstep,thisrho-thisrhostep), whatindex(thisz,thisrho),nrefl)

for i in range(3000):
    print i
    costhetai, costhetaf = getfinalphoton()
    histo.Fill(costhetai, costhetaf)
histin=histo.ProjectionX()
histout=histo.ProjectionY()
histo.GetXaxis().SetTitle("cos(theta_in)")
histo.GetYaxis().SetTitle("cos(theta_out)")
c1.cd(1)
histo.Draw("COL Z")
c1.cd(2)
histin.Draw()
c1.cd(3)
histout.Draw()
c1.Update()
raw_input("RET to exit")
