import ROOT
import numpy
import math
from array import array
import sys
import numpy as np
import ast
import glob

plotType=0
#eventNum=0
NEvents=193
gridOn=True
fitOn=True
plotValue=""
plotValue1=-1

if (len(sys.argv)>1):
    eventNum=int(sys.argv[1])
if (len(sys.argv)>2):
    plotType= int(sys.argv[2])
if (len(sys.argv)>3): 
    try:
        plotValue = ast.literal_eval(sys.argv[3])
    except:
        plotValue = sys.argv[3]
if (len(sys.argv)>4):
    try:
        plotValue1 = ast.literal_eval(sys.argv[4])
    except:
        plotValue1 = -1
fileDifs=[]
for rootfile in glob.glob("*allevents.root"):
    fileDifs.append(rootfile[:11])

#fileDifs = ["PTXI______"]

#fitFileName = "PTXIABGaCL_100000steps_1cones.root"
#fitFile = ROOT.TFile(fitFileName)
#bayesTree = ROOT.gDirectory.Get("bayesTree")
#fitFile = "out_n193_step30000.root"
#inFile = "testCC_new.root"

def generateGrid(ntheta,nphi):
    r = 840
    thetagrid = []
    phigrid = []
    xgrid = []
    ygrid = []
    zgrid = []
    for theta in np.arange(0, np.pi, np.pi/ntheta):
        for phi in np.arange(0, 2*np.pi, 2*np.pi/nphi):
            thetagrid.append(theta)
            phigrid.append(phi)
            xgrid.append(r*np.sin(theta)*np.cos(phi))
            ygrid.append(r*np.sin(theta)*np.sin(phi))
            zgrid.append(r*np.cos(theta))
    xgridA = array('d', xgrid)
    ygridA = array('d', ygrid)
    zgridA = array('d', zgrid)
    thetagridA = array('d', thetagrid)
    phigridA = array('d', phigrid)
  
    return xgridA, ygridA, zgridA, thetagridA, phigridA

def spherical (x, y, z):
    r = np.sqrt(x*x+y*y+z*z)
    theta = np.arccos(z/r)
    phiTemp = np.arctan(y/x)
    if (x < 0):
        phiTemp = phiTemp + np.pi
    if (phiTemp < 0):
        phiTemp = phiTemp + 2*np.pi
    phi = phiTemp
    return r, theta, phi


def getPMTData(infile, fitData, eventNumN):
    dataFile = ROOT.TFile(infile)
    tPMT = ROOT.gDirectory.Get("tPMT")
    nEvent = tPMT.GetEntry(eventNumN)
    mr_x = fitData["mr_x"]
    mr_y = fitData["mr_y"]
    mr_z = fitData["mr_z"]
    mt = fitData["mt"]
#    lightSpeed = fitData["mlightSpeed"]
    lightSpeed = 21.5
    x=[]
    y=[]
    z=[]
    phi=[]
    theta=[]
    r=[]
    nPMTs = tPMT.nPMTs
#    print nPMTs
    for i in range(nPMTs):
        Qi_x = tPMT.pmtX[i] - mr_x
        Qi_y = tPMT.pmtY[i] - mr_y
        Qi_z = tPMT.pmtZ[i] - mr_z
        Qi = np.sqrt(Qi_x*Qi_x+Qi_y*Qi_y+Qi_z*Qi_z)
        t_i = tPMT.pmtT[i]
        correctedTime = t_i - mt - Qi/lightSpeed
#        print correctedTime
        if np.abs(correctedTime)<5:
            x.append(tPMT.pmtX[i])
            y.append(tPMT.pmtY[i])
            z.append(tPMT.pmtZ[i])
            rnew, thetanew, phinew = spherical(x[-1],y[-1],z[-1])
            r.append(rnew)
            theta.append(thetanew)
            phi.append(phinew)
    xArray = array("d",x)
    yArray = array("d",y)
    zArray = array("d",z)
    rArray = array("d",r)
    phiArray = array("d",phi)
    thetaArray = array("d",theta)
    return xArray, yArray, zArray, phiArray, thetaArray, len(xArray)


def getFitData(fitFile,fast=True, eventNumN=193):
    nEvent = bayesTree.GetEntry(eventNumN)
    nCones=bayesTree.nCones
    mu_x=[]
    mu_y=[]
    mu_z=[]
    fu_angle=[]
    mconeIntensity=[]
    mflatIntensity=[]
    mleslieIntensity=[]
    munknownIntensity=[]
    for cone in range(nCones):
        mu_x.append(bayesTree.mu_x[cone])
        mu_y.append(bayesTree.mu_y[cone])
        mu_z.append(bayesTree.mu_z[cone])
        fu_angle.append(bayesTree.fu_angle[cone])
        mconeIntensity.append(bayesTree.mconeIntensity[cone])
#    u_x=[]
#    u_y=[]
#    u_z=[]
#    coneIntensity=[]
#    if (nCones>0):
#      u_x_0=[]
#      u_y_0=[]
#      u_z_0=[]
#      coneIntensity_0=[]
#      if (nCones>1):
#        u_x_1=[]
#        u_y_1=[]
#        u_z_1=[]
#        coneIntensity_1=[]
#        if (nCones>2):
#          u_x_2=[]
#          u_y_2=[]
#          u_z_2=[]
#          coneIntensity_2=[]
#          if (nCones>3):
#            u_x_3=[]
#            u_y_3=[]
#            u_z_3=[]
#            coneIntensity_3=[]
#            if (nCones>4):
#              u_x_4=[]
#              u_y_4=[]
#              u_z_4=[]
#              coneIntensity_4=[]
#    flatIntensity=[]
#    _step = []
#    lnP=[]
#    t = []
#    r_x = []
#    r_y = []
#    r_z = []
#    r = []
#    alpha = []
#    beta = []
#    gamma = []
#    trms = []
#    angle = []
#    coneAngle = []
#    lightSpeed = []
#        u_x.append([])
#        u_y.append([])
#        u_z.append([])
#        coneIntensity.append([])
#
#    if (not fast):
#        for step in range(bayesTree.steps):
#            _step.append(bayesTree.step[step])
#            t.append(bayesTree.t[step])
#            r_x.append(bayesTree.r_x[step])
#            r_y.append(bayesTree.r_y[step])
#            r_z.append(bayesTree.r_z[step])
#            r.append(bayesTree.r[step])
#            alpha.append(bayesTree.alpha[step])
#            beta.append(bayesTree.beta[step])
#            gamma.append(bayesTree.gamma[step])
#            trms.append(bayesTree.trms[step])
#            angle.append(bayesTree.angle[step])
#            coneAngle.append(bayesTree.coneAngle[step])
#            lightSpeed.append(bayesTree.lightSpeed[step])
#            lnPtest=bayesTree.lnP[step]
#            if (lnPtest>-1930):
#              lnP.append(lnPtest)
#            flatIntensity.append(bayesTree.flatIntensity[step])
#            if (nCones>0):
#              coneIntensity_0.append(bayesTree.coneIntensity_0[step])
#              u_x_0.append(bayesTree.u_x_0[step])
#              u_y_0.append(bayesTree.u_y_0[step])
#              u_z_0.append(bayesTree.u_z_0[step])
#              if (nCones>1):
#                coneIntensity_1.append(bayesTree.coneIntensity_1[step])
#                u_x_1.append(bayesTree.u_x_1[step])
#                u_y_1.append(bayesTree.u_y_1[step])
#                u_z_1.append(bayesTree.u_z_1[step])
#                if (nCones>2):
#                  coneIntensity_2.append(bayesTree.coneIntensity_2[step])
#                  u_x_2.append(bayesTree.u_x_2[step])
#                  u_y_2.append(bayesTree.u_y_2[step])
#                  u_z_2.append(bayesTree.u_z_2[step])
#                  if (nCones>4):
#                    coneIntensity_4.append(bayesTree.coneIntensity_4[step])
#                    u_x_4.append(bayesTree.u_x_4[step])
#                    u_y_4.append(bayesTree.u_y_4[step])
#                    u_z_4.append(bayesTree.u_z_4[step])
#                    if (nCones>4):
#                      coneIntensity_4.append(bayesTree.coneIntensity_4[step])
#                      u_x_4.append(bayesTree.u_x_4[step])
#                      u_y_4.append(bayesTree.u_y_4[step])
#                      u_z_4.append(bayesTree.u_z_4[step])
##            for cone in range(nCones):
##              u_x[cone].append(bayesTree.u_x[step*nCones + cone])
##              u_y[cone].append(bayesTree.u_y[step*nCones + cone])
##              u_z[cone].append(bayesTree.u_z[step*nCones + cone])
##              coneIntensity[cone].append(bayesTree.coneIntensity[step*nCones + cone])
#              
    outputTree={
#      'steps' : bayesTree.steps,
      'nCones' : bayesTree.nCones,
      'fr_x' : bayesTree.fr_x,
      'fr_y' : bayesTree.fr_y,
      'fr_z' : bayesTree.fr_z,
      'fr' : bayesTree.fr,
      'malpha' : bayesTree.malpha,
      'mbeta' : bayesTree.mbeta,
      'mgamma' : bayesTree.mgamma,
      'mtrms' : bayesTree.mtrms,
      'mangle' : bayesTree.mangle,
      'mr_x' : bayesTree.mr_x,
      'mr_y' : bayesTree.mr_y,
      'mr_z' : bayesTree.mr_z,
      'mr' : bayesTree.mr,
      'mt' : bayesTree.mt,
      'mconeAngle' : bayesTree.mconeAngle,
      'mlightSpeed' : bayesTree.mlightSpeed,
      'sr_x' : bayesTree.sr_x,
      'sr_y' : bayesTree.sr_y,
      'sr_z' : bayesTree.sr_z,
      'sr' : bayesTree.sr,

      'mu_x' : mu_x, 
      'mu_y' : mu_y, 
      'mu_z' : mu_z, 
      'fu_angle' : fu_angle,
      'mconeIntensity' : mconeIntensity,
      'mflatIntensity' : mflatIntensity,
      'mleslieIntensity' : mleslieIntensity,
      'munknownIntensity' : munknownIntensity,
#      'u_x' : u_x,
#      'u_y' : u_y,
#      'u_z' : u_z,
#      'coneIntensity' : coneIntensity,
#      'flatIntensity' : flatIntensity,
#
#      'step' : _step,
#      't' : t,
#      'r_x' : r_x,
#      'r_y' : r_y,
#      'r_z' : r_z,
#      'r' : r,
#      'alpha' : alpha,
#      'beta' : beta,
#      'gamma' : gamma,
#      'trms' : trms,
#      'angle' : angle,
#      'coneAngle' : coneAngle,
#      'lightSpeed' : lightSpeed,
#      'lnP' : lnP,
#
    }
#    if (nCones>0):
#      outputTree.update({'u_x_0':u_x_0, 'u_y_0':u_y_0, 'u_z_0':u_z_0, 'coneIntensity_0':coneIntensity_0})
#    if (nCones>1):
#      outputTree.update({'u_x_1':u_x_1, 'u_y_1':u_y_1, 'u_z_1':u_z_1, 'coneIntensity_1':coneIntensity_1})
#    if (nCones>2):
#      outputTree.update({'u_x_2':u_x_2, 'u_y_2':u_y_2, 'u_z_2':u_z_2, 'coneIntensity_2':coneIntensity_2})
#    if (nCones>3):
#      outputTree.update({'u_x_3':u_x_3, 'u_y_3':u_y_3, 'u_z_3':u_z_3, 'coneIntensity_3':coneIntensity_3})
#    if (nCones>4):
#      outputTree.update({'u_x_4':u_x_4, 'u_y_4':u_y_4, 'u_z_4':u_z_4, 'coneIntensity_4':coneIntensity_4})
#
    return outputTree


def plotCone(fitData):
    mr_x = fitData["mr_x"]
    mr_y = fitData["mr_y"]
    mr_z = fitData["mr_z"]
    mconeAngle = fitData["mconeAngle"]
    axArrayA=[]
    ayArrayA=[]
    azArrayA=[]
    cone_xArrayA=[]
    cone_yArrayA=[]
    cone_zArrayA=[]
    outerRing_xArrayA=[]
    outerRing_yArrayA=[]
    outerRing_zArrayA=[]
    outerRing_thetaArrayA=[]
    outerRing_phiArrayA=[]
    for cone in range(fitData["nCones"]):
        mu_x0 = fitData["mu_x"][cone]
        mu_y0 = fitData["mu_y"][cone]
        mu_z0 = fitData["mu_z"][cone]
        mu_x = mu_x0 / np.sqrt(mu_x0*mu_x0+mu_y0*mu_y0+mu_z0*mu_z0)
        mu_y = mu_y0 / np.sqrt(mu_x0*mu_x0+mu_y0*mu_y0+mu_z0*mu_z0)
        mu_z = mu_z0 / np.sqrt(mu_x0*mu_x0+mu_y0*mu_y0+mu_z0*mu_z0)
   
    #axis starts now
        axisX = []
        axisY = []
        axisZ = []
        ax=mr_x
        ay=mr_y
        az=mr_z
        count=0
        while(ax*ax+ay*ay+az*az<840*840 and count<1680):
            axisX.append(ax)
            axisY.append(ay)
            axisZ.append(az)
            ax=ax+mu_x
            ay=ay+mu_y
            az=az+mu_z
            count=count+1
        axArray=array("d",axisX)
        ayArray=array("d",axisY)
        azArray=array("d",axisZ)

        # cone starts now
        perpToU_x1 = -mu_y/np.sqrt(mu_x*mu_x+mu_y*mu_y)*np.sin(mconeAngle)
        perpToU_y1 = mu_x/np.sqrt(mu_x*mu_x+mu_y*mu_y)*np.sin(mconeAngle)
        #perpToU_z1 = 0
        perpToU_x2 = -mu_z*perpToU_y1 # + mu_y*perpToU_z1
        perpToU_y2 = mu_z*perpToU_x1 # - mu_x*perpToU_z1
        perpToU_z2 = mu_x*perpToU_y1 - mu_y*perpToU_x1
        shortU_x = mu_x*np.cos(mconeAngle)
        shortU_y = mu_y*np.cos(mconeAngle)
        shortU_z = mu_z*np.cos(mconeAngle)
#        print shortU_x*perpToU_x1+shortU_y*perpToU_y1
#        print shortU_x*perpToU_x2+shortU_y*perpToU_y2+shortU_z*perpToU_z2
#        print perpToU_x1*perpToU_x2+perpToU_y1*perpToU_y2
#        print mu_x*mu_x+mu_y*mu_y+mu_z*mu_z 
#        print (shortU_x*shortU_x+shortU_y*shortU_y+shortU_z*shortU_z)/(np.cos(mconeAngle)*np.cos(mconeAngle))
#        print (perpToU_x1*perpToU_x1+perpToU_y1*perpToU_y1)/(np.sin(mconeAngle)*np.sin(mconeAngle))
#        print (perpToU_x2*perpToU_x2+perpToU_y2*perpToU_y2+perpToU_z2*perpToU_z2)/(np.sin(mconeAngle)*np.sin(mconeAngle))
        ringNum = 100
        coneStart_x=[]
        coneStart_y=[]
        coneStart_z=[]
        cone_x = []
        cone_y = []
        cone_z = []
        outerRing_x = []
        outerRing_y = []
        outerRing_z = []
        outerRing_theta=[]
        outerRing_phi=[]
        for i in range(ringNum):
            coneStart_x.append(shortU_x+perpToU_x1*np.sin(2*np.pi/ringNum*i)+perpToU_x2*np.cos(2*np.pi/ringNum*i))
            coneStart_y.append(shortU_y+perpToU_y1*np.sin(2*np.pi/ringNum*i)+perpToU_y2*np.cos(2*np.pi/ringNum*i))
            coneStart_z.append(shortU_z+perpToU_z2*np.cos(2*np.pi/ringNum*i))
            conePos_x = mr_x
            conePos_y = mr_y
            conePos_z = mr_z
            while (conePos_x*conePos_x+conePos_y*conePos_y+conePos_z*conePos_z <840*840):
                conePos_x = conePos_x + coneStart_x[i]
                conePos_y = conePos_y + coneStart_y[i]
                conePos_z = conePos_z + coneStart_z[i]
                cone_x.append(conePos_x)
                cone_y.append(conePos_y)
                cone_z.append(conePos_z)
            outerRing_x.append(conePos_x)
            outerRing_y.append(conePos_y)
            outerRing_z.append(conePos_z)
            oRr, oRtheta, oRphi = spherical(conePos_x, conePos_y, conePos_z)
            outerRing_theta.append(oRtheta)
            outerRing_phi.append(oRphi)
        cone_xArray = array("d",cone_x)
        cone_yArray = array("d",cone_y)
        cone_zArray = array("d",cone_z)
        outerRing_xArray=array("d",outerRing_x)
        outerRing_yArray=array("d",outerRing_y)
        outerRing_zArray=array("d",outerRing_z)
        outerRing_thetaArray=array("d",outerRing_theta)
        outerRing_phiArray=array("d",outerRing_phi)
#    Ax=axArray[-1]-mr_x
#    Ay=axArray[-1]-mr_y
#    Az=axArray[-1]-mr_z
#    lengthAxis=np.sqrt(Ax*Ax+Ay*Ay+Az*Az)
#    for i in range(len(outerRing_x)):
#        oRx=outerRing_x[i]-mr_x
#        oRy=outerRing_y[i]-mr_y
#        oRz=outerRing_z[i]-mr_z
#        lengthOR=np.sqrt((oRx)*(oRx)+(oRy)*(oRy)+(oRz)*(oRz))
#        cosine = (Ax*oRx+Ay*oRy+Az*oRz)/(lengthOR*lengthAxis)
#        print lengthAxis, lengthOR, cosine
#        print cosine
        axArrayA.append(axArray)
        ayArrayA.append(ayArray)
        azArrayA.append(azArray)
        cone_xArrayA.append(cone_xArray)
        cone_yArrayA.append(cone_yArray)
        cone_zArrayA.append(cone_zArray)
        outerRing_xArrayA.append(outerRing_xArray)
        outerRing_yArrayA.append(outerRing_yArray)
        outerRing_zArrayA.append(outerRing_zArray)
        outerRing_thetaArrayA.append(outerRing_thetaArray)
        outerRing_phiArrayA.append(outerRing_phiArray)
    outputCone={
        "axArrayA" : axArrayA,
        "ayArrayA" : ayArrayA,
        "azArrayA" : azArrayA,
        "cone_xArrayA" : cone_xArrayA,
        "cone_yArrayA" : cone_yArrayA,
        "cone_zArrayA" : cone_zArrayA,
        "outerRing_xArrayA" : outerRing_xArrayA,
        "outerRing_yArrayA" : outerRing_yArrayA,
        "outerRing_zArrayA" : outerRing_zArrayA,
        "outerRing_thetaArrayA" : outerRing_thetaArrayA,
        "outerRing_phiArrayA" : outerRing_phiArrayA
    } 

    return outputCone

ROOT.gROOT.Reset()
c0=ROOT.TCanvas("c1","c1",600, 600)

#fileDifs = ["PTXI______"]
files=[
#"data/193ev_lesliesurf_lesliecomparison_b14/PTDXI_FV____100000steps_1cones_allevents.root",
#"data/193ev_leslie_lesliecomparison_b14/PTDXI_FV____100000steps_1cones_allevents.root",
"data/193ev_bulkleslies_b14/PTDXI_FV____100000steps_1cones_allevents.root",
"data/193ev_surfleslies_b14/PTDXI_FV____100000steps_1cones_allevents.root",
#"data/193ev_surfleslies/PTDXI_FV____100000steps_1cones_allevents.root",
#"data/3000ev_8B_lesliecomparison_nob14/PTDXI_FV____100000steps_1cones_allevents.root",
#"data/3000ev_8B_lesliecomparison_nob14_cosAE/PTDXI_FV____100000steps_1cones_allevents.root",
"data/495ev_8B_b14/PTDXI_FV____100000steps_1cones_allevents.root",
]
cone_flat=[]
colors=[ROOT.kRed,ROOT.kBlue,ROOT.kGreen,ROOT.kBlack,ROOT.kMagenta]
for eachfile in files:
    cone_flat.append(ROOT.TGraph())
#    allFileName = dif + "_100000steps_3cones_allevents.root"
    allFile = ROOT.TFile(eachfile, "read")
    allEventTree = ROOT.gDirectory.Get("allEventTree")
    nEvent = allEventTree.GetEntry(0)
    variableDict={
        "mt" : {"data" : allEventTree.mt, "dmin":0, "dmax":200},
        "nPMTsTotal" : {"data" : allEventTree.nPMTsTotal, "dmin":0, "dmax":200},
        "nPMTsGood" : {"data" : allEventTree.nPMTsGood, "dmin":0, "dmax":200},
        "trueE" : {"data" : allEventTree.trueE, "dmin":0, "dmax":20},
        "mr" : {"data" : allEventTree.mr, "dmin":0, "dmax":900},
        "mr_x" : {"data" : allEventTree.mr_x, "dmin":-600, "dmax":600},
        "fr" : {"data" : allEventTree.fr, "dmin":600, "dmax":600},
        "fu_angle_0" : {"data" : allEventTree.fu_angle_0, "dmin":0, "dmax":3},
        "fu_cos_0": {"data" : np.cos(allEventTree.fu_angle_0),"dmin":-1,"dmax":1},
        "fu_angle_1" : {"data" : allEventTree.fu_angle_1, "dmin":0, "dmax":3},
        "fu_cos_1": {"data" : np.cos(allEventTree.fu_angle_1),"dmin":-1,"dmax":1},
        "fu_angle_2" : {"data" : allEventTree.fu_angle_2, "dmin":0, "dmax":3},
        "fu_cos_2": {"data" : np.cos(allEventTree.fu_angle_2),"dmin":-1,"dmax":1},
        "mu_x_0" : {"data" : allEventTree.mu_x_0, "dmin":-1, "dmax":1},
        "mu_y_0" : {"data" : allEventTree.mu_y_0, "dmin":-1, "dmax":1},
        "mu_z_0" : {"data" : allEventTree.mu_z_0, "dmin":-1, "dmax":1},
        "fr_x" : {"data" : allEventTree.fr_x, "dmin":-100, "dmax":100},
        "fr_y" : {"data" : allEventTree.fr_y, "dmin":-100, "dmax":100},
        "mangle" : {"data" : allEventTree.mangle, "dmin":0, "dmax":60},
        "mconeAngle" : {"data" : allEventTree.mconeAngle, "dmin":0, "dmax":1},
        "mlightSpeed" : {"data" : allEventTree.mlightSpeed, "dmin":0, "dmax":30},
        "dr" : {"data" : np.sqrt(np.power(allEventTree.fr_x,2)+np.power(allEventTree.fr_y,2)+np.power(allEventTree.fr_z,2)), "dmin":0, "dmax":300},
        "malpha" : {"data" : allEventTree.malpha, "dmin":0, "dmax":1},
        "mtrms" : {"data" : allEventTree.mtrms, "dmin":0, "dmax":15},
        "mconeIntensity_0" : {"data" : allEventTree.mconeIntensity_0, "dmin":0, "dmax":1},
        "mflatIntensity" : {"data" : allEventTree.mflatIntensity, "dmin":0, "dmax":1},
        "mleslieIntensity" : {"data" : allEventTree.mleslieIntensity, "dmin":0, "dmax":1},
        "munknownIntensity" : {"data" : allEventTree.munknownIntensity, "dmin":0, "dmax":1},
        "beta14" : {"data" : allEventTree.beta14, "dmin":0, "dmax":2}
    }
    fxdivmx=np.zeros(len(variableDict["mr_x"]["data"]),dtype=float)
    flat_div_cone_plus_flat=np.zeros(len(variableDict["mr_x"]["data"]),dtype=float)
    coneness=np.zeros(len(variableDict["mr_x"]["data"]),dtype=float)
    leslieness=np.zeros(len(variableDict["mr_x"]["data"]),dtype=float)
    for index in range(len(variableDict["mr_x"]["data"])):
        fxdivmx[index]=(variableDict["fr_x"]["data"][index]/variableDict["mr_x"]["data"][index])
        flat_div_cone_plus_flat[index]=allEventTree.mflatIntensity[index]/(1-allEventTree.munknownIntensity[index])
#        coneness[index]=1-allEventTree.mflatIntensity[index]-allEventTree.munknownIntensity[index]
        coneness[index]=np.sum([allEventTree.mconeIntensity_0[index],allEventTree.mconeIntensity_1[index],allEventTree.mconeIntensity_2[index]])
        leslieness[index]=allEventTree.mleslieIntensity[index]
        if len(cone_flat)==2 and variableDict["mr"]["data"][index]<400:
            pass
#            print index
    variableDict.update({"fr_x_div_mr_x" : {"data" : fxdivmx, "dmin":-5, "dmax":5}})
    variableDict.update({"flat_div_cone_plus_flat" : {"data" : flat_div_cone_plus_flat, "dmin":0, "dmax":1}})
    variableDict.update({"coneness" : {"data" : coneness, "dmin":0, "dmax":1}})
    variableDict.update({"leslieness" : {"data" : leslieness, "dmin":0, "dmax":1}})
#    for variable in ["mt","nPMTsTotal","nPMTsGood","mr","malpha","mr_x","fr_x_div_mr_x","fr","fr_y","fu_angle_0","fu_cos_0","fu_angle_1","fu_cos_1","fu_angle_2","fu_cos_2","mu_x_0","mu_y_0","mu_z_0","fr_x","mangle","mconeAngle","mlightSpeed","mconeIntensity_0","mflatIntensity","munknownIntensity","mlightSpeed","dr","mtrms","flat_div_cone_plus_flat","coneness"]:
#        title = variable + " histogram " + dif
#        nbins=100
#        variablehist=ROOT.TH1F("varhist", title,nbins,variableDict[variable]["dmin"],variableDict[variable]["dmax"])
#        variablehist.Clear()
#        variablehist.SetTitle(title)
#        for item in variableDict[variable]["data"]:
#            variablehist.Fill(item)
#        variablehist.DrawNormalized()  
#        if variable is "mr":
#            xlinetemp = np.arange(0, 606, 6)
#            xline = array("d", xlinetemp)
#            yline = array("d", 3.0*(3.0/2.0)/nbins*np.power(xlinetemp/600., 2))
#            print np.sum(yline)
#            line = ROOT.TGraph(len(xline), xline, yline)
#            line.Draw("SAME C")      
#        c0.Update()
#        c0.SaveAs(title+".png")
    cone_flat[-1] = ROOT.TGraph(len(variableDict["coneness"]["data"]),array("d",variableDict["coneness"]["data"]),array("d",variableDict["beta14"]["data"]))
    for i in np.arange(len(variableDict["mr"]["data"])-1,0,-1):
        if variableDict["trueE"]["data"][i]<9 and variableDict["trueE"]["data"][i]>1:
            cone_flat[-1].RemovePoint(int(i))
        elif variableDict["mconeIntensity_0"]["data"][i]<.4:
            print i
    cone_flat[-1].SetMarkerStyle(20)
    cone_flat[-1].GetXaxis().SetLimits(-0,1.0)
    cone_flat[-1].GetXaxis().SetTitle("cone intensity")
    cone_flat[-1].GetYaxis().SetRangeUser(-0.1,1.0)
    cone_flat[-1].GetYaxis().SetTitle("beta14")
    cone_flat[-1].SetMarkerColor(colors[len(cone_flat)-1])
    if len(cone_flat) is 1:
        cone_flat[-1].Draw("AP")
    else:
        cone_flat[-1].Draw("P SAME")
    c0.Update()
    allFile.Close()
c0.SaveAs("beta14.png")
     





#if (plotType==1):
#    fitData = getFitData(fitFile, fast=True, eventNumN=eventNum)
#    xArray, yArray, zArray, phiArray, thetaArray, nPMTs = getPMTData("testCC_new.root",fitData, eventNum)
##    grid2D.SetMarkerStyle(20)
#    c1=ROOT.TCanvas("c1","c1",600, 600)
#    pmt2D = ROOT.TGraph2D(len(xArray),xArray,yArray,zArray)
#    pmt2D.SetMarkerStyle(20)
#    pmt2D.Draw("AP")
#    if (gridOn>0):
#        xgrid,ygrid,zgrid,thetagrid,phigrid = generateGrid(193,193)
#        grid2D = ROOT.TGraph2D(len(xgrid),xgrid,ygrid,zgrid)
#        grid2D.SetMarkerColor(ROOT.kRed)
#        grid2D.Draw("P SAME")
#    if (fitOn>0):
#        coneData = plotCone(fitData)
#        cone2D=[]
#        axis2D=[]
#        ring2D=[]
#        for cone in range(fitData["nCones"]):
#            cone2D.append(ROOT.TGraph2D(len(coneData["cone_xArrayA"][cone]),coneData["cone_xArrayA"][cone],coneData["cone_yArrayA"][cone],coneData["cone_zArrayA"][cone]))
#            cone2D[cone].SetMarkerColor(ROOT.kGreen+int(fitData["mconeIntensity"][cone]*5))
#            cone2D[cone].Draw("P SAME")
#            axis2D.append(ROOT.TGraph2D(len(coneData["axArrayA"][cone]),coneData["axArrayA"][cone],coneData["ayArrayA"][cone],coneData["azArrayA"][cone]))
#            axis2D[cone].SetMarkerColor(ROOT.kBlue+fitData['nCones']-cone-1)
#            axis2D[cone].Draw("P SAME")
#            ring2D.append(ROOT.TGraph2D(len(coneData["outerRing_xArrayA"][cone]),coneData["outerRing_xArrayA"][cone],coneData["outerRing_yArrayA"][cone], coneData["outerRing_zArrayA"][cone]))
#            ring2D[cone].SetMarkerColor(ROOT.kGreen+int(fitData["mconeIntensity"][cone]*5))
#            ring2D[cone].SetMarkerStyle(20)
#            ring2D[cone].Draw("P SAME")
#        
#    c1.Update()
#        
#if (plotType == 0 or plotType == 2):
#    c0=ROOT.TCanvas("c1","c1",600, 600)
#    fitData = getFitData(fitFile, fast=True, eventNumN=eventNum)
#    xArray, yArray, zArray, phiArray, thetaArray, nPMTs = getPMTData("testCC_new.root",fitData, eventNum)
#    allGraph = ROOT.TMultiGraph()
#    print nPMTs
#    pmt1D = ROOT.TGraph(nPMTs,phiArray,thetaArray)
#    pmt1D.SetMarkerStyle(20)
#    if (fitOn):
#        cone1D=[]
#        axis1D=[]
#        ring1D=[]
#        coneData = plotCone(fitData)
#        for cone in range(fitData["nCones"]):
#            ring1D.append(ROOT.TGraph(len(coneData["outerRing_thetaArrayA"][cone]),coneData["outerRing_phiArrayA"][cone],coneData["outerRing_thetaArrayA"][cone]))
#            ring1D[cone].SetMarkerStyle(20)
#            ring1D[cone].SetMarkerColor(ROOT.kGreen+int(fitData["mconeIntensity"][cone]*5))
##            ring1D[cone].SetMarkerColor(ROOT.kGreen+fitData["nCones"]-cone-1)
#            allGraph.Add(ring1D[cone])
#    allGraph.Add(pmt1D)
#    allGraph.Draw("AP")
#    c0.Update()
#    c0.SaveAs("plots/ev"+str(eventNum)+"_2dPlot.png")
#
#if (plotType == 2):
#    c2=ROOT.TCanvas("c2","c2",600,800)
#    c2.Divide(1,2)
#    fitData = getFitData(fitFile, fast=False, eventNumN=eventNum)
#    print "checking lnP.  Steps:", len(fitData["lnP"])
#    oknum=0
#    for i in range(len(fitData["lnP"])):
#        if (math.isnan(fitData["lnP"][i])):
#            print i, "nan"
#        if (float('Inf')==fitData["lnP"][i]):
#            print i, "inf"
#        if (float('Inf')==fitData["lnP"][i]):
#            print i, "-inf"
#        if (fitData["lnP"][i] < 19300 and fitData["lnP"][i]>-19300):
#            oknum += 1
#        else:
#            print i, fitData["lnP"][i], "outside range"
#    print "ok:", oknum
#    xArray, yArray, zArray, phiArray, thetaArray, nPMTs = getPMTData("testCC_new.root",fitData, eventNum)
#    step = array('d', fitData["step"])
#    lightSpeed = array('d', fitData["lightSpeed"])
#    nCones=fitData["nCones"]
#    saveConeIndex=False
#    if plotValue == "all":
#        plotValueArray = ["r_x", "r_y", "r_z","flatIntensity","lnP", "coneAngle", "lightSpeed"]
##       plotValueArray = plotValueArray + ["u_x","u_y","u_z","coneIntensity"]
#        if (nCones>0):
#          plotValueArray = plotValueArray + ["u_x_0","u_y_0","u_z_0","coneIntensity_0"]
#        if (nCones>1):
#          plotValueArray = plotValueArray + ["u_x_1","u_y_1","u_z_1","coneIntensity_1"]
#        if (nCones>2):
#          plotValueArray = plotValueArray + ["u_x_2","u_y_2","u_z_2","coneIntensity_2"]
#        if (nCones>3):
#          plotValueArray = plotValueArray + ["u_x_3","u_y_3","u_z_3","coneIntensity_3"]
#        if (nCones>4):
#          plotValueArray = plotValueArray + ["u_x_4","u_y_4","u_z_4","coneIntensity_4"]
#    elif (type(plotValue)==str):
#        plotValueArray = [plotValue]
#    if type(plotValue1) is int:
#        plotValue1Array = [plotValue1]
#    else:
#        plotValue1Array = plotValue1
#    for plotValuei in plotValueArray:
#        saveConeIndex=False
#        if plotValuei in ["u_x", "u_y", "u_z"]:
#            saveConeIndex=True
#            if (plotValue1Array[0] < 0):
#                plotValueConesThisVariable = range(fitData["nCones"])
#            else:
#                plotValueConesThisVariable = plotValue1Array
#        elif plotValuei in ["coneIntensity"]:
#            saveConeIndex=True
#            if (plotValue1Array[0] < 0):
#                plotValueConesThisVariable = range(fitData["nCones"])
#            else:
#                plotValueConesThisVariable = plotValue1Array
#        else:
#            plotValueConesThisVariable = [-1]
#        if (plotValuei == ""):
#            plotValuei="lnP"
#        for plotValueCone in plotValueConesThisVariable:
#            if plotValuei in ["u_x", "u_y", "u_z", "coneIntensity"]:
#                plotArray = array('d', fitData[plotValuei][plotValueCone])
#            else:
#                plotArray = array('d', fitData[plotValuei])
#
#            MCMCplot = ROOT.TGraph(fitData["steps"], step, plotArray)
#            histo=ROOT.TH1F("histo","histo",193,np.min(plotArray),np.max(plotArray))
#            for datum in plotArray:
#                histo.Fill(datum)
#           
#            c2.cd(1)
#            MCMCplot.Draw("AP")
#            if saveConeIndex:
#                MCMCplot.SetTitle("Event " + str(eventNum)+ " " + plotValuei + " " + str(plotValueCone))
#            else:
#                MCMCplot.SetTitle("Event " + str(eventNum)+ " " + plotValuei)
#            c2.cd(2)
#            histo.Draw()
#            c2.Update()
#            if saveConeIndex:
#                c2.SaveAs("plots/ev"+str(eventNum)+plotValuei+str(plotValueCone)+".png")
#            else:
#                c2.SaveAs("plots/ev"+str(eventNum)+plotValuei+".png")
#
#raw_input("RET to exit")
