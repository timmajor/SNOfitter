import ROOT
import math
from array import array
import sys
import numpy as np
import ast
import glob
import time


dataFile = ROOT.TFile("/home/tim/sno_code/fittingcode/testCC_new.root")
tPMT = ROOT.gDirectory.Get("tPMT")

def getPMTData(eventNum, lightSpeed):
    nEvent = tPMT.GetEntry(eventNum)
#    lightSpeed = 21.5
    nPMTs = tPMT.nPMTs
    t=[]
    for i in range(nPMTs):
        Qi_x = tPMT.pmtX[i] - tPMT.trueX
        Qi_y = tPMT.pmtY[i] - tPMT.trueY
        Qi_z = tPMT.pmtZ[i] - tPMT.trueZ
        Qi = np.sqrt(Qi_x*Qi_x+Qi_y*Qi_y+Qi_z*Qi_z)
        correctedTime = tPMT.pmtT[i] - tPMT.trueT - Qi/lightSpeed
        t.append(correctedTime)
    return t, nPMTs


ROOT.gROOT.Reset()

c1 = ROOT.TCanvas("c1")

f1 = ROOT.TF1("f","[2]/(sqrt(3.14159)*[1])*exp(-(x-[0])*(x-[0])/(2*[1]*[1]))+[3]",50, 250)
#weirdevents=[]
#trmshist = ROOT.TH1F("trmshist","trmshist",450,0,15)
trmss=[]
trmsrmss=[]
speeds=[]
lightspeed=22.5

#for lightspeed in np.arange(20,25,0.5):
for lightspeed in [22.0]:
#    for event in range(tPMT.GetEntries()):
    for event in range(60):
        thist = ROOT.TH1F("thist","thist",150,50,150)
        t, nPMTs = getPMTData(event, lightspeed)
    
        for element in t:
            thist.Fill(element)
    #    thist.Draw()
    
#        f1.SetParameters(115,3,.3,.0001)
        peak=thist.GetBinCenter(thist.GetMaximumBin())
        f1.SetParameters(peak,2,3,.0001)
        f1.SetParLimits(1,0.5,3)
        f1.SetParLimits(2,0.01,1000)
        f1.SetParLimits(3,0,1)
    
        thist.Draw()
        thist.Fit("f","L","",peak-10,peak+10)
        trms = f1.GetParameter(1)
    #    trmss.append(trms)
        c1.Update()
        time.sleep(1)
    #    print trms
#        trmshist.Fill(trms)
        del thist
    #print trmss
#    trmshist.Draw()
#    f1.SetParameters(2.5,.6,10,.0001)
#    f1.SetParameters(2,.5,4000,.0001)
#    trmshist.Fit("f","L")
#    trmss.append(f1.GetParameter(0))
#    trmsrmss.append(f1.GetParameter(1))
#    speeds.append(lightspeed)

#print speeds, trmss, trmsrmss

#speedArray=array("d",speeds)
#trmsArray=array("d",trmss)
#trmserrorArray=array("d",trmsrmss)
#speederrorArray=array("d",np.zeros(10))
#tg1=ROOT.TGraphErrors(len(speeds), speedArray,trmsArray,speederrorArray,trmserrorArray)

#tg1.Draw("A*")
#if (plotType==1):
#    fitData = getFitData(fitFile, fast=True)
#    xArray, yArray, zArray, phiArray, thetaArray, nPMTs = getPMTData("/home/tim/sno_code/fittingcode/testCC_new.root",fitData)
##    grid2D.SetMarkerStyle(20)
#    c1=ROOT.TCanvas("c1","c1",600, 600)
#    pmt2D = ROOT.TGraph2D(len(xArray),xArray,yArray,zArray)
#    pmt2D.SetMarkerStyle(20)
#    pmt2D.Draw("AP")
#    if (gridOn>0):
#        xgrid,ygrid,zgrid,thetagrid,phigrid = generateGrid(100,100)
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
#if (plotType == 0):
#    c0=ROOT.TCanvas("c1","c1",600, 600)
#    fitData = getFitData(fitFile, fast=True)
#    xArray, yArray, zArray, phiArray, thetaArray, nPMTs = getPMTData("/home/tim/sno_code/fittingcode/testCC_new.root",fitData)
#    allGraph = ROOT.TMultiGraph()
##    print nPMTs
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
#    c0.SaveAs("1event/ev"+str(eventNum)+"_2dPlot.png")
#
#if (plotType == 2):
#    c2=ROOT.TCanvas("c2","c2",600,800)
#    c2.Divide(1,2)
#    fitData = getFitData(fitFile, fast=False)
#    print "checking lnP.  Steps:", len(fitData["lnP"])
#    oknum=0
#    for i in range(len(fitData["lnP"])):
#        if (math.isnan(fitData["lnP"][i])):
#            print i, "nan"
#        if (float('Inf')==fitData["lnP"][i]):
#            print i, "inf"
#        if (float('Inf')==fitData["lnP"][i]):
#            print i, "-inf"
#        if (fitData["lnP"][i] < 1000 and fitData["lnP"][i]>-1000):
#            oknum += 1
#        else:
#            print i, fitData["lnP"][i], "outside range"
#    print "ok:", oknum
#    xArray, yArray, zArray, phiArray, thetaArray, nPMTs = getPMTData("/home/tim/sno_code/fittingcode/testCC_new.root",fitData)
#    step = array('d', fitData["step"])
#    lightSpeed = array('d', fitData["lightSpeed"])
#    nCones=fitData["nCones"]
#    saveConeIndex=False
#    if plotValue == "all":
#        plotValueArray = ["r_x", "r_y", "r_z","alpha","flatIntensity","unknownIntensity","lnP", "coneAngle", "lightSpeed","trms","t"]
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
#            histo=ROOT.TH1F("histo","histo",50,np.min(plotArray),np.max(plotArray))
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
#                c2.SaveAs("1event/ev"+str(eventNum)+plotValuei+str(plotValueCone)+".png")
#            else:
#                c2.SaveAs("1event/ev"+str(eventNum)+plotValuei+".png")
raw_input("RET to exit")

