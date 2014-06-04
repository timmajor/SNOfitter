#!/usr/bin/python

import os
import subprocess
import glob
import shutil
import fileinput
import re

# A detailed description so I can remember what run was for
runDescription='''
Compare beta14 to coneLikelihood.
8B events
This time put back in cos(areaelement)
Leslie time distribution piecewise on mc
Leslie position center of sqrt at -.07
Fixed light speed, alpha, coneAngle.
'''
# A short version for directory name
shortRunDescription="testingallPMT"
numberofsteps="100000"
numberofevents="1"
numberofcones="1"
#inputfile="AVsurface_readable_b14.root"
#inputfile="AVbulk_readable_b14.root"
inputfile="8B_readable_pmtN.root"
#inputfile="testCC_new.root"

needDoubles=True
if inputfile is "testCC_new.root":
    needDoubles=False

os.chdir("/home/tim/sno_code/fittingcode")
shutil.copy("BayesGaussianFloat4master.C","BayesGaussianFloat4.C")
shutil.copy("conemacromaster.mac","conemacro.mac")

    

for line in fileinput.input("conemacro.mac", inplace=True):
    line=line.replace("<NUMBEROFEVENTS>",numberofevents)
    line=line.replace("<INPUTFILE>",inputfile)
    line=line.replace("\n","")
    print line
        
for line in fileinput.input("BayesGaussianFloat4.C", inplace=True):
    if needDoubles:
        line=line.replace("Float_t", "Double_t")
    line=line.replace("<NUMBEROFSTEPS>",numberofsteps)
    line=line.replace("<NUMBEROFCONES>",numberofcones)
    line=line.replace("\n","")
    print line
        
subprocess.call("root -l -q conemacro.mac",shell=True)

newDir=""
if os.path.exists("data/"+shortRunDescription):
    increment=0
    while os.path.exists("data/"+shortRunDescription+str(increment)):
        increment+=1
    newDir="data/"+shortRunDescription+str(increment)
else:
    newDir="data/"+shortRunDescription
os.makedirs(newDir)

os.chdir("data/")
myRootFile=""
for rootfile in glob.glob("*.root"):
    myRootFile=rootfile
    shutil.move(rootfile, "../"+newDir)

tempstring1=re.sub('steps_'+str(numberofcones)+'cones.root','',myRootFile)
numberofsteps=re.sub('[^0-9]*_','',tempstring1)

os.chdir("../masterscripts")
for eachfile in glob.glob("*"):
    shutil.copy(eachfile, "../"+newDir)
    # This adjusts the scripts to have the right number of events
    for line in fileinput.input("../"+newDir+"/"+eachfile, inplace=True):
        line=line.replace("<NUMBEROFEVENTS>",numberofevents)
        line=line.replace("<NUMBEROFSTEPS>",numberofsteps)
        line=line.replace("<NUMBEROFCONES>",numberofcones)
        line=line.replace("<INPUTFILE>",inputfile)
        line=line.replace("\n","")
        print line
        

os.chdir("../"+newDir)
shutil.move("../../BayesGaussianFloat4.C","./")
shutil.move("../../conemacro.mac","./")

rundesc=open("rundescription.txt", 'w')
rundesc.write(runDescription)
rundesc.close()

os.makedirs("1event")

subprocess.call("python allEvents.py",shell=True)
#subprocess.call("python plotAll.py",shell=True)
subprocess.call("./manyplots.sh",shell=True)
#subprocess.call("python plotFit.py 0 0",shell=True)
