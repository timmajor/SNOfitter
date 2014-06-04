#include "TROOT.h"
#include "TMath.h"
#include "TObject.h"
#include "TH1F.h"
#include "TFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include <cstring>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

Double_t* getPMTXs(){
  streampos size;
  Double_t * memblock;
  ifstream file ("pmtdata/x.pmt", ios::in|ios::binary|ios::ate);
  if (file.is_open()){
    size = file.tellg();
    memblock = new Double_t [size/4];
    file.seekg (0, ios::beg);
    file.read (reinterpret_cast<char*> (memblock), size);
    file.close();
  }
  return memblock;
}

Double_t* getPMTYs(){
  streampos size;
  Double_t * memblock;
  ifstream file ("pmtdata/y.pmt", ios::in|ios::binary|ios::ate);
  if (file.is_open()){
    size = file.tellg();
    memblock = new Double_t [size/4];
    file.seekg (0, ios::beg);
    file.read (reinterpret_cast<char*>(memblock), size);
    file.close();
  }
  return memblock;
}

Double_t* getPMTZs(){
  streampos size;
  Double_t * memblock;
  ifstream file ("pmtdata/z.pmt", ios::in|ios::binary|ios::ate);
  if (file.is_open()){
    size = file.tellg();
    memblock = new Double_t [size/4];
    file.seekg (0, ios::beg);
    file.read (reinterpret_cast<char*>(memblock), size);
    file.close();
  }
  return memblock;
}

Int_t* getPMTTypes(){
  streampos size;
  Int_t * memblock;
  ifstream file ("pmtdata/n.pmt", ios::in|ios::binary|ios::ate);
  if (file.is_open()){
    size = file.tellg();
    memblock = new Int_t [size/4];
    file.seekg (0, ios::beg);
    file.read (reinterpret_cast<char*>(memblock), size);
    file.close();
  }
  return memblock;
}

void testfunc(){
  cout<< "tested" <<endl;
}

Double_t GetT0(Int_t *nPMTs, Int_t *pmtN, Int_t *pmtType, Float_t *pmtT, Float_t *pmtX, Float_t *pmtY, Float_t *pmtZ, Double_t x, Double_t y, Double_t z) {
  TH1F h_times("h_times","h_times",80,0,400);
  for (Int_t i=0; i<*nPMTs; i++) {
    Double_t Qi_x = pmtX[pmtN[i]] - x;
    Double_t Qi_y = pmtY[pmtN[i]] - y;
    Double_t Qi_z = pmtZ[pmtN[i]] - z;
    Double_t Qi = sqrt(Qi_x*Qi_x+Qi_y*Qi_y+Qi_z*Qi_z);
    Double_t t_i = pmtT[i];
    Bool_t goodPMT = (t_i>50 && t_i<300);
    if (goodPMT and pmtType[pmtN[i]]==1) {
      h_times.Fill(t_i - Qi/21.5);
    }
  }
  cout<<h_times.GetBinCenter(h_times.GetMaximumBin())<<endl;;
  return h_times.GetBinCenter(h_times.GetMaximumBin());
}


void GetNewStep(Double_t t, Double_t r_x, Double_t r_y, Double_t r_z,
		Double_t *coneIntensity, Double_t flatIntensity, Double_t leslieIntensity, Double_t unknownIntensity, Double_t *u_x, Double_t *u_y, Double_t *u_z,
		Double_t alpha, Double_t trms, 
		Double_t coneAngle, Double_t lightSpeed, 
		Double_t &new_t, Double_t &new_r_x, Double_t &new_r_y, Double_t &new_r_z,
		Double_t *new_coneIntensity, Double_t &new_flatIntensity, Double_t &new_leslieIntensity, Double_t &new_unknownIntensity, Double_t *new_u_x, Double_t *new_u_y, Double_t *new_u_z,
		Double_t &new_alpha, Double_t &new_trms, 
		Double_t &new_coneAngle, Double_t &new_lightSpeed, Int_t nCones,
                Bool_t *floating) {

//  Bool_t useT = floating[0];
//  Bool_t useP = floating[1];
  Bool_t floatDirection = floating[2];
  Bool_t floatPosition = floating[3];
  Bool_t floatIntensity = floating[4];
  Bool_t floatAlpha = floating[5];
  Bool_t useFlat = floating[6];
  Bool_t useLeslie = floating[7];
  Bool_t floatConeAngle = floating[9];
  Bool_t floatLightSpeed = floating[10];

  new_t = t + gRandom->Gaus(0,1.0);
  if (floatPosition) {
    new_r_x = r_x + gRandom->Gaus(0,20);
    new_r_y = r_y + gRandom->Gaus(0,20);
    new_r_z = r_z + gRandom->Gaus(0,20);
  }
  else {
//    new_t = t;
    new_r_x = r_x;
    new_r_y = r_y;
    new_r_z = r_z;
  }
  if (floatDirection) {
    for (Int_t cone=0; cone<nCones; cone++) {
      Double_t rotX = gRandom->Gaus(0,0.1);
      Double_t rotY = gRandom->Gaus(0,0.1);
      Double_t rotZ = gRandom->Gaus(0,0.1);
/*      new_u_x[cone] = cos(rotX)*u_x[cone] + sin(rotX)*cos(rotZ)*u_y[cone] + sin(rotX)*sin(rotZ)*u_z[cone];
      new_u_y[cone] = -sin(rotX)*u_x[cone] + cos(rotX)*cos(rotZ)*u_y[cone] + cos(rotX)*sin(rotZ)*u_z[cone];
      new_u_z[cone] = -sin(rotZ)*u_y[cone] + cos(rotZ)*u_z[cone];
*/
// Assumes small angle rotations
      new_u_x[cone] =      u_x[cone] + rotZ*u_y[cone] - rotY*u_z[cone];
      new_u_y[cone] =-rotZ*u_x[cone] +      u_y[cone] + rotX*u_z[cone];
      new_u_z[cone] = rotY*u_x[cone] - rotX*u_y[cone] +      u_z[cone];

      Double_t u = sqrt(pow(new_u_x[cone],2)+pow(new_u_y[cone],2)+pow(new_u_z[cone],2));

      new_u_x[cone] /= u;
      new_u_y[cone] /= u;
      new_u_z[cone] /= u;
    }
  }
  else {
    for (Int_t cone=0; cone<nCones; cone++) {
      new_u_x[cone] = u_x[cone];
      new_u_y[cone] = u_y[cone];
      new_u_z[cone] = u_z[cone];
    }
  }
  if (floatIntensity) {
    Double_t total = 0;
    for (Int_t cone=0; cone<nCones; cone++) {
      new_coneIntensity[cone] = coneIntensity[cone] + gRandom->Gaus(0,2.0);
      total += new_coneIntensity[cone];
    }
    if (useFlat){
      new_flatIntensity = flatIntensity + gRandom->Gaus(0,2.0);
    } 
    else {
      new_flatIntensity = 0;
    }
    if (useLeslie){
      new_leslieIntensity = leslieIntensity + gRandom->Gaus(0,2.0);
    }
    else {
      new_leslieIntensity = 0;
    }
    total += new_flatIntensity;
    total += new_leslieIntensity;
    new_unknownIntensity = unknownIntensity + gRandom->Gaus(0,2.0);
//    new_unknownIntensity = 1-total;
    
  }
  else {
    for (Int_t cone=0; cone<nCones; cone++) {
      new_coneIntensity[cone] = coneIntensity[cone];
    }
    new_flatIntensity = flatIntensity;
    new_leslieIntensity = leslieIntensity;
    new_unknownIntensity = unknownIntensity;
  }
  if (floatAlpha) {
    new_alpha = alpha + gRandom->Gaus(0,0.01);
  }
  else {
    new_alpha = alpha;
  }
  new_trms = trms + gRandom->Gaus(0,0.03);
  if (floatConeAngle){
    new_coneAngle = coneAngle + gRandom->Gaus(0,0.001);
  }
  else{
    new_coneAngle = coneAngle;
  }
  if (floatLightSpeed){
  new_lightSpeed = lightSpeed + gRandom->Gaus(0,0.02);
  //new_lightSpeed = lightSpeed + gRandom->Gaus(0,0.002);
  }
  else{
    new_lightSpeed = lightSpeed;
  }
}


Double_t CalculateLikelihood(Int_t *nPMTs, Int_t *pmtN, Int_t *pmtType, Float_t *pmtT, Float_t *pmtX, Float_t *pmtY, Float_t *pmtZ,
			     Double_t t, Double_t x, Double_t y, Double_t z,
			     Double_t *coneIntensity, Double_t flatIntensity, Double_t leslieIntensity, Double_t unknownIntensity, Double_t *ux, Double_t *uy, Double_t *uz,
			     Double_t alpha, Double_t trms, 
			     Double_t coneAngle, Double_t lightSpeed, Int_t nCones, Bool_t *floating) {
  Bool_t useFlat = floating[6];
  Bool_t useLeslie = floating[7];
  Double_t pmtArea = 1000; //cm^2: This is a ballpark.
  Double_t lnP = 0;
  Double_t Qi_x;
  Double_t Qi_y;
  Double_t Qi_z;
  Double_t Qi;
  Double_t t_i;
  Double_t correctedTime;
  Double_t promptLikelihood;
  Double_t flatLikelihood;
  Double_t leslieLikelihood;
  Double_t unknownLikelihood;
  Double_t coneLikelihood[3];
  Double_t thisCosAngle;
  Double_t areaElementCos;
  Double_t normparam;
  Double_t gausparam;
  Double_t expparam;
  Double_t absexpparam;
  Double_t powparam;
  Double_t powexpparam;
  Double_t cosDiff;

  Double_t r = sqrt(x*x + y*y + z*z);
  Double_t lesliex, lesliey, lesliez;
  if (x*x+y*y+z*z>0){
    lesliex = x * 600.0/r;
    lesliey = y * 600.0/r;
    lesliez = z * 600.0/r;
  }
  else{
    lesliex=600;
    lesliey=0;
    lesliez=0;
  }
  Double_t thisPMTX;
  Double_t thisPMTY;
  Double_t thisPMTZ;
  Double_t thisPMTR;
  Double_t Qli_x;
  Double_t Qli_y;
  Double_t Qli_z;
  Double_t Qli;
  Double_t leslieCorrectedTime;
  Double_t leslieAreaElementCos;
  Double_t leslieCosAngle;
  Double_t leslieTimeLikelihood;
  Double_t eToTheMinusLambda[10];
  Double_t productOfeToTheMinusLambdas;
  Double_t likelihood;
  Double_t pdfProb[10];
  Double_t timeLikelihoods[10];
  Double_t pdfDenom;

  Bool_t thisPMThit[9728], useThisPMT[9728];
  Double_t spacedTimes[9728];

  for (Int_t i=0; i<9728; i++){
    spacedTimes[i]=0.0;
    thisPMThit[i]=0;
    if (pmtType[i]==1){
      useThisPMT[i]=1;
    }
    else{
      useThisPMT[i]=0;
    }
  }
  for (Int_t i=0; i<*nPMTs; i++) {
    spacedTimes[pmtN[i]]=pmtT[i];
    thisPMThit[pmtN[i]]=1;
  }
  for (Int_t i=0; i<9728; i++) {
    if (useThisPMT[i]){  //Check that this is a "normal" PMT
//    if (useThisPMT[i] && thisPMThit[i]){  //Old style
      thisPMTX=pmtX[i];
      thisPMTY=pmtY[i];
      thisPMTZ=pmtZ[i];
      thisPMTR=sqrt(thisPMTX*thisPMTX + thisPMTY*thisPMTY + thisPMTZ*thisPMTZ);
      Qi_x = thisPMTX - x;
      Qi_y = thisPMTY - y;
      Qi_z = thisPMTZ - z;
      Qi = sqrt(Qi_x*Qi_x+Qi_y*Qi_y+Qi_z*Qi_z);
      areaElementCos = (Qi_x*thisPMTX + Qi_y*thisPMTY + Qi_z*thisPMTZ) / (thisPMTR * Qi );
      Qli_x = thisPMTX - lesliex;
      Qli_y = thisPMTY - lesliey;
      Qli_z = thisPMTZ - lesliez;
      Qli = sqrt(Qli_x*Qli_x+Qli_y*Qli_y+Qli_z*Qli_z);
      leslieAreaElementCos = (Qli_x*thisPMTX + Qli_y*thisPMTY + Qli_z*thisPMTZ) / (thisPMTR * Qli );

      for (Int_t cone=0; cone<nCones; cone++) {
        coneLikelihood[cone]=coneIntensity[cone];
//        coneLikelihood[cone]=1.0;
      }
      unknownLikelihood = unknownIntensity;
      flatLikelihood = flatIntensity;
      leslieLikelihood = leslieIntensity;
//      unknownLikelihood = 1.0;
//      flatLikelihood = 1.0;
//      leslieLikelihood = 1.0; 

      for (Int_t cone=0; cone<nCones; cone++) {
        thisCosAngle = (Qi_x*ux[cone] + Qi_y*uy[cone] + Qi_z*uz[cone])/Qi;
        cosDiff=thisCosAngle-cos(coneAngle);

        normparam = (13.24*alpha*alpha - 27.52*alpha +18.49)*(370*cos(coneAngle)+171)/441;
        gausparam = 2.83*alpha -3.21;
        expparam = (-2.6*alpha*alpha + 4.2*alpha -1)*(-1.05*cos(coneAngle) +1.46)/.694;
        absexpparam = -10.8*alpha*alpha + 17.2*alpha -5.5;
        powparam = 3.07*alpha*alpha - 3.52*alpha -2.21;
        powexpparam = .37;
        
        coneLikelihood[cone] = normparam * exp(
          gausparam * cosDiff * cosDiff
          + expparam * cosDiff
          + absexpparam * abs(thisCosAngle)
          + powparam * pow (abs(cosDiff), powexpparam)
        );
        coneLikelihood[cone] /= (2*3.1415927);
        coneLikelihood[cone] *= pmtArea/(Qi*Qi);

        if (areaElementCos>.574){
          coneLikelihood[cone] *= areaElementCos;
        }
        else {
          coneLikelihood[cone] *= 0;
        }
      }
      leslieCosAngle = (Qli_x*lesliex + Qli_y*lesliey + Qli_z*lesliez)/(Qli*600);
      leslieLikelihood *= 1.33579*sqrt(abs(leslieCosAngle+.07));
     
      //solid angle of a PMT is pmtArea * cosAreaElement / r^2
      flatLikelihood /= (4.0*3.1415927);
      unknownLikelihood /= (4.0*3.1415927); 
      leslieLikelihood /= (2*3.1415927); // pdf is in cos(theta) not theta, so 2pi instead of 4pi

      flatLikelihood *= pmtArea/(Qi*Qi); 
      leslieLikelihood *= pmtArea/(Qli*Qli);
      unknownLikelihood *= pmtArea/(thisPMTR*thisPMTR);

      if (areaElementCos>.574){ //because PMT response is kind of a theta function
        flatLikelihood *= areaElementCos;
      }
      else{
        flatLikelihood *= 0;
      }
      if (leslieAreaElementCos>.574){
        leslieLikelihood *= leslieAreaElementCos;
      }
      else{
        leslieLikelihood *= 0;
      }
      // unknown is always from center, so cosUnknownAreaElement===1


      //In this order: [cone[0], ... , cones[nCones-1], unknown, flat, leslie]
      pdfDenom=0;
      for (Int_t cone=0; cone<nCones; cone++){
        eToTheMinusLambda[cone] = exp(-coneLikelihood[cone]);
        pdfProb[cone] = coneIntensity[cone];
        pdfDenom += coneIntensity[cone];
      }
      eToTheMinusLambda[nCones] = exp(-unknownLikelihood);
      pdfProb[nCones] = unknownIntensity;
      pdfDenom += unknownIntensity;
      if (useFlat){
        eToTheMinusLambda[nCones+useFlat] = exp(-flatLikelihood);
        pdfProb[nCones+useFlat] = flatIntensity;
        pdfDenom += flatIntensity;
      }
      if (useLeslie){
        eToTheMinusLambda[nCones+useFlat+useLeslie] = exp(-leslieLikelihood);
        pdfProb[nCones+useFlat+useLeslie] = leslieIntensity;
        pdfDenom += leslieIntensity;
      }
      productOfeToTheMinusLambdas = 1.0;
      for (Int_t pdf=0; pdf<nCones+useFlat+useLeslie; pdf++){
        productOfeToTheMinusLambdas *= eToTheMinusLambda[pdf];
        pdfProb[pdf] /= pdfDenom;
      }
        
      if (thisPMThit[i]){

        t_i = spacedTimes[i];
        correctedTime = t_i - t - Qi/lightSpeed;
        leslieCorrectedTime = t_i - t - Qli/lightSpeed;
        promptLikelihood = exp(-pow(correctedTime/trms,2)/2)/(sqrt(2.*3.14159)*trms);
        leslieTimeLikelihood = .0965491*exp(-(leslieCorrectedTime*leslieCorrectedTime)/(2*2.42725*2.42725)*(leslieCorrectedTime<2.95655)+(leslieCorrectedTime>2.95655)*(-0.143011*leslieCorrectedTime+-0.345818)*(leslieCorrectedTime<14.4765)+(leslieCorrectedTime>14.4765)*(-0.0457376*leslieCorrectedTime-1.72445));

        for (Int_t cone=0; cone<nCones; cone++){
          timeLikelihoods[cone]=promptLikelihood; 
        }
        timeLikelihoods[nCones] = 1.0/400.0;
        if (useFlat){
          timeLikelihoods[nCones+useFlat]=promptLikelihood;
        }
        if (useLeslie){
          timeLikelihoods[nCones+useFlat+useLeslie]=leslieTimeLikelihood;
        }
        
        likelihood=0;
        for (Int_t pdf=0; pdf<nCones+useFlat+useLeslie; pdf++){
//          likelihood += timeLikelihoods[pdf] * (1-eToTheMinusLambda[pdf]) * productOfeToTheMinusLambdas / eToTheMinusLambda[pdf] * pdfProb[pdf];
          likelihood += timeLikelihoods[pdf] * (1-eToTheMinusLambda[pdf]) * pdfProb[pdf];
        }

        lnP += log(likelihood);
      }
      else { //if PMT not hit
        lnP += log(productOfeToTheMinusLambdas);
      }
    }
  }

  //cout << endl;
  //gSystem->Sleep(1000);
  lnP += -pow((coneAngle-.7446)/.017,2.0)/2.0;
  lnP += -pow((trms-1.707)/0.17,2.0)/2.0;
  lnP += -pow((alpha-.751)/0.01,2.0)/2.0;
  lnP += -pow((lightSpeed-21.5)/0.3,2.0)/2.0;
  //lnP += -pow((alpha-0.781)/0.046,2.0)/2.0;
  //lnP += -pow((beta-0.079)/0.008,2.0)/2.0;
  //lnP += -pow((gamma-0.0299)/0.0010,2.0)/2.0;
  //lnP += log(x*x+y*y+z*z);
  
  return lnP;
}

void BayesGauss(TFile *inFile, Int_t start_event, Int_t n_events, Bool_t *floating) {
  TTree *tPMT = (TTree*)inFile->Get("tPMT");
//  Double_t *fit_x = (Double_t*) tPMT->GetLeaf("fit_x")->GetValuePointer();
//  Double_t *fit_y = (Double_t*) tPMT->GetLeaf("fit_y")->GetValuePointer();
//  Double_t *fit_z = (Double_t*) tPMT->GetLeaf("fit_z")->GetValuePointer();
//  Double_t *fit_ux = (Double_t*) tPMT->GetLeaf("fit_ux")->GetValuePointer();
//  Double_t *fit_uy = (Double_t*) tPMT->GetLeaf("fit_uy")->GetValuePointer();
//  Double_t *fit_uz = (Double_t*) tPMT->GetLeaf("fit_uz")->GetValuePointer();

  Int_t *nPMTs = (Int_t*) tPMT->GetLeaf("nPMTs")->GetValuePointer();
  Int_t *pmtN = (Int_t*) tPMT->GetLeaf("pmtN")->GetValuePointer();
  Double_t *pmtX = getPMTXs();
  Double_t *pmtY = getPMTYs();
  Double_t *pmtZ = getPMTZs();
  Int_t *pmtType = getPMTTypes();
//  Float_t *pmtX = (Float_t*) tPMT->GetLeaf("pmtX")->GetValuePointer();
//  Float_t *pmtY = (Float_t*) tPMT->GetLeaf("pmtY")->GetValuePointer();
//  Float_t *pmtZ = (Float_t*) tPMT->GetLeaf("pmtZ")->GetValuePointer();
  Float_t *pmtT = (Float_t*) tPMT->GetLeaf("pmtT")->GetValuePointer();
  Float_t *trueX = (Float_t*) tPMT->GetLeaf("trueX")->GetValuePointer();
  Float_t *trueY = (Float_t*) tPMT->GetLeaf("trueY")->GetValuePointer();
  Float_t *trueZ = (Float_t*) tPMT->GetLeaf("trueZ")->GetValuePointer();
  Float_t *trueUX = (Float_t*) tPMT->GetLeaf("trueUX")->GetValuePointer();
  Float_t *trueUY = (Float_t*) tPMT->GetLeaf("trueUY")->GetValuePointer();
  Float_t *trueUZ = (Float_t*) tPMT->GetLeaf("trueUZ")->GetValuePointer();

  static const Int_t steps = <NUMBEROFSTEPS>;
  static const Int_t maxCones = <NUMBEROFCONES>;

  Int_t nCones = (Int_t) maxCones;
  Int_t maxSteps = (Int_t) steps;
  Int_t step[steps];
  Double_t t[steps];
  Double_t r_x[steps];
  Double_t r_y[steps];
  Double_t r_z[steps];
  Double_t r[steps];
  Double_t coneIntensity[steps][maxCones];
  Double_t u_x[steps][maxCones];
  Double_t u_y[steps][maxCones];
  Double_t u_z[steps][maxCones];
  Double_t flatIntensity[steps];
  Double_t leslieIntensity[steps];
  Double_t unknownIntensity[steps];

  Double_t coneIntensity_0[steps];
  Double_t u_x_0[steps];
  Double_t u_y_0[steps];
  Double_t u_z_0[steps];
  Double_t coneIntensity_1[steps];
  Double_t u_x_1[steps], u_y_1[steps], u_z_1[steps];
  Double_t coneIntensity_2[steps], u_x_2[steps], u_y_2[steps], u_z_2[steps];
  Double_t coneIntensity_3[steps], u_x_3[steps], u_y_3[steps], u_z_3[steps];
  Double_t coneIntensity_4[steps], u_x_4[steps], u_y_4[steps], u_z_4[steps];

  Double_t alpha[steps];
  Double_t trms[steps];
  Double_t coneAngle[steps];
  Double_t lightSpeed[steps];  
  Double_t lnP[steps];
  Double_t fr_x=0, fr_y=0, fr_z=0, fu_angle[maxCones], fr=0, mt=0, mr_x=0, mr_y=0, mr_z=0, mr=0, mconeIntensity[maxCones], mflatIntensity=0, mleslieIntensity=0, munknownIntensity=0, mu_x[maxCones], mu_y[maxCones], mu_z[maxCones], malpha=0,  mtrms=0, mconeAngle=0, mlightSpeed=0;
  Double_t lr=0, lr_x=0, lr_y=0, lr_z=0, lu_x[maxCones], lu_y[maxCones], lu_z[maxCones], lconeIntensity[maxCones];
  Double_t sr_x=0, sr_y=0, sr_z=0, sr=0;

  char filename[255] = "data/";
  if (floating[0] && floating[1]){
    strcat(filename,"PT");
  }else if (floating[0]){
    strcat(filename,"T_");
  }else if (floating[1]){
    strcat(filename,"P_");
  }else{
    strcat(filename,"__");
  }
  if (floating[2]){
    strcat(filename,"D");
  }else{
    strcat(filename,"_");
  }
  if (floating[3]){
    strcat(filename,"X");
  }else{
    strcat(filename,"_");
  }
  if (floating[4]){
    strcat(filename,"I");
  }else{
    strcat(filename,"_");
  }
  if (floating[5]){
    strcat(filename,"A");
  }else{
    strcat(filename,"_");
  }
  if (floating[6]){
    strcat(filename,"F");
  }else{
    strcat(filename,"_");
  }
  if (floating[7]){
    strcat(filename,"V");
  }else{
    strcat(filename,"_");
  }
  if (floating[8]){
    strcat(filename,"a");
  }else{
    strcat(filename,"_");
  }
  if (floating[9]){
    strcat(filename,"C");
  }else{
    strcat(filename,"_");
  }
  if (floating[10]){
    strcat(filename,"L");
  }else{
    strcat(filename,"_");
  }
  strcat(filename,"_");
  char stepschar[20];
  sprintf(stepschar, "%d", steps);
  strcat(filename, stepschar);
  strcat(filename, "steps_");
  char coneschar[20];
  sprintf(coneschar, "%d", maxCones);
  strcat(filename, coneschar);
  strcat(filename, "cones.root");

  TFile *bayesFile = new TFile(filename,"RECREATE");
  TTree *bayesTree = new TTree("bayesTree","bayesTree");
  //Note the positions, etc. are floats
  bayesTree->Branch("nCones",&nCones,"nCones/I");
  bayesTree->Branch("steps",&maxSteps,"steps/I");
  bayesTree->Branch("step",step,"step[steps]/I");
  bayesTree->Branch("t",t,"t[steps]/D");
  bayesTree->Branch("r_x",r_x,"r_x[steps]/D");
  bayesTree->Branch("r_y",r_y,"r_y[steps]/D");
  bayesTree->Branch("r_z",r_z,"r_z[steps]/D");
  bayesTree->Branch("r",r,"r[steps]/D");
//  bayesTree->Branch("coneIntensity",coneIntensity,"coneIntensity[80000][3]/D");
  bayesTree->Branch("flatIntensity",flatIntensity,"flatIntensity[steps]/D");
  bayesTree->Branch("leslieIntensity",leslieIntensity,"leslieIntensity[steps]/D");
  bayesTree->Branch("unknownIntensity",unknownIntensity,"unknownIntensity[steps]/D");
//  bayesTree->Branch("u_x",u_x,"u_x[80000][3]/D");
//  bayesTree->Branch("u_y",u_y,"u_y[80000][3]/D");
//  bayesTree->Branch("u_z",u_z,"u_z[80000][3]/D");

  if (maxCones>0){
    bayesTree->Branch("coneIntensity_0",coneIntensity_0,"coneIntensity_0[steps]/D");
    bayesTree->Branch("u_x_0",u_x_0,"u_x_0[steps]/D");
    bayesTree->Branch("u_y_0",u_y_0,"u_y_0[steps]/D");
    bayesTree->Branch("u_z_0",u_z_0,"u_z_0[steps]/D");
  }
  if (maxCones>1){
    bayesTree->Branch("coneIntensity_1",coneIntensity_1,"coneIntensity_1[steps]/D");
    bayesTree->Branch("u_x_1",u_x_1,"u_x_1[steps]/D");
    bayesTree->Branch("u_y_1",u_y_1,"u_y_1[steps]/D");
    bayesTree->Branch("u_z_1",u_z_1,"u_z_1[steps]/D");
  }
  if (maxCones>2){
    bayesTree->Branch("coneIntensity_2",coneIntensity_2,"coneIntensity_2[steps]/D");
    bayesTree->Branch("u_x_2",u_x_2,"u_x_2[steps]/D");
    bayesTree->Branch("u_y_2",u_y_2,"u_y_2[steps]/D");
    bayesTree->Branch("u_z_2",u_z_2,"u_z_2[steps]/D");
  }
  if (maxCones>3){
    bayesTree->Branch("coneIntensity_3",coneIntensity_3,"coneIntensity_3[steps]/D");
    bayesTree->Branch("u_x_3",u_x_3,"u_x_3[steps]/D");
    bayesTree->Branch("u_y_3",u_y_3,"u_y_3[steps]/D");
    bayesTree->Branch("u_z_3",u_z_3,"u_z_3[steps]/D");
  }
  if (maxCones>4){
    bayesTree->Branch("coneIntensity_4",coneIntensity_4,"coneIntensity_4[steps]/D");
    bayesTree->Branch("u_x_4",u_x_4,"u_x_4[steps]/D");
    bayesTree->Branch("u_y_4",u_y_4,"u_y_4[steps]/D");
    bayesTree->Branch("u_z_4",u_z_4,"u_z_4[steps]/D");
  }

  bayesTree->Branch("alpha",alpha,"alpha[steps]/D");
  bayesTree->Branch("trms",trms,"trms[steps]/D");
  bayesTree->Branch("coneAngle",coneAngle,"coneAngle[steps]/D");
  bayesTree->Branch("lightSpeed",lightSpeed,"lightSpeed[steps]/D");
  bayesTree->Branch("lnP",lnP,"lnP[steps]/D");
  bayesTree->Branch("fr_x",&fr_x,"fr_x/D");
  bayesTree->Branch("fr_y",&fr_y,"fr_y/D");
  bayesTree->Branch("fr_z",&fr_z,"fr_z/D");
  bayesTree->Branch("fr",&fr,"fr/D");
  bayesTree->Branch("fu_angle",fu_angle,"fu_angle[nCones]/D");
  bayesTree->Branch("mt",&mt,"mt/D");
  bayesTree->Branch("mr_x",&mr_x,"mr_x/D");
  bayesTree->Branch("mr_y",&mr_y,"mr_y/D");
  bayesTree->Branch("mr_z",&mr_z,"mr_z/D");
  bayesTree->Branch("mr",&mr,"mr/D");
  bayesTree->Branch("lr_x",&lr_x,"lr_x/D");
  bayesTree->Branch("lr_y",&lr_y,"lr_y/D");
  bayesTree->Branch("lr_z",&lr_z,"lz_x/D");
  bayesTree->Branch("lr",&lr,"lr/D");
  bayesTree->Branch("mconeIntensity",mconeIntensity,"mconeIntensity[nCones]/D");
  bayesTree->Branch("mflatIntensity",&mflatIntensity,"mflatIntensity/D");
  bayesTree->Branch("mleslieIntensity",&mleslieIntensity,"mleslieIntensity/D");
  bayesTree->Branch("munknownIntensity",&munknownIntensity,"munknownIntensity/D");
  bayesTree->Branch("lconeIntensity",lconeIntensity,"lconeIntensity[nCones]/D");
  bayesTree->Branch("mu_x",mu_x,"mu_x[nCones]/D");
  bayesTree->Branch("mu_y",mu_y,"mu_y[nCones]/D");
  bayesTree->Branch("mu_z",mu_z,"mu_z[nCones]/D");
  bayesTree->Branch("lu_x",lu_x,"lu_x[nCones]/D");
  bayesTree->Branch("lu_y",lu_y,"lu_y[nCones]/D");
  bayesTree->Branch("lu_z",lu_z,"lu_z[nCones]/D");
  bayesTree->Branch("malpha",&malpha,"malpha/D");
  bayesTree->Branch("mtrms",&mtrms,"mtrms/D");
  bayesTree->Branch("mconeAngle",&mconeAngle,"mconeAngle/D");
  bayesTree->Branch("mlightSpeed",&mlightSpeed,"mlightSpeed/D");
  bayesTree->Branch("sr_x",&sr_x,"sr_x/D");
  bayesTree->Branch("sr_y",&sr_y,"sr_y/D");
  bayesTree->Branch("sr_z",&sr_z,"sr_z/D");
  bayesTree->Branch("sr",&sr,"sr/D");

  gRandom->SetSeed();

  Long64_t nentries = tPMT->GetEntries();
  cout << nentries << " entries total" << endl;

  //Cycles over events
  for (Int_t k=0; k<start_event; k++) {
    tPMT->GetEntry(k);
    bayesTree->Fill();
  }
//  Double_t u_length=0.;
  for (Int_t k=start_event; k<start_event+n_events; k++) {
    tPMT->GetEntry(k);
    cout << k << endl;
//    cout << *nPMTs << endl;
//    nCones = 3;
    step[0] = 0; 
    r_x[0] = 0;
    r_y[0] = 0;
    r_z[0] = 0;
//    r_x[0] = *trueX;
//    r_y[0] = *trueY;
//    r_z[0] = *trueZ;

    r[0] = sqrt(pow(r_x[0],2)+pow(r_y[0],2)+pow(r_z[0],2));
    t[0] = GetT0(nPMTs, pmtN, pmtType, pmtT, pmtX, pmtY, pmtZ, r_x[0], r_y[0], r_z[0]);

    coneIntensity[0][0] = 80;
    coneIntensity[0][1] = 5;
    flatIntensity[0]=5;
    leslieIntensity[0]=0.;
    unknownIntensity[0]=5;
//    u_x[0][0] = *trueUX;
//    u_y[0][0] = *trueUY;
//    u_z[0][0] = *trueUZ;
    
    u_x[0][0] = 1;
    u_y[0][0] = 0;
    u_z[0][0] = 0;
/*
    for (Int_t i=0; i<*nPMTs; i++) {
      u_x[0][0] += pmtX[i];
      u_y[0][0] += pmtY[i];
      u_z[0][0] += pmtZ[i];
    }
    u_length=sqrt(u_x[0][0]*u_x[0][0]+u_y[0][0]*u_y[0][0]+u_z[0][0]*u_z[0][0]);
    u_x[0][0] /=u_length;
    u_y[0][0] /=u_length;
    u_z[0][0] /=u_length;
*/
    for (Int_t cone=1; cone<nCones; cone++) {
//      coneIntensity[0][cone] = (1-coneIntensity[0][0]-flatIntensity[0]-leslieIntensity[0]-unknownIntensity[0])/(Double_t)(nCones-1);
      u_x[0][cone] = -(u_x[0][0]);
      u_y[0][cone] = -(u_y[0][0]);
      u_z[0][cone] = -(u_z[0][0]);
    }
    alpha[0] = .751;
    trms[0] = 1.71;
    coneAngle[0] = 0.7446;
    lightSpeed[0] = 21.5;
    //lightSpeed[0] = 0.10;
    lnP[0] = CalculateLikelihood(nPMTs, pmtN, pmtType, pmtT, pmtX, pmtY, pmtZ, t[0], r_x[0], r_y[0], r_z[0], coneIntensity[0], flatIntensity[0], leslieIntensity[0], unknownIntensity[0], u_x[0], u_y[0], u_z[0], alpha[0], trms[0], coneAngle[0], lightSpeed[0], nCones, floating);
//    lnP[0] = CalculateLikelihood(nUsedPMTs, pmtTuse, pmtXuse, pmtYuse, pmtZuse, t[0], r_x[0], r_y[0], r_z[0], coneIntensity[0], u_x[0], u_y[0], u_z[0], alpha[0], beta[0], gamma[0], trms[0], angle[0], coneAngle[0], lightSpeed[0], nCones);

    //Guesses fits
    for (Int_t i=0; i<steps-1; i++) {
      if (i%1000 == 0){
        cout << "step " << i << endl;
      }
      step[i+1] = i+1;
      GetNewStep(t[i], r_x[i], r_y[i], r_z[i], coneIntensity[i], flatIntensity[i], leslieIntensity[i], unknownIntensity[i], u_x[i], u_y[i], u_z[i], alpha[i], trms[i], coneAngle[i], lightSpeed[i], t[i+1], r_x[i+1], r_y[i+1], r_z[i+1], coneIntensity[i+1], flatIntensity[i+1], leslieIntensity[i+1], unknownIntensity[i+1], u_x[i+1], u_y[i+1], u_z[i+1], alpha[i+1], trms[i+1], coneAngle[i+1], lightSpeed[i+1], nCones, floating);
      r[i+1] = sqrt(pow(r_x[i+1],2)+pow(r_y[i+1],2)+pow(r_z[i+1],2));
      Bool_t isGood = (alpha[i+1]<1 && alpha[i+1]>0 && /* (pow(r_x[i+1],2)+pow(r_y[i+1],2)+pow(r_z[i+1],2))<600*600 &&*/ trms[i+1]>0 && coneAngle[i+1]>0 && coneAngle[i+1]<1 && lightSpeed[i+1]>0 && lightSpeed[i+1]<30);
      for (Int_t cone=0; cone<nCones; cone++) {
	isGood = isGood && (coneIntensity[i+1][cone]>0/* && coneIntensity[i+1][cone]<1*/);
      }
      isGood = isGood && (flatIntensity[i+1]>-1e-15/* && flatIntensity[i+1]<1*/);
      isGood = isGood && (leslieIntensity[i+1]>-1e-15/* && leslieIntensity[i+1]<1*/);
      isGood = isGood && (unknownIntensity[i+1]>-1e-15/* && unknownIntensity[i+1]<1*/);
      if (isGood) {
	lnP[i+1] = CalculateLikelihood(nPMTs, pmtN, pmtType, pmtT, pmtX, pmtY, pmtZ, t[i+1], r_x[i+1], r_y[i+1], r_z[i+1], coneIntensity[i+1], flatIntensity[i+1], leslieIntensity[i+1], unknownIntensity[i+1], u_x[i+1], u_y[i+1], u_z[i+1], alpha[i+1], trms[i+1], coneAngle[i+1], lightSpeed[i+1], nCones, floating);
//	lnP[i+1] = CalculateLikelihood(nUsedPMTs, pmtTuse, pmtXuse, pmtYuse, pmtZuse, t[i+1], r_x[i+1], r_y[i+1], r_z[i+1], coneIntensity[i+1], u_x[i+1], u_y[i+1], u_z[i+1], alpha[i+1], beta[i+1], gamma[i+1], trms[i+1], angle[i+1], coneAngle[i+1], lightSpeed[i+1], nCones);
      }

      //Go to more likely position, probability of going to less likely
      Double_t logProb = TMath::Log(gRandom->Uniform(0,1));
      //cout << i << "\t" << lnP[i+1] << "\t" << lnP[i] << "\t" << lnP[i+1]-lnP[i] << "\t" << logProb << "\t";
      //cout << u_x[i][0] << "\t" << u_y[i][0] << "\t" << u_z[i][0] << "\t" << u_x[i+1][0] << "\t" << u_y[i+1][0] << "\t" << u_z[i+1][0] << "\t";
      //cout << u_x[i][1] << "\t" << u_y[i][1] << "\t" << u_z[i][1] << "\t" << u_x[i+1][1] << "\t" << u_y[i+1][1] << "\t" << u_z[i+1][1] << "\t";
      //cout << endl;
      if ((lnP[i+1]-lnP[i])<logProb || r[i+1]>824 || !isGood) {
	t[i+1] = t[i];
	r_x[i+1] = r_x[i];
	r_y[i+1] = r_y[i];
	r_z[i+1] = r_z[i];
	r[i+1] = r[i];
	for (Int_t cone=0; cone<nCones; cone++) {
	  coneIntensity[i+1][cone] = coneIntensity[i][cone];
	  u_x[i+1][cone] = u_x[i][cone];
	  u_y[i+1][cone] = u_y[i][cone];
	  u_z[i+1][cone] = u_z[i][cone];
	}
	flatIntensity[i+1] = flatIntensity[i];
	leslieIntensity[i+1] = leslieIntensity[i];
	unknownIntensity[i+1] = unknownIntensity[i];
	alpha[i+1] = alpha[i];
	trms[i+1] = trms[i];
	coneAngle[i+1] = coneAngle[i];
	lightSpeed[i+1] = lightSpeed[i];
	lnP[i+1] = lnP[i];
      }
      if (nCones>0){
        coneIntensity_0[i] = coneIntensity[i][0];
        u_x_0[i] = u_x[i][0];
        u_y_0[i] = u_y[i][0];
        u_z_0[i] = u_z[i][0];
      }
      if (nCones>1){
        coneIntensity_1[i] = coneIntensity[i][1];
        u_x_1[i] = u_x[i][1];
        u_y_1[i] = u_y[i][1];
        u_z_1[i] = u_z[i][1];
      }
      if (nCones>2){
        coneIntensity_2[i] = coneIntensity[i][2];
        u_x_2[i] = u_x[i][2];
        u_y_2[i] = u_y[i][2];
        u_z_2[i] = u_z[i][2];
      }
      if (nCones>3){
        coneIntensity_3[i] = coneIntensity[i][3];
        u_x_3[i] = u_x[i][3];
        u_y_3[i] = u_y[i][3];
        u_z_3[i] = u_z[i][3];
      }
      if (nCones>4){
        coneIntensity_4[i] = coneIntensity[i][4];
        u_x_4[i] = u_x[i][4];
        u_y_4[i] = u_y[i][4];
        u_z_4[i] = u_z[i][4];
      }
    }
    mt = 0;
    mr_x = 0;
    mr_y = 0;
    mr_z = 0;
    mr = 0;
    lr_x = 0;
    lr_y = 0;
    lr_z = 0;
    lr = 0;
    sr_x = 0;
    sr_y = 0;
    sr_z = 0;
    sr = 0;
    for (Int_t cone=0; cone<nCones; cone++) {
      mconeIntensity[cone] = 0;
      mu_x[cone] = 0;
      mu_y[cone] = 0;
      mu_z[cone] = 0;
      lconeIntensity[cone] = 0;
      lu_x[cone] = 0;
      lu_y[cone] = 0;
      lu_z[cone] = 0;
    }
    mflatIntensity = 0;
    mleslieIntensity = 0;
    munknownIntensity = 0;
    malpha = 0;
    mtrms = 0;
    mconeAngle = 0;
    mlightSpeed = 0;
    for (Int_t i = (Int_t)(0.5*steps); i<steps; i++) {
      mt += t[i];
      mr_x += r_x[i];
      mr_y += r_y[i];
      mr_z += r_z[i];
      mr += r[i];
      sr_x += r_x[i]*r_x[i];
      sr_y += r_y[i]*r_y[i];
      sr_z += r_z[i]*r_z[i];
      sr += r[i]*r[i];
      for (Int_t cone=0; cone<nCones; cone++) {
	mconeIntensity[cone] += coneIntensity[i][cone];
	mu_x[cone] += u_x[i][cone];
	mu_y[cone] += u_y[i][cone];
	mu_z[cone] += u_z[i][cone];
      }
      mflatIntensity += flatIntensity[i];
      mleslieIntensity += leslieIntensity[i];
      munknownIntensity += unknownIntensity[i];
      malpha += alpha[i];
      mtrms += trms[i];
      mconeAngle += coneAngle[i];
      mlightSpeed += lightSpeed[i];
    }
    mt /= (Double_t)(0.5*steps);
    mr_x /= (Double_t)(0.5*steps);
    mr_y /= (Double_t)(0.5*steps); 
    mr_z /= (Double_t)(0.5*steps);
    mr /= (Double_t)(0.5*steps);
    sr_x = sr_x/(Double_t)(0.5*steps)-mr_x*mr_x;
    sr_y = sr_y/(Double_t)(0.5*steps)-mr_y*mr_y; 
    sr_z = sr_z/(Double_t)(0.5*steps)-mr_z*mr_z;
    sr = sr/(Double_t)(0.5*steps)-mr*mr;
    for (Int_t cone=0; cone<nCones; cone++) {
      mconeIntensity[cone] /= (Double_t)(0.5*steps);
      mu_x[cone] /= (Double_t)(0.5*steps);
      mu_y[cone] /= (Double_t)(0.5*steps);
      mu_z[cone] /= (Double_t)(0.5*steps);
    }
    mflatIntensity /= (Double_t)(0.5*steps);
    mleslieIntensity /= (Double_t)(0.5*steps);
    munknownIntensity /= (Double_t)(0.5*steps);
    malpha /= (Double_t)(0.5*steps);
    mtrms /= (Double_t)(0.5*steps);
    mconeAngle /= (Double_t)(0.5*steps);
    mlightSpeed /= (Double_t)(0.5*steps);

    for (Int_t i = (Int_t)(steps-1000); i<steps; i++) {
      lr_x += r_x[i];
      lr_y += r_y[i];
      lr_z += r_z[i];
      lr += r[i];
      for (Int_t cone=0; cone<nCones; cone++) {
	lconeIntensity[cone] += coneIntensity[i][cone];
	lu_x[cone] += u_x[i][cone];
	lu_y[cone] += u_y[i][cone];
	lu_z[cone] += u_z[i][cone];
      }
    }
    lr_x /= (Double_t)(1000);
    lr_y /= (Double_t)(1000); 
    lr_z /= (Double_t)(1000);
    lr /= (Double_t)(1000);
    for (Int_t cone=0; cone<nCones; cone++) {
      lconeIntensity[cone] /= (Double_t)(1000);
      lu_x[cone] /= (Double_t)(1000);
      lu_y[cone] /= (Double_t)(1000);
      lu_z[cone] /= (Double_t)(1000);
    }

    fr_x = mr_x - *trueX;
    fr_y = mr_y - *trueY;
    fr_z = mr_z - *trueZ;
    for (Int_t cone=0; cone<nCones; cone++) {
      fu_angle[cone] = (mu_x[cone]*(*trueUX))+(mu_y[cone]*(*trueUY))+(mu_z[cone]*(*trueUZ));
    }
    fr = mr - sqrt(pow(*trueX,2)+pow(*trueY,2)+pow(*trueZ,2));
    bayesTree->Fill();
  }
  bayesTree->Write();
  bayesFile->Close();
}


