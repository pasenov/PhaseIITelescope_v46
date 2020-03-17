//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file electromagnetic/TestEm18/include/EventAction.hh
/// \brief Definition of the EventAction class
//
// $Id: EventAction.hh 82401 2014-06-18 14:43:54Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4EmCalculator.hh"
#include "G4Step.hh"

#include "Randomize.hh"
#include <iomanip>

class DetectorConstruction;
class RunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
  public:
    EventAction(DetectorConstruction*, RunAction*);
   ~EventAction();

  public:
    virtual void BeginOfEventAction(const G4Event* evt);
    virtual void   EndOfEventAction(const G4Event*);

    
    void AddEnergyDeposit1(G4double edep)                     {fEnergyDeposit1  += edep;};
    void AddEnergyDeposit2(G4double edep)                     {fEnergyDeposit2  += edep;};
    void AddPrimaryTrackLength(G4double track)                {fPrimaryTrackLength  += track;};
    void AddSecondary1(G4double ekin)                         {fEnergySecondary1  += ekin;};
    void AddSecondary2(G4double ekin)                         {fEnergySecondary2  += ekin;};
    void AddSecondaryxPolarization(G4double spolarization)    {fSecondaryxPolarization  += spolarization;};
    void AddSecondaryyPolarization(G4double spolarization)    {fSecondaryyPolarization  += spolarization;};
    void AddSecondaryzPolarization(G4double spolarization)    {fSecondaryzPolarization  += spolarization;};
    void AddSecondaryTrackLength(G4double track)              {fSecondaryTrackLength  += track;};
    void AddSecondaryDetTrackLength(G4double track)           {fSecondaryDetTrackLength  += track;};
    void AddTertiary1(G4double ekin)                          {fEnergyTertiary1  += ekin;};
    void AddTertiary2(G4double ekin)                          {fEnergyTertiary2  += ekin;};
    void AddTertiaryxPolarization(G4double spolarization)     {fTertiaryxPolarization  += spolarization;};
    void AddTertiaryyPolarization(G4double spolarization)     {fTertiaryyPolarization  += spolarization;};
    void AddTertiaryzPolarization(G4double spolarization)     {fTertiaryzPolarization  += spolarization;};
    void AddTertiaryTrackLength(G4double track)               {fTertiaryTrackLength  += track;};
    void TrackCheck(G4double trackcheck)                      {fTrack1 = fTrack2; fTrack2 = trackcheck; 
                                                               if (fTrack1 > fTrack2) {
                                                                   G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
                                                                   analysisManager->FillH1(31, fTrack1);
                                                                   //G4cout 
     									//<< "\n Last track length of secondary created in detector calculated = " 
     									//<< G4BestUnit(fTrack1, "Length")
     									//<< G4endl;
                                                               }};
    void AddEnergyStrip1a(G4double edep, G4int i)  {fEnergyStrip1a[i] += edep;};
    void AddEnergyStrip1b(G4double edep, G4int i)  {fEnergyStrip1b[i] += edep;};
    void AddEnergyStrip2a(G4double edep, G4int i)  {fEnergyStrip2a[i] += edep;
									     /*G4cout 
     										<< "\n Energy added to strip, edep = " 
     										<< G4BestUnit(edep, "Energy")
     										<< G4endl;*/ };
    void AddEnergyStrip2b(G4double edep, G4int i)  {fEnergyStrip2b[i] += edep;
									     /*G4cout 
     										<< "\n Energy added to strip, edep = " 
     										<< G4BestUnit(edep, "Energy")
     										<< G4endl; */};

    void AddEnergyPixel(G4double edep, G4int i, G4int j, G4int k, G4int l)  {fEnergyPixel[i][j][k][l] += edep; 
									     /*G4cout 
     										<< "\n Energy added to pixel, edep = " 
     										<< G4BestUnit(edep, "Energy")
     										<< "\n Module number = " 
	     									<< i
     										<< G4endl; */};

    void AddMomentumDirection1(G4double dirx, G4double diry, G4double dirz)  {fMomDir1x += dirx; fMomDir1y += diry; fMomDir1z += dirz;};
    void AddMomentumDirection2(G4double dirx, G4double diry, G4double dirz)  {fMomDir2x += dirx; fMomDir2y += diry; fMomDir2z += dirz;};

    void AddPointDet1entReal(G4ThreeVector vec)  {fPointDet1entReal += vec;};
    void AddPointDet2entReal(G4ThreeVector vec)  {fPointDet2entReal += vec;};

    void AddPointPix1entReal(G4ThreeVector vec)  {fPointPix1entReal += vec;};
    void AddPointPix2entReal(G4ThreeVector vec)  {fPointPix2entReal += vec;};
    void AddPointPix3entReal(G4ThreeVector vec)  {fPointPix3entReal += vec;};
    void AddPointPix4entReal(G4ThreeVector vec)  {fPointPix4entReal += vec;};
    void AddPointPix5entReal(G4ThreeVector vec)  {fPointPix5entReal += vec;};
    void AddPointPix6entReal(G4ThreeVector vec)  {fPointPix6entReal += vec;};
    void AddPointPix7entReal(G4ThreeVector vec)  {fPointPix7entReal += vec;};
    void AddPointPix8entReal(G4ThreeVector vec)  {fPointPix8entReal += vec;};

    void AddA1(G4ThreeVector vec)  {A1 += vec;};
    void AddA(G4ThreeVector vec)   {A += vec;};
    void AddA2(G4ThreeVector vec)  {A2 += vec;};
    void AddB1(G4ThreeVector vec)  {B1 += vec;};
    void AddB(G4ThreeVector vec)   {B += vec;};
    void AddB2(G4ThreeVector vec)  {B2 += vec;};
    void AddC1(G4ThreeVector vec)  {C1 += vec;};
    void AddC(G4ThreeVector vec)   {C += vec;};
    void AddC2(G4ThreeVector vec)  {C2 += vec;};
    void AddD1(G4ThreeVector vec)  {D1 += vec;};
    void AddD(G4ThreeVector vec)   {D += vec;};
    void AddD2(G4ThreeVector vec)  {D2 += vec;};
    void AddE1(G4ThreeVector vec)  {E1 += vec;};
    void AddE(G4ThreeVector vec)   {E += vec;};
    void AddE2(G4ThreeVector vec)  {E2 += vec;};
    void AddF1(G4ThreeVector vec)  {F1 += vec;};
    void AddF(G4ThreeVector vec)   {F += vec;};
    void AddF2(G4ThreeVector vec)  {F2 += vec;};
    void AddG1(G4ThreeVector vec)  {G1 += vec;};
    void AddG(G4ThreeVector vec)   {G += vec;};
    void AddG2(G4ThreeVector vec)  {G2 += vec;};
    void AddH1(G4ThreeVector vec)  {H1 += vec;};
    void AddH(G4ThreeVector vec)   {H += vec;};
    void AddH2(G4ThreeVector vec)  {H2 += vec;};
    void AddI1(G4ThreeVector vec)  {I1 += vec;};
    void AddI(G4ThreeVector vec)   {I += vec;};
    void AddI2(G4ThreeVector vec)  {I2 += vec;};
    void AddJ1(G4ThreeVector vec)  {J1 += vec;};
    void AddJ(G4ThreeVector vec)   {J += vec;};
    void AddJ2(G4ThreeVector vec)  {J2 += vec;};

    void AddAReal(G4ThreeVector vec)   {AReal += vec;};
    void AddBReal(G4ThreeVector vec)   {BReal += vec;};
    void AddCReal(G4ThreeVector vec)   {CReal += vec;};
    void AddDReal(G4ThreeVector vec)   {DReal += vec;};
    void AddEReal(G4ThreeVector vec)   {EReal += vec;};
    void AddFReal(G4ThreeVector vec)   {FReal += vec;};
    void AddGReal(G4ThreeVector vec)   {GReal += vec;};
    void AddHReal(G4ThreeVector vec)   {HReal += vec;};
    void AddIReal(G4ThreeVector vec)   {IReal += vec;};
    void AddJReal(G4ThreeVector vec)   {JReal += vec;};

    void AddEnergyDepositL1(G4double edep)                     {fEnergyDepositL1  += edep;};
    void AddEnergyDepositL2(G4double edep)                     {fEnergyDepositL2  += edep;};
    void AddEnergyDepositL3(G4double edep)                     {fEnergyDepositL3  += edep;};
    void AddEnergyDepositL4(G4double edep)                     {fEnergyDepositL4  += edep;};
    void AddEnergyDepositL5(G4double edep)                     {fEnergyDepositL5  += edep;};
    void AddEnergyDepositL6(G4double edep)                     {fEnergyDepositL6  += edep;};
    void AddEnergyDepositL7(G4double edep)                     {fEnergyDepositL7  += edep;};
    void AddEnergyDepositL8(G4double edep)                     {fEnergyDepositL8  += edep;};
    void AddSecondaryL1(G4double ekin)                         {fEnergySecondaryL1  += ekin;};
    void AddSecondaryL2(G4double ekin)                         {fEnergySecondaryL2  += ekin;};
    void AddSecondaryL3(G4double ekin)                         {fEnergySecondaryL3  += ekin;};
    void AddSecondaryL4(G4double ekin)                         {fEnergySecondaryL4  += ekin;};
    void AddSecondaryL5(G4double ekin)                         {fEnergySecondaryL5  += ekin;};
    void AddSecondaryL6(G4double ekin)                         {fEnergySecondaryL6  += ekin;};
    void AddSecondaryL7(G4double ekin)                         {fEnergySecondaryL7  += ekin;};
    void AddSecondaryL8(G4double ekin)                         {fEnergySecondaryL8  += ekin;};

    void AddEnergyDepositL1Active(G4double edep)               {fEnergyDepositL1Active  += edep;};
        
  private:
    DetectorConstruction*	fDetectorConstruction;
    RunAction*    fRunAction;
    
    G4double      fEnergyDeposit1;
    G4double      fEnergyDeposit2;
    G4double      fPrimaryTrackLength;
    G4double      fEnergySecondary1;
    G4double      fEnergySecondary2;       
    G4double      fSecondaryxPolarization;  
    G4double      fSecondaryyPolarization;  
    G4double      fSecondaryzPolarization; 
    G4double      fSecondaryTrackLength; 
    G4double      fSecondaryDetTrackLength;
    G4double      fEnergyTertiary1;
    G4double      fEnergyTertiary2;   
    G4double      fTertiaryxPolarization;   
    G4double      fTertiaryyPolarization; 
    G4double      fTertiaryzPolarization; 
    G4double      fTertiaryTrackLength;

    G4double      fEnergyDepositL1;
    G4double      fEnergyDepositL2;
    G4double      fEnergyDepositL3;
    G4double      fEnergyDepositL4;
    G4double      fEnergyDepositL5;
    G4double      fEnergyDepositL6;
    G4double      fEnergyDepositL7;
    G4double      fEnergyDepositL8;

    G4double      fEnergyDepositL1Active;

    G4double      fEnergySecondaryL1;
    G4double      fEnergySecondaryL2;
    G4double      fEnergySecondaryL3;
    G4double      fEnergySecondaryL4;
    G4double      fEnergySecondaryL5;
    G4double      fEnergySecondaryL6;
    G4double      fEnergySecondaryL7;
    G4double      fEnergySecondaryL8;

    G4double      fEnergyStrip1a[1017];
    G4double      fEnergyStrip1b[1017];
    G4double      fEnergyStrip2a[1017];
    G4double      fEnergyStrip2b[1017];
    G4double      fWeightStrip1a[1017];
    G4double      fWeightStrip1b[1017];
    G4double      fWeightStrip2a[1017];
    G4double      fWeightStrip2b[1017];
    G4double      fMomDir1x, fMomDir1y, fMomDir1z;
    G4double      fMomDir2x, fMomDir2y, fMomDir2z;
    G4double      fEnergyPixel[17][17][53][81];
    G4double      fWeightPixel[17][17][53][81];

    G4double      enPix;
    G4double      fEnergyPix;

    G4double 	  fStripCenterNRX, fStripCenterNRY, fStripCenterNRZ; //Strip center position if the DUT wasn't rotated
    G4double 	  fStripCenterX, fStripCenterY, fStripCenterZ;

    G4double 	  fPixelCenterNRX, fPixelCenterNRY, fPixelCenterNRZ, fPixelCenterNTZ, fPixelCenterNRX1, fPixelCenterNRY1, fPixelCenterNRZ1; //Pixel center position if the DUT wasn't rotated
    G4double 	  fPixelCenterX, fPixelCenterY, fPixelCenterZ;

    G4double	  poszBPIX12, poszBPIX34, poszBPIX56, poszBPIX78, poszBPIX910, poszBPIX1112, poszBPIX1314, poszBPIX1516;

    G4double	  fDiff1x, fDiff1y, fDiff2x, fDiff2y;

    G4double      Det1SizeZ, Det2SizeZ, Dist, Strip1Depth, Strip1Length, Strip2Depth, Strip2Length, StripDist, StripWidth, StripPitch, posEndArm1, posBeginningArm2, BPIXSizeX, BPIXSizeY, BPIXSizeZ, pixelDepth, pixelX, pixelY, ROChor, ROCvert, XangleDUT, XangleBPIX, YangleBPIX, ElField1, ElField2;

    G4int 	  totalNbStrips, stripNo, totalNbHitStrips;
    G4int         fCharge1, fCharge2, fChargeL1, fCS1a, fCS1b, fCS2a, fCS2b, fChargeStrip1, fChargeStrip2, fChargePix;
    G4int         fHitSensor1, fHitSensor2;
    G4int	  fHitPixelDet[17];

    G4int         fNbHitsStrip1a[1017];
    G4int         fNbHitsStrip1b[1017];
    G4int         fNbHitsStrip2a[1017];
    G4int         fNbHitsStrip2b[1017];
    G4int         fChargeStrip1a[1017];
    G4int         fChargeStrip1b[1017];
    G4int         fChargeStrip2a[1017];
    G4int         fChargeStrip2b[1017];
    G4int         fChargePixel[17][17][53][81];
    G4int         fChargeModule[17];
    G4int         fNbHitsPixel[17][17][53][81];

    G4int         fHitsMultiplicity1b;
    G4int         fHitsMultiplicity2b;
    G4int         fHitsMultiplicityPix[17];
    G4int         fClusterSizeXPixMin, fClusterSizeXPixMax, fClusterSizeXPixROCOfMin, fClusterSizeXPixROCOfMax;
    G4int         fClusterSizeYPixMin, fClusterSizeYPixMax, fClusterSizeYPixROCOfMin, fClusterSizeYPixROCOfMax;
    G4int         fClusterSizeXPix[17];
    G4int         fClusterSizeYPix[17];

    G4int         mod, row, col;
    G4int         fClusterOccupancy[17][161][417];

    G4int         fPixThreshold;

    G4ThreeVector fPointDet1ent, fPointDet1entReal, fPointDet2ent, fPointDet2entReal, fPointPix1ent, fPointPix1entReal, fPointPix2ent, fPointPix2entReal, fPointPix3ent, fPointPix3entReal, fPointPix4ent, fPointPix4entReal, fPointPix5ent, fPointPix5entReal, fPointPix6ent, fPointPix6entReal, fPointPix7ent, fPointPix7entReal, fPointPix8ent, fPointPix8entReal;

    G4ThreeVector A1, A, A2, At, B1, B, B2, Bt, C1, C, C2, Ct, D1, D, D2, Dt, E1, E, E2, Et, F1, F, F2, Ft, G1, G, G2, Gt, H1, H, H2, Ht, I1, I, I2, It, J1, J, J2, Jt, Rav, Af, Bf, Cf, Df, Gf, Hf, If, Jf;
    G4ThreeVector AReal, BReal, CReal, DReal, EReal, FReal, GReal, HReal, IReal, JReal;

    G4ThreeVector A1loc, Aloc, A2loc, Atloc, B1loc, Bloc, B2loc, Btloc, C1loc, Cloc, C2loc, Ctloc, D1loc, Dloc, D2loc, Dtloc, E1loc, Eloc, E2loc, Etloc, F1loc, Floc, F2loc, Ftloc, G1loc, Gloc, G2loc, Gtloc, H1loc, Hloc, H2loc, Htloc, I1loc, Iloc, I2loc, Itloc, J1loc, Jloc, J2loc, Jtloc, Ravloc, Afloc, Bfloc, Cfloc, Dfloc, Gfloc, Hfloc, Ifloc, Jfloc, Af1, Bf1, Cf1, Df1, Gf1, Hf1, If1, Jf1;
    G4ThreeVector ARealloc, BRealloc, CRealloc, DRealloc, ERealloc, FRealloc, GRealloc, HRealloc, IRealloc, JRealloc;

    G4double M[3][3], MA[3][3], MB[3][3], MC[3][3], MD[3][3], MG[3][3], MH[3][3], MI[3][3], MJ[3][3];
    G4double p1, eig1, eig2, eig3, traceM, q, p2, p, detBa, r, phi;
    G4double Iu[3][3], Ba[3][3], Ga[3][3];

    G4double AtAx, AtAy, BtBx, BtBy, CtCx, CtCy, DtDx, DtDy, EtEx, EtEy, FtFx, FtFy, GtGx, GtGy, HtHx, HtHy, ItIx, ItIy, JtJx, JtJy;
    G4double AtAxloc, AtAyloc, BtBxloc, BtByloc, CtCxloc, CtCyloc, DtDxloc, DtDyloc, EtExloc, EtEyloc, FtFxloc, FtFyloc, GtGxloc, GtGyloc, HtHxloc, HtHyloc, ItIxloc, ItIyloc, JtJxloc, JtJyloc;
    G4double tA, tB, tC, tD, tE, tF, tG, tH, tI, tJ;
    G4double Aztr, Bztr, Cztr, Dztr, Gztr, Hztr, Iztr, Jztr;
    G4double Afztr, Bfztr, Cfztr, Dfztr, Gfztr, Hfztr, Ifztr, Jfztr;
    G4double Ga11, Ga00, lambdaVectorX, lambdaVectorY, lambdaVectorZ;

    //primary track's deflection angle in radians
    G4double ftheta;

    G4double lambda;
    G4double lambdaVector[3];

    G4double      fTrack1;
    G4double      fTrack2;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
