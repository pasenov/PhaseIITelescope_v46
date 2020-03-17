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
// $Id: EventAction.cc 82401 2014-06-18 14:43:54Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4EmCalculator.hh"
#include "G4Step.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"
#include <iomanip>
#include <iostream>

#include <stdio.h>
#include <stdlib.h>
//#include <Eigen/Eigenvalues>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(DetectorConstruction* DA, RunAction* RA)
:G4UserEventAction(), fDetectorConstruction(DA), fRunAction(RA)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{
 //geometry parameters
 Det1SizeZ=fDetectorConstruction->GetSize1();
 Det2SizeZ=fDetectorConstruction->GetSize2();
 Dist=fDetectorConstruction->GetDist();
 Strip1Depth=fDetectorConstruction->GetStrip1Depth();
 Strip1Length=fDetectorConstruction->GetStrip1Length();
 Strip2Depth=fDetectorConstruction->GetStrip2Depth();
 Strip2Length=fDetectorConstruction->GetStrip2Length();
 StripDist=fDetectorConstruction->GetStripDist();
 StripWidth=fDetectorConstruction->GetStripWidth();
 StripPitch = StripDist + StripWidth;
 posEndArm1=-(fDetectorConstruction->Getpos_EndArm1Abs());
 posBeginningArm2=fDetectorConstruction->Getpos_BeginningArm2Abs();
 poszBPIX12=fDetectorConstruction->GetZPosBPIX12();
 poszBPIX34=fDetectorConstruction->GetZPosBPIX34();
 poszBPIX56=fDetectorConstruction->GetZPosBPIX56();
 poszBPIX78=fDetectorConstruction->GetZPosBPIX78();
 poszBPIX910=fDetectorConstruction->GetZPosBPIX910();
 poszBPIX1112=fDetectorConstruction->GetZPosBPIX1112();
 poszBPIX1314=fDetectorConstruction->GetZPosBPIX1314();
 poszBPIX1516=fDetectorConstruction->GetZPosBPIX1516();
 BPIXSizeX=fDetectorConstruction->GetSizeXBPIX();
 BPIXSizeY=fDetectorConstruction->GetSizeYBPIX();
 BPIXSizeZ=fDetectorConstruction->GetSizeZBPIX();
 pixelX=fDetectorConstruction->GetPixelPitchX();
 pixelY=fDetectorConstruction->GetPixelPitchY();
 pixelDepth=fDetectorConstruction->GetPixelDepth();
 XangleDUT=fDetectorConstruction->GetDUTangleX();
 XangleBPIX=fDetectorConstruction->GetBPIXangleX();
 YangleBPIX=fDetectorConstruction->GetBPIXangleY();
 ElField1=fDetectorConstruction->GetElField1();
 ElField2=fDetectorConstruction->GetElField2();
 totalNbStrips = 1016;
 ROChor = 80*pixelX + pixelX;
 ROCvert = 52*pixelY + 2*pixelY;
 fPixThreshold = 1700; //1700


 // initialisation per event
 fEnergyDeposit1  =  fEnergyDeposit2  = fEnergySecondary1 = fEnergySecondary2 = fEnergyTertiary1 = fEnergyTertiary2 = 0.;
 fEnergyDepositL1  =  fEnergyDepositL2  =  fEnergyDepositL3  =  fEnergyDepositL4  =  fEnergyDepositL5  =  fEnergyDepositL6  =  fEnergyDepositL7  =  fEnergyDepositL8  = fEnergySecondaryL1 = fEnergySecondaryL2 = fEnergySecondaryL3 = fEnergySecondaryL4 = fEnergySecondaryL5 = fEnergySecondaryL6 = fEnergySecondaryL7 = fEnergySecondaryL8 = fEnergyDepositL1Active = 0.;

 for (G4int j = 0; j <= 1016; j = j + 1) {
    fEnergyStrip1a[j] = 0.;
    fEnergyStrip1b[j] = 0.;
    fEnergyStrip2a[j] = 0.;
    fEnergyStrip2b[j] = 0.;
    fWeightStrip1a[j] = 0.;
    fWeightStrip1b[j] = 0.;
    fWeightStrip2a[j] = 0.;
    fWeightStrip2b[j] = 0.;
    fNbHitsStrip1a[j] = 0;
    fNbHitsStrip1b[j] = 0;
    fNbHitsStrip2a[j] = 0;
    fNbHitsStrip2b[j] = 0;
    fChargeStrip1a[j] = 0;
    fChargeStrip1b[j] = 0;
    fChargeStrip2a[j] = 0;
    fChargeStrip2b[j] = 0;
 }

 totalNbHitStrips = 0;

 fCharge1 = fCharge2 = fChargeL1 = fCS1a = fCS1b = fCS2a = fCS2b = fChargeStrip1 = fChargeStrip2 = fChargePix = 0;

 for (G4int je1 = 0; je1 <= 16; je1 = je1 + 1) {
    fHitsMultiplicityPix[je1] = 0;
    fClusterSizeXPix[je1] = 0;
    fClusterSizeYPix[je1] = 0;
    for (G4int je2 = 0; je2 <= 16; je2 = je2 + 1) {
 	for (G4int je3 = 0; je3 <= 52; je3 = je3 + 1) {
	   for (G4int je4 = 0; je4 <= 80; je4 = je4 + 1) {
		fEnergyPixel[je1][je2][je3][je4] = 0.;
		fChargePixel[je1][je2][je3][je4] = 0;
		fNbHitsPixel[je1][je2][je3][je4] = 0;
 	   }
 	}
    }
 }

 for (G4int je5 = 0; je5 <= 16; je5 = je5 + 1) {
    for (G4int je6 = 0; je6 <= 16; je6 = je6 + 1) {
       for (G4int je7 = 0; je7 <= 52; je7 = je7 + 1) {
          for (G4int je8 = 0; je8 <= 80; je8 = je8 + 1) {
                fWeightPixel[je5][je6][je7][je8] = 0.;
	  }
       }
    }
 }

 fChargeStrip1 = fChargeStrip2 = 0;
 fCS1a = fCS1b = fCS2a = fCS2b = 0;
 for (G4int pix = 0; pix <= 16; pix = pix + 1)  {
    fChargeModule[pix] = 0;
 }

 fPrimaryTrackLength = fSecondaryTrackLength = fSecondaryDetTrackLength = fTertiaryTrackLength = 0.;
 fSecondaryxPolarization = fSecondaryyPolarization = fSecondaryzPolarization = fTertiaryxPolarization = fTertiaryyPolarization = fTertiaryzPolarization = 0.;
 fHitsMultiplicity1b = fHitsMultiplicity2b = 0;
 fHitSensor1 = fHitSensor2 = 0;
 for (G4int pixd = 0; pixd <= 16; pixd = pixd + 1)  {
    fHitPixelDet[pixd] = 0;
 }

 fMomDir1x = fMomDir1y = fMomDir1z =  fMomDir2x = fMomDir2y = fMomDir2z = 0.;

 mod = row = col = -10;

 for (G4int k = 0; k < 17; k = k + 1)  {
    for (G4int i = 0; i < 161; i = i + 1)  {
       for (G4int j = 0; j < 417; j = j + 1)  {
	  fClusterOccupancy[k][i][j] = 0;
       }
    }
 }

 A1 = A = A2 = At = AReal = B1 = B = B2 = Bt = BReal = C1 = C = C2 = Ct = CReal = D1 = D = D2 = Dt = DReal = E1 = E = E2 = Et = EReal = F1 = F = F2 = Ft = FReal = G1 = G = G2 = Gt = GReal = H1 = H = H2 = Ht = HReal = I1 = I = I2 = It = IReal = J1 = J = J2 = Jt = JReal = Rav = Af = Bf = Cf = Df = Gf = Hf = If = Jf = G4ThreeVector(0, 0, 0);

 A1loc = Aloc = A2loc = Atloc = ARealloc = B1loc = Bloc = B2loc = Btloc = BRealloc = C1loc = Cloc = C2loc = Ctloc = CRealloc = D1loc = Dloc = D2loc = Dtloc = DRealloc = E1loc = Eloc = E2loc = Etloc = ERealloc = F1loc = Floc = F2loc = Ftloc = FRealloc = G1loc = Gloc = G2loc = Gtloc = GRealloc = H1loc = Hloc = H2loc = Htloc = HRealloc = I1loc = Iloc = I2loc = Itloc = IRealloc = J1loc = Jloc = J2loc = Jtloc = JRealloc = Ravloc = Afloc = Bfloc = Cfloc = Dfloc = Gfloc = Hfloc = Ifloc = Jfloc = Af1 = Bf1 = Cf1 = Df1 = Gf1 = Hf1 = If1 = Jf1 = G4ThreeVector(0, 0, 0);

 fTrack1 = fTrack2 = 0.;

 fClusterSizeXPixMin = fClusterSizeXPixMax = fClusterSizeXPixROCOfMin = fClusterSizeXPixROCOfMax = fClusterSizeYPixMin = fClusterSizeYPixMax = fClusterSizeYPixROCOfMin = fClusterSizeYPixROCOfMax = 0;

 // get event ID
 G4int evtNb = evt->GetEventID();
 //G4cout 
    //<< "\n Event ID = " 
    //<< evtNb
    //<< G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*)
{
 fRunAction->AddEnergyDeposit1(fEnergyDeposit1);
 fRunAction->AddEnergyDeposit2(fEnergyDeposit2);
 fRunAction->AddEnergyDepositL1(fEnergyDepositL1);

 G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 analysisManager->SetVerboseLevel(1);
 analysisManager->FillH1(1, fEnergyDeposit1);
 analysisManager->FillH1(2, fEnergyDeposit2);
 analysisManager->FillH1(3, fEnergySecondary1);
 analysisManager->FillH1(4, fEnergySecondary2);
 analysisManager->FillH1(5, fEnergyTertiary1);
 analysisManager->FillH1(6, fEnergyTertiary2);
 analysisManager->FillH1(7, fEnergyDeposit1+fEnergySecondary1);
 analysisManager->FillH1(8, fEnergyDeposit2+fEnergySecondary2);
 analysisManager->FillH1(21, fSecondaryxPolarization);
 analysisManager->FillH1(22, fSecondaryyPolarization);
 analysisManager->FillH1(23, fSecondaryzPolarization);
 analysisManager->FillH1(24, fTertiaryxPolarization);
 analysisManager->FillH1(25, fTertiaryyPolarization);
 analysisManager->FillH1(26, fTertiaryzPolarization);
 analysisManager->FillH1(27, fPrimaryTrackLength);
 analysisManager->FillH1(28, fSecondaryTrackLength);
 analysisManager->FillH1(29, fTertiaryTrackLength);
 analysisManager->FillH1(30, fSecondaryDetTrackLength);
 analysisManager->FillH1(31, fTrack2);

 analysisManager->FillH1(301, fEnergyDepositL1);
 analysisManager->FillH1(302, fEnergySecondaryL1);
 analysisManager->FillH1(303, fEnergyDepositL1+fEnergySecondaryL1);

 analysisManager->FillH1(304, fEnergyDepositL2);
 analysisManager->FillH1(305, fEnergySecondaryL2);
 analysisManager->FillH1(306, fEnergyDepositL2+fEnergySecondaryL2);

 analysisManager->FillH1(307, fEnergyDepositL3);
 analysisManager->FillH1(308, fEnergySecondaryL3);
 analysisManager->FillH1(309, fEnergyDepositL3+fEnergySecondaryL3);

 analysisManager->FillH1(310, fEnergyDepositL4);
 analysisManager->FillH1(311, fEnergySecondaryL4);
 analysisManager->FillH1(312, fEnergyDepositL4+fEnergySecondaryL4);

 analysisManager->FillH1(313, fEnergyDepositL5);
 analysisManager->FillH1(314, fEnergySecondaryL5);
 analysisManager->FillH1(315, fEnergyDepositL5+fEnergySecondaryL5);

 analysisManager->FillH1(316, fEnergyDepositL6);
 analysisManager->FillH1(317, fEnergySecondaryL6);
 analysisManager->FillH1(318, fEnergyDepositL6+fEnergySecondaryL6);

 analysisManager->FillH1(319, fEnergyDepositL7);
 analysisManager->FillH1(320, fEnergySecondaryL7);
 analysisManager->FillH1(321, fEnergyDepositL7+fEnergySecondaryL7);

 analysisManager->FillH1(322, fEnergyDepositL8);
 analysisManager->FillH1(323, fEnergySecondaryL8);
 analysisManager->FillH1(324, fEnergyDepositL8+fEnergySecondaryL8);




 for (G4int i = 52; i < 60; i = i + 1) {
   analysisManager->FillH1(i, fEnergyStrip1a[i+453]);
 }
 for (G4int j = 60; j < 68; j = j + 1) {
   analysisManager->FillH1(j, fEnergyStrip1b[j+445]);
 }
 for (G4int l = 68; l < 76; l = l + 1) {
   analysisManager->FillH1(l, fEnergyStrip2a[l+437]);
 }
 for (G4int n = 76; n < 84; n = n + 1) {
   analysisManager->FillH1(n, fEnergyStrip2b[n+429]);
 }

 // fill ntuple and charge per strip histograms
 for (G4int k = 0; k < 1017; k++) {
   /*analysisManager->FillNtupleDColumn(k, fEnergyStrip1a[k]);
   analysisManager->FillNtupleDColumn(k+1017, fEnergyStrip1b[k]);
   analysisManager->FillNtupleDColumn(k+2*1017, fEnergyStrip2a[k]);
   analysisManager->FillNtupleDColumn(k+3*1017, fEnergyStrip2b[k]); */

   //charge collected per strip in number of electrons [e]
   fChargeStrip1a[k] = (fEnergyStrip1a[k])/(3.67*eV);
   fChargeStrip1b[k] = (fEnergyStrip1b[k])/(3.67*eV);
   fChargeStrip2a[k] = (fEnergyStrip2a[k])/(3.67*eV);
   fChargeStrip2b[k] = (fEnergyStrip2b[k])/(3.67*eV);

   fCS1a = fChargeStrip1a[k];
   //G4cout 
     //<< "\n fCS1a = "
     //<< fCS1a;
   fCS1b = fChargeStrip1b[k];
   //G4cout 
     //<< "\n fCS1b = "
     //<< fCS1b;
   fCS2a = fChargeStrip2a[k];
   //G4cout 
     //<< "\n fCS2a = "
     //<< fCS2a;
   fCS2b = fChargeStrip2b[k];
   //G4cout 
     //<< "\n fCS2b = "
     //<< fCS2b;
 
   if (fCS1a >= 5000)  {
      fRunAction->AddNbHitsStrip1a(k);
      fHitSensor1 = 1;
      //G4cout 
         //<< "\n Strip hit: Sensor: 1 Row: a Strip: "
    	 //<< k
    	 //<< G4endl;
      //G4cout 
        //<< "\n NbHitsStrip1a advanced"
        //<< G4endl;

      fChargeStrip1 = fChargeStrip1 + 5000;
      totalNbHitStrips++;
   }
   if (fCS1b >= 5000)  {
      fRunAction->AddNbHitsStrip1b(k);
      fHitSensor1 = 1;
      fHitsMultiplicity1b ++;
      //G4cout 
         //<< "\n Strip hit: Sensor: 1 Row: b Strip: "
    	 //<< k
    	 //<< G4endl;
      //G4cout 
        //<< "\n NbHitsStrip1b advanced"
        //<< G4endl;

      fChargeStrip1 = fChargeStrip1 + 5000;
      totalNbHitStrips++;
   }
   if (fCS2a >= 5000)  {
      fRunAction->AddNbHitsStrip2a(k);
      fHitSensor2 = 1;
      //G4cout 
         //<< "\n Strip hit: Sensor: 2 Row: a Strip: "
    	 //<< k
    	 //<< G4endl;
      //G4cout 
        //<< "\n NbHitsStrip2a advanced"
        //<< G4endl;

      fChargeStrip2 = fChargeStrip2 + 5000;
      totalNbHitStrips++;
   }
   if (fCS2b >= 5000)  {
      fRunAction->AddNbHitsStrip2b(k);
      fHitSensor2 = 1;
      fHitsMultiplicity2b ++;
      //G4cout 
         //<< "\n Strip hit: Sensor: 2 Row: b Strip: "
    	 //<< k
    	 //<< G4endl;
      //G4cout 
        //<< "\n NbHitsStrip2b advanced"
        //<< G4endl;

      fChargeStrip2 = fChargeStrip2 + 5000;
      totalNbHitStrips++;
   }
 }

 //G4cout 
    //<< "\n totalNbHitStrips = "
    //<< totalNbHitStrips
    //<< G4endl;

 //charge weights and entrance positions of primaries in strip sensors per event
 for (G4int k = 0; k < 1017; k++) {
   fCS1a = fChargeStrip1a[k];
   //G4cout 
     //<< "\n fCS1a = "
     //<< fCS1a;
   fCS1b = fChargeStrip1b[k];
   //G4cout 
     //<< "\n fCS1b = "
     //<< fCS1b;
   fCS2a = fChargeStrip2a[k];
   //G4cout 
     //<< "\n fCS2a = "
     //<< fCS2a;
   fCS2b = fChargeStrip2b[k];
   //G4cout 
     //<< "\n fCS2b = "
     //<< fCS2b;

  if (fCS1a >= 5000)  {
      fWeightStrip1a[k] = 5000/fChargeStrip1;

      fStripCenterNRX = Strip1Length/2;
      fStripCenterNRY = (k-508-0.5)*StripPitch;
      fStripCenterNRZ = -Dist/2 - Det1SizeZ;

      fStripCenterX = fStripCenterNRX;
      fStripCenterY = fStripCenterNRY*cos(XangleDUT) + fStripCenterNRZ*sin(XangleDUT);
      fStripCenterZ = -fStripCenterNRY*sin(XangleDUT) + fStripCenterNRZ*cos(XangleDUT);

      E = E + G4ThreeVector((fStripCenterX*5000/fChargeStrip1), (fStripCenterY*5000/fChargeStrip1), (fStripCenterZ*5000/fChargeStrip1));

      //G4cout 
         //<< "\n fStripCenterX = "
         //<< fStripCenterX;
      //G4cout 
         //<< "\n fStripCenterY = "
         //<< fStripCenterY;
      //G4cout 
         //<< "\n fStripCenterZ = "
         //<< fStripCenterZ;

      //G4cout 
         //<< "\n fCS1a = "
         //<< fCS1a;
      //G4cout 
         //<< "\n fChargeStrip1 = "
         //<< fChargeStrip1;
      //G4cout 
         //<< "\n fWeightStrip1a[k] = "
         //<< fWeightStrip1a[k];

      //G4cout 
         //<< "\n fPointDet1ent = "
         //<< fPointDet1ent;
   }
   if (fCS1b >= 5000)  {
      fWeightStrip1b[k] = 5000/fChargeStrip1;
      fStripCenterNRX = -Strip1Length/2;
      fStripCenterNRY = (k-508-0.5)*StripPitch;
      fStripCenterNRZ = -Dist/2 - Det1SizeZ;

      fStripCenterX = fStripCenterNRX;
      fStripCenterY = fStripCenterNRY*cos(XangleDUT) + fStripCenterNRZ*sin(XangleDUT);
      fStripCenterZ = -fStripCenterNRY*sin(XangleDUT) + fStripCenterNRZ*cos(XangleDUT);

      E = E + G4ThreeVector((fStripCenterX*5000/fChargeStrip1), (fStripCenterY*5000/fChargeStrip1), (fStripCenterZ*5000/fChargeStrip1));

      //G4cout 
         //<< "\n fStripCenterX = "
         //<< fStripCenterX;
      //G4cout 
         //<< "\n fStripCenterY = "
         //<< fStripCenterY;
      //G4cout 
         //<< "\n fStripCenterZ = "
         //<< fStripCenterZ;

      //G4cout 
         //<< "\n fCS1b = "
         //<< fCS1b;
      //G4cout 
         //<< "\n fChargeStrip1 = "
         //<< fChargeStrip1;
      //G4cout 
         //<< "\n fWeightStrip1b[k] = "
         //<< fWeightStrip1b[k];

      //G4cout 
         //<< "\n fPointDet1ent = "
         //<< fPointDet1ent;
   }
   if (fCS2a >= 5000)  {
      fWeightStrip2a[k] = 5000/fChargeStrip2;

      fStripCenterNRX = Strip2Length/2;
      fStripCenterNRY = (k-508-0.5)*StripPitch;
      fStripCenterNRZ = Dist/2;

      fStripCenterX = fStripCenterNRX;
      fStripCenterY = fStripCenterNRY*cos(XangleDUT) + fStripCenterNRZ*sin(XangleDUT);
      fStripCenterZ = -fStripCenterNRY*sin(XangleDUT) + fStripCenterNRZ*cos(XangleDUT);

      F = F + G4ThreeVector((fStripCenterX*5000/fChargeStrip2), (fStripCenterY*5000/fChargeStrip2), (fStripCenterZ*5000/fChargeStrip2));

      //G4cout 
         //<< "\n fStripCenterX = "
         //<< fStripCenterX;
      //G4cout 
         //<< "\n fStripCenterY = "
         //<< fStripCenterY;
      //G4cout 
         //<< "\n fStripCenterZ = "
         //<< fStripCenterZ;

      //G4cout 
         //<< "\n fCS2a = "
         //<< fCS2a;
      //G4cout 
         //<< "\n fChargeStrip2 = "
         //<< fChargeStrip2;
      //G4cout 
         //<< "\n fWeightStrip2a[k] = "
         //<< fWeightStrip2a[k];

      //G4cout 
         //<< "\n fPointDet2ent = "
         //<< fPointDet2ent;
   }
   if (fCS2b >= 5000)  {
      fWeightStrip2b[k] = 5000/fChargeStrip2;

      fStripCenterNRX = -Strip2Length/2;
      fStripCenterNRY = (k-508-0.5)*StripPitch;
      fStripCenterNRZ = Dist/2;


      fStripCenterX = fStripCenterNRX;
      fStripCenterY = fStripCenterNRY*cos(XangleDUT) + fStripCenterNRZ*sin(XangleDUT);
      fStripCenterZ = -fStripCenterNRY*sin(XangleDUT) + fStripCenterNRZ*cos(XangleDUT);

      F = F + G4ThreeVector((fStripCenterX*5000/fChargeStrip2), (fStripCenterY*5000/fChargeStrip2), (fStripCenterZ*5000/fChargeStrip2));

      //G4cout 
         //<< "\n fStripCenterX = "
         //<< fStripCenterX;
      //G4cout 
         //<< "\n fStripCenterY = "
         //<< fStripCenterY;
      //G4cout 
         //<< "\n fStripCenterZ = "
         //<< fStripCenterZ;

      //G4cout 
         //<< "\n fCS2b = "
         //<< fCS2b;
      //G4cout 
         //<< "\n fChargeStrip2 = "
         //<< fChargeStrip2;
      //G4cout 
         //<< "\n fWeightStrip2b[k] = "
         //<< fWeightStrip2b[k];

      //G4cout 
         //<< "\n fPointDet2ent = "
         //<< fPointDet2ent;
   }
 }

 if ((fHitSensor1 == 1) && (fHitSensor2 == 1))  {
      fRunAction->AddNbCoinc();
 }

 // fill ntuple and charge per pixel histograms
 G4int kec = 0;
 for (G4int ke1 = 0; ke1 < 17; ke1++) {
    for (G4int ke2 = 0; ke2 < 17; ke2++) {
 	for (G4int ke3 = 0; ke3 < 53; ke3++) {
 	   for (G4int ke4 = 0; ke4 < 81; ke4++) {
		//analysisManager->FillNtupleDColumn(kec+8*1017, fEnergyPixel[ke1][ke2][ke3][ke4]);

   		//charge collected per pixel in number of electrons [e]

		enPix = fEnergyPixel[ke1][ke2][ke3][ke4];
		if (enPix > 0)   {
		  /*G4cout
		    << "\n enPix = "
		    << G4BestUnit(enPix, "Energy")
		    << G4endl;*/
		}

   		fChargePixel[ke1][ke2][ke3][ke4] = enPix/(3.67*eV);

   		fChargePix = fChargePixel[ke1][ke2][ke3][ke4];
   		/*G4cout 
     		  << "\n fChargePix = "
     		  << fChargePix
     		  << " e" 
                  << G4endl; */
 
   		if (fChargePix >= fPixThreshold)  {
      		   fRunAction->AddNbHitsPixel(ke1, ke2, ke3, ke4);
      		   /*G4cout 
        	      << "\n NbHitsPixel advanced. fChargePix = "
		      << fChargePix
     		      << " e" 
        	      << G4endl;*/
		   kec++;
   		}
 	   }
 	}
    }
 }

 //analysisManager->AddNtupleRow();

 analysisManager->FillH1(189, fHitsMultiplicity1b);
 //G4cout 
    //<< "\n hit multiplicity of sensor 1: "
    //<< fHitsMultiplicity1b
    //<< G4endl;
 analysisManager->FillH1(190, fHitsMultiplicity2b);
 //G4cout 
    //<< "\n hit multiplicity of sensor 2: "
    //<< fHitsMultiplicity2b
    //<< G4endl;


 fCharge1 = fEnergyDeposit1/(3.67*eV);
 fCharge2 = fEnergyDeposit2/(3.67*eV);

 analysisManager->FillH1(155, fCharge1);
 analysisManager->FillH1(156, fCharge2);


 for (G4int ke1 = 0; ke1 < 17; ke1++) {
   for (G4int ke2 = 0; ke2 < 17; ke2++) {
     for (G4int ke3 = 0; ke3 < 53; ke3++) {
       for (G4int ke4 = 0; ke4 < 81; ke4++) {
         //charge collected per pixel in number of electrons [e]
         fChargePixel[ke1][ke2][ke3][ke4] = (fEnergyPixel[ke1][ke2][ke3][ke4])/(3.67*eV);
	 fEnergyPix = fEnergyPixel[ke1][ke2][ke3][ke4];
	 /*G4cout
	    << "\n fEnergyPix = "
	    << G4BestUnit(fEnergyPix, "Energy")
	    << G4endl; */

         fChargePix = fChargePixel[ke1][ke2][ke3][ke4];
	 if (fChargePix > 0)  {
	    /*G4cout
	       << "\n fChargePix = "
	       << fChargePix
    	       << "\n Module: "
	       << ke1
	       << "\n ROC: "
    	       << ke2
    	       << "\n row: "
    	       << ke3
    	       << "\n col: "
    	       << ke4
	       << G4endl; */
	 }
 
         if (fChargePix >= fPixThreshold)  {
	   fChargeModule[ke1] = fChargeModule[ke1] + fChargePix;
 	   /*G4cout 
    	      << "\n Pixel hit: Module: "
	      << ke1
	      << "ROC: "
    	      << ke2
    	      << " row: "
    	      << ke3
    	      << " col: "
    	      << ke4
    	      << " fChargeModule[ke1] = "
    	      << fChargeModule[ke1]
    	      << G4endl; */
         }
       }
     }
   }
 }

 fChargeL1 = fChargeModule[2];
 if (fChargeL1 > 0)   {
    analysisManager->FillH1(357, fChargeL1);
 }

 //charge weights and entrance positions of primaries in BPIX modules per event
 for (G4int ke1 = 0; ke1 < 17; ke1++) {
  for (G4int ke2 = 0; ke2 < 17; ke2++) {
    for (G4int ke3 = 0; ke3 < 53; ke3++) {
      for (G4int ke4 = 0; ke4 < 81; ke4++) {
        fChargePix = fChargePixel[ke1][ke2][ke3][ke4];
	if (fChargePix > 0)   {
	}
  
        if (fChargePix >= fPixThreshold)  {
 	  fHitPixelDet[ke1] = 1;
          fWeightPixel[ke1][ke2][ke3][ke4] = fChargePix/fChargeModule[ke1];

	  if ((ke2 == 1) || (ke2 == 9))  {
	    if (ke3 == 1)  {
	      fPixelCenterNRY = 3*54*pixelY + 53*pixelY;
	    }
	    if ((ke3 >= 2) && (ke3<=51))  {
	      fPixelCenterNRY = 3*54*pixelY + (53-ke3)*pixelY + pixelY/2;
	    }
	    if (ke3 == 52)  {
	      fPixelCenterNRY = 3*54*pixelY + pixelY;
	    }
	  }
	  if ((ke2 == 2) || (ke2 == 10))  {
	    if (ke3 == 1)  {
	      fPixelCenterNRY = 2*54*pixelY + 53*pixelY;
	    }
	    if ((ke3 >= 2) && (ke3<=51))  {
	      fPixelCenterNRY = 2*54*pixelY + (53-ke3)*pixelY + pixelY/2;
	    }
	    if (ke3 == 52)  {
	      fPixelCenterNRY = 2*54*pixelY + pixelY;
	    }
	  }
	  if ((ke2 == 3) || (ke2 == 11))  {
	    if (ke3 == 1)  {
	      fPixelCenterNRY = 1*54*pixelY + 53*pixelY;
	    }
	    if ((ke3 >= 2) && (ke3<=51))  {
	      fPixelCenterNRY = 1*54*pixelY + (53-ke3)*pixelY + pixelY/2;
	    }
	    if (ke3 == 52)  {
	      fPixelCenterNRY = 1*54*pixelY + pixelY;
	    }
	  }
	  if ((ke2 == 4) || (ke2 == 12))  {
	    if (ke3 == 1)  {
	      fPixelCenterNRY = 0*54*pixelY + 53*pixelY;
	    }
	    if ((ke3 >= 2) && (ke3<=51))  {
	      fPixelCenterNRY = 0*54*pixelY + (53-ke3)*pixelY + pixelY/2;
	    }
	    if (ke3 == 52)  {
	      fPixelCenterNRY = 0*54*pixelY + pixelY;
	    }
	  }
	  if ((ke2 == 5) || (ke2 == 13))  {
	    if (ke3 == 1)  {
	      fPixelCenterNRY = -0*54*pixelY - pixelY;
	    }
	    if ((ke3 >= 2) && (ke3<=51))  {
	      fPixelCenterNRY = -0*54*pixelY - ke3*pixelY - pixelY/2;
	    }
	    if (ke3 == 52)  {
	      fPixelCenterNRY = -0*54*pixelY - 53*pixelY;
	    }
	  }
	  if ((ke2 == 6) || (ke2 == 14))  {
	    if (ke3 == 1)  {
	      fPixelCenterNRY = -1*54*pixelY - pixelY;
	    }
	    if ((ke3 >= 2) && (ke3<=51))  {
	      fPixelCenterNRY = -1*54*pixelY - ke3*pixelY - pixelY/2;
	    }
	    if (ke3 == 52)  {
	      fPixelCenterNRY = -1*54*pixelY - 53*pixelY;
	    }
	  }
	  if ((ke2 == 7) || (ke2 == 15))  {
	    if (ke3 == 1)  {
	      fPixelCenterNRY = -2*54*pixelY - pixelY;
	    }
	    if ((ke3 >= 2) && (ke3<=51))  {
	      fPixelCenterNRY = -2*54*pixelY - ke3*pixelY - pixelY/2;
	    }
	    if (ke3 == 52)  {
	      fPixelCenterNRY = -2*54*pixelY - 53*pixelY;
	    }
	  }
	  if ((ke2 == 8) || (ke2 == 16))  {
	    if (ke3 == 1)  {
	      fPixelCenterNRY = -3*54*pixelY - pixelY;
	    }
	    if ((ke3 >= 2) && (ke3<=51))  {
	      fPixelCenterNRY = -3*54*pixelY - ke3*pixelY - pixelY/2;
	    }
	    if (ke3 == 52)  {
	      fPixelCenterNRY = -3*54*pixelY - 53*pixelY;
	    }
	  }

	  if ((ke2>=1) && (ke2<=8))  {
	    if (ke4 == 80) {
	      fPixelCenterNRX = pixelX;
            }
	    if ((ke4 >= 1) && (ke4 <= 79)) {
	      fPixelCenterNRX = (81 - ke4)*pixelX + pixelX/2;
            }
	  }
	  if ((ke2>=9) && (ke2<=16))  {
	    if (ke4 == 80) {
	      fPixelCenterNRX = -80*pixelX;
            }
	    if ((ke4 >= 1) && (ke4 <= 79)) {
	      fPixelCenterNRX = -(ke4 - 1)*pixelX - pixelX/2;
            }
	  }

	  if ((ke1 == 1) || (ke1 == 3) || (ke1 == 5) || (ke1 == 7) || (ke1 == 9) || (ke1 == 11) || (ke1 == 13) || (ke1 == 15))  {
	     fPixelCenterNRX = fPixelCenterNRX + BPIXSizeX/2;
	  }
	  if ((ke1 == 2) || (ke1 == 4) || (ke1 == 6) || (ke1 == 8) || (ke1 == 10) || (ke1 == 12) || (ke1 == 14) || (ke1 == 16))  {
	     fPixelCenterNRX = fPixelCenterNRX - BPIXSizeX/2;
	  }

	  fPixelCenterNRZ = BPIXSizeZ/2;

  	  fPixelCenterNRX1 = fPixelCenterNRX;
  	  fPixelCenterNRY1 = fPixelCenterNRY*cos(XangleBPIX) - fPixelCenterNRZ*sin(XangleBPIX);
   	  fPixelCenterNRZ1 = fPixelCenterNRY*sin(XangleBPIX) + fPixelCenterNRZ*cos(XangleBPIX);

  	  fPixelCenterX = fPixelCenterNRX1*cos(YangleBPIX) + fPixelCenterNRZ1*sin(YangleBPIX);
  	  fPixelCenterY = fPixelCenterNRY1;
  	  fPixelCenterNTZ = -fPixelCenterNRX1*sin(YangleBPIX) + fPixelCenterNRZ1*cos(YangleBPIX);

	  if ((ke1 == 1) || (ke1 == 2))  {
	     fPixelCenterZ = fPixelCenterNTZ + poszBPIX12;
	  }
	  if ((ke1 == 3) || (ke1 == 4))  {
	     fPixelCenterZ = fPixelCenterNTZ + poszBPIX34;
	  }
	  if ((ke1 == 5) || (ke1 == 6))  {
	     fPixelCenterZ = fPixelCenterNTZ + poszBPIX56;
	  }
	  if ((ke1 == 7) || (ke1 == 8))  {
	     fPixelCenterZ = fPixelCenterNTZ + poszBPIX78;
	  }
	  if ((ke1 == 9) || (ke1 == 10))  {
	     fPixelCenterZ = fPixelCenterNTZ + poszBPIX910;
	  }
	  if ((ke1 == 11) || (ke1 == 12))  {
	     fPixelCenterZ = fPixelCenterNTZ + poszBPIX1112;
	  }
	  if ((ke1 == 13) || (ke1 == 14))  {
	     fPixelCenterZ = fPixelCenterNTZ + poszBPIX1314;
	  }
	  if ((ke1 == 15) || (ke1 == 16))  {
	     fPixelCenterZ = fPixelCenterNTZ + poszBPIX1516;
	  }

	  if ((ke1 == 1) || (ke1 == 2)) {
             A = A + G4ThreeVector((fPixelCenterX*fChargePix/fChargeModule[ke1]), (fPixelCenterY*fChargePix/fChargeModule[ke1]), (fPixelCenterZ*fChargePix/fChargeModule[ke1]));
	     /*G4cout
		<< "\n A = "
	        << G4BestUnit(A, "Length")
	        << G4endl; */
	  }
	  if ((ke1 == 3) || (ke1 == 4))  {
           B = B + G4ThreeVector((fPixelCenterX*fChargePix/fChargeModule[ke1]), (fPixelCenterY*fChargePix/fChargeModule[ke1]), (fPixelCenterZ*fChargePix/fChargeModule[ke1]));
	     /*G4cout
		<< "\n B = "
	        << G4BestUnit(B, "Length")
	        << G4endl; */
	  }
	  if ((ke1 == 5) || (ke1 == 6))  {
	   C = C + G4ThreeVector((fPixelCenterX*fChargePix/fChargeModule[ke1]), (fPixelCenterY*fChargePix/fChargeModule[ke1]), (fPixelCenterZ*fChargePix/fChargeModule[ke1]));
	     /*G4cout
		<< "\n C = "
	        << G4BestUnit(C, "Length")
	        << G4endl; */
	  }
	  if ((ke1 == 7) || (ke1 == 8))  {
             D = D + G4ThreeVector((fPixelCenterX*fChargePix/fChargeModule[ke1]), (fPixelCenterY*fChargePix/fChargeModule[ke1]), (fPixelCenterZ*fChargePix/fChargeModule[ke1]));
	     /*G4cout
		<< "\n D = "
	        << G4BestUnit(D, "Length")
	        << G4endl; */
	  }
	  if ((ke1 == 9) || (ke1 == 10))  {
             G = G + G4ThreeVector((fPixelCenterX*fChargePix/fChargeModule[ke1]), (fPixelCenterY*fChargePix/fChargeModule[ke1]), (fPixelCenterZ*fChargePix/fChargeModule[ke1]));
	     /*G4cout
		<< "\n G = "
	        << G4BestUnit(G, "Length")
	        << G4endl; */
	  }
	  if ((ke1 == 11) || (ke1 == 12))  {
             H = H + G4ThreeVector((fPixelCenterX*fChargePix/fChargeModule[ke1]), (fPixelCenterY*fChargePix/fChargeModule[ke1]), (fPixelCenterZ*fChargePix/fChargeModule[ke1]));
	     /*G4cout
		<< "\n H = "
	        << G4BestUnit(H, "Length")
	        << G4endl; */
	  }
	  if ((ke1 == 13) || (ke1 == 14))  {
             I = I + G4ThreeVector((fPixelCenterX*fChargePix/fChargeModule[ke1]), (fPixelCenterY*fChargePix/fChargeModule[ke1]), (fPixelCenterZ*fChargePix/fChargeModule[ke1]));
	     /*G4cout
		<< "\n I = "
	        << G4BestUnit(I, "Length")
	        << G4endl; */
	  }
	  if ((ke1 == 15) || (ke1 == 16))  {
             J = J + G4ThreeVector((fPixelCenterX*fChargePix/fChargeModule[ke1]), (fPixelCenterY*fChargePix/fChargeModule[ke1]), (fPixelCenterZ*fChargePix/fChargeModule[ke1]));
	     /*G4cout
		<< "\n J = "
	        << G4BestUnit(J, "Length")
	        << G4endl; */
	  }
        }
      }
    }
  }
 }

 for (G4int i2 = 157; i2 < 165; i2 = i2 + 1) {
   analysisManager->FillH1(i2, fChargeStrip1a[i2+348]);
 }
 for (G4int p2a = 165; p2a < 173; p2a = p2a + 1) {
   analysisManager->FillH1(p2a, fChargeStrip1b[p2a+340]);
 }
 for (G4int l2 = 173; l2 < 181; l2 = l2 + 1) {
   analysisManager->FillH1(l2, fChargeStrip2a[l2+332]);
 }
 for (G4int n2 = 181; n2 < 189; n2 = n2 + 1) {
   analysisManager->FillH1(n2, fChargeStrip2b[n2+324]);
 }

 ftheta = acos(fMomDir1x*fMomDir2x + fMomDir1y*fMomDir2y + fMomDir1z*fMomDir2z);
 analysisManager->FillH1(116, ftheta);

 //G4cout
    //<< "\n The deflection angle is: "
    //<< ftheta
    //<< " rad"
    //<< G4endl;


 for (G4int kea = 0; kea < 17; kea++) {
    for (G4int ke2 = 0; ke2 < 17; ke2++) {
      for (G4int ke3 = 0; ke3 < 53; ke3++) {
        for (G4int ke4 = 0; ke4 < 81; ke4++) {
          //charge collected per pixel in number of electrons [e]
          fChargePixel[kea][ke2][ke3][ke4] = (fEnergyPixel[kea][ke2][ke3][ke4])/(3.67*eV);

          fChargePix = fChargePixel[kea][ke2][ke3][ke4];
 
          if (fChargePix >= fPixThreshold)  {
	    fChargeModule[kea] = fChargeModule[kea] + fChargePix;
          }
	}
     }
   }
 }

 for (G4int kea = 0; kea < 17; kea++) {
    for (G4int ke2 = 0; ke2 < 17; ke2++) {
      for (G4int ke3 = 0; ke3 < 53; ke3++) {
        for (G4int ke4 = 0; ke4 < 81; ke4++) {
          fChargePix = fChargePixel[kea][ke2][ke3][ke4];
 
          if (fChargePix >= fPixThreshold)  {
	     mod = kea;
	     if ((ke2 >= 9) && (ke2 <= 16))   {
	        row = 80 + ke4;
	        col = (16 - ke2)*52 + 53 - ke3;
	     }
	     if ((ke2 >= 1) && (ke2 <= 8))   {
	        row = ke4;
	        col = (8 - ke2)*52 + 53 - ke3;
	     }
	     fClusterOccupancy[mod][row][col]++;
   	     if (kea > 0)  {
      	        analysisManager->FillH2((kea+20), col, row);
   	     }
          }
	}
     }
   }

 }

 for (G4int kea = 0; kea < 17; kea++) {
    for (G4int ke2 = 0; ke2 < 17; ke2++) {
      for (G4int ke3 = 0; ke3 < 53; ke3++) {
        for (G4int ke4 = 0; ke4 < 81; ke4++) {
          //charge collected per pixel in number of electrons [e]

          fChargePix = fChargePixel[kea][ke2][ke3][ke4];
 
          if (fChargePix >= fPixThreshold)  {
	    fHitsMultiplicityPix[kea] = fHitsMultiplicityPix[kea] + 1;

	    if (fClusterSizeXPixMin == 0)  {
	       fClusterSizeXPixMin = ke4;
	       fClusterSizeXPixROCOfMin = ke2;
	    }
	    fClusterSizeXPixMax = ke4;
	    fClusterSizeXPixROCOfMax = ke2;

	    if (fClusterSizeYPixMin == 0)  {
	       fClusterSizeYPixMin = ke3;
	       fClusterSizeYPixROCOfMin = ke2;
	    }
	    fClusterSizeYPixMax = ke3;
	    fClusterSizeYPixROCOfMax = ke2;
          }
	}
     }
   }

   if (fClusterSizeXPixROCOfMax > fClusterSizeXPixROCOfMin)  {
      fClusterSizeXPix[kea] = fClusterSizeXPixMax + 81 - fClusterSizeXPixMin;
   }
   else  {
      fClusterSizeXPix[kea] = fClusterSizeXPixMax + 1 - fClusterSizeXPixMin;
   }

   if (fClusterSizeYPixROCOfMax > fClusterSizeYPixROCOfMin)  {
      fClusterSizeYPix[kea] = fClusterSizeYPixMax + 53 +52*(fClusterSizeYPixROCOfMax - fClusterSizeYPixROCOfMin - 1) - fClusterSizeYPixMin;
   }
   else  {
      fClusterSizeYPix[kea] = fClusterSizeYPixMax + 1 - fClusterSizeYPixMin;
   }
 }
 
 if (fClusterSizeXPixMin > 0)  {
   analysisManager->FillH1(325, fClusterSizeXPix[1]);
   analysisManager->FillH1(326, fClusterSizeXPix[2]);
   analysisManager->FillH1(327, fClusterSizeXPix[3]);
   analysisManager->FillH1(328, fClusterSizeXPix[4]);
   analysisManager->FillH1(329, fClusterSizeXPix[5]);
   analysisManager->FillH1(330, fClusterSizeXPix[6]);
   analysisManager->FillH1(331, fClusterSizeXPix[7]);
   analysisManager->FillH1(332, fClusterSizeXPix[8]);
   analysisManager->FillH1(333, fClusterSizeXPix[9]);
   analysisManager->FillH1(334, fClusterSizeXPix[10]);
   analysisManager->FillH1(335, fClusterSizeXPix[11]);
   analysisManager->FillH1(336, fClusterSizeXPix[12]);
   analysisManager->FillH1(337, fClusterSizeXPix[13]);
   analysisManager->FillH1(338, fClusterSizeXPix[14]);
   analysisManager->FillH1(339, fClusterSizeXPix[15]);
   analysisManager->FillH1(340, fClusterSizeXPix[16]);

   analysisManager->FillH1(341, fClusterSizeYPix[1]);
   analysisManager->FillH1(342, fClusterSizeYPix[2]);
   analysisManager->FillH1(343, fClusterSizeYPix[3]);
   analysisManager->FillH1(344, fClusterSizeYPix[4]);
   analysisManager->FillH1(345, fClusterSizeYPix[5]);
   analysisManager->FillH1(346, fClusterSizeYPix[6]);
   analysisManager->FillH1(347, fClusterSizeYPix[7]);
   analysisManager->FillH1(348, fClusterSizeYPix[8]);
   analysisManager->FillH1(349, fClusterSizeYPix[9]);
   analysisManager->FillH1(350, fClusterSizeYPix[10]);
   analysisManager->FillH1(351, fClusterSizeYPix[11]);
   analysisManager->FillH1(352, fClusterSizeYPix[12]);
   analysisManager->FillH1(353, fClusterSizeYPix[13]);
   analysisManager->FillH1(354, fClusterSizeYPix[14]);
   analysisManager->FillH1(355, fClusterSizeYPix[15]);
   analysisManager->FillH1(356, fClusterSizeYPix[16]);

   if (fHitsMultiplicityPix[1] > 0)   {
      analysisManager->FillH1(197, fHitsMultiplicityPix[1]);
   }
   if (fHitsMultiplicityPix[2] > 0)   {
      analysisManager->FillH1(198, fHitsMultiplicityPix[2]);
   }
   if (fHitsMultiplicityPix[3] > 0)   {
      analysisManager->FillH1(199, fHitsMultiplicityPix[3]);
   }
   if (fHitsMultiplicityPix[4] > 0)   {
      analysisManager->FillH1(200, fHitsMultiplicityPix[4]);
   }
   if (fHitsMultiplicityPix[5] > 0)   {
      analysisManager->FillH1(201, fHitsMultiplicityPix[5]);
   }
   if (fHitsMultiplicityPix[6] > 0)   {
      analysisManager->FillH1(202, fHitsMultiplicityPix[6]);
   }
   if (fHitsMultiplicityPix[7] > 0)   {
      analysisManager->FillH1(203, fHitsMultiplicityPix[7]);
   }
   if (fHitsMultiplicityPix[8] > 0)   {
      analysisManager->FillH1(204, fHitsMultiplicityPix[8]);
   }
   if (fHitsMultiplicityPix[9] > 0)   {
      analysisManager->FillH1(205, fHitsMultiplicityPix[9]);
   }
   if (fHitsMultiplicityPix[10] > 0)   {
      analysisManager->FillH1(206, fHitsMultiplicityPix[10]);
   }
   if (fHitsMultiplicityPix[11] > 0)   {
      analysisManager->FillH1(207, fHitsMultiplicityPix[11]);
   }
   if (fHitsMultiplicityPix[12] > 0)   {
      analysisManager->FillH1(208, fHitsMultiplicityPix[12]);
   }
   if (fHitsMultiplicityPix[13] > 0)   {
      analysisManager->FillH1(209, fHitsMultiplicityPix[13]);
   }
   if (fHitsMultiplicityPix[14] > 0)   {
      analysisManager->FillH1(210, fHitsMultiplicityPix[14]);
   }
   if (fHitsMultiplicityPix[15] > 0)   {
      analysisManager->FillH1(211, fHitsMultiplicityPix[15]);
   }
   if (fHitsMultiplicityPix[16] > 0)   {
      analysisManager->FillH1(212, fHitsMultiplicityPix[16]);
   }
 }

 Rav = (A+B+C+D+G+H+I+J)/8;
 At = A - Rav;
 Bt = B - Rav;
 Ct = C - Rav;
 Dt = D - Rav;
 Gt = G - Rav;
 Ht = H - Rav;
 It = I - Rav;
 Jt = J - Rav;

 M[0][0] = ((At.x())*(At.x()) + (Bt.x())*(Bt.x()) + (Ct.x())*(Ct.x()) + (Dt.x())*(Dt.x()) + (Gt.x())*(Gt.x()) + (Ht.x())*(Ht.x()) + (It.x())*(It.x()) + (Jt.x())*(Jt.x()))/8;
 M[0][1] = M[1][0] = ((At.x())*(At.y()) + (Bt.x())*(Bt.y()) + (Ct.x())*(Ct.y()) + (Dt.x())*(Dt.y()) + (Gt.x())*(Gt.y()) + (Ht.x())*(Ht.y()) + (It.x())*(It.y()) + (Jt.x())*(Jt.y()))/8;
 M[0][2] = M[2][0] = ((At.x())*(At.z()) + (Bt.x())*(Bt.z()) + (Ct.x())*(Ct.z()) + (Dt.x())*(Dt.z()) + (Gt.x())*(Gt.z()) + (Ht.x())*(Ht.z()) + (It.x())*(It.z()) + (Jt.x())*(Jt.z()))/8;
 M[1][1] = ((At.y())*(At.y()) + (Bt.y())*(Bt.y()) + (Ct.y())*(Ct.y()) + (Dt.y())*(Dt.y()) + (Gt.y())*(Gt.y()) + (Ht.y())*(Ht.y()) + (It.y())*(It.y()) + (Jt.y())*(Jt.y()))/8;
 M[1][2] = M[2][1] = ((At.y())*(At.z()) + (Bt.y())*(Bt.z()) + (Ct.y())*(Ct.z()) + (Dt.y())*(Dt.z()) + (Gt.y())*(Gt.z()) + (Ht.y())*(Ht.z()) + (It.y())*(It.z()) + (Jt.y())*(Jt.z()))/8;
 M[2][2] = ((At.z())*(At.z()) + (Bt.z())*(Bt.z()) + (Ct.z())*(Ct.z()) + (Dt.z())*(Dt.z()) + (Gt.z())*(Gt.z()) + (Ht.z())*(Ht.z()) + + (It.z())*(It.z()) + (Jt.z())*(Jt.z()))/8;

 //G4cout 
    //<< "\n M[0][0]: " 
    //<< M[0][0]
    //<< "\n M[0][1]: " 
    //<< M[0][1]
    //<< "\n M[0][2]: " 
    //<< M[0][2]
    //<< "\n M[1][0]: " 
    //<< M[1][0]
    //<< "\n M[1][1]: " 
    //<< M[1][1]
    //<< "\n M[1][2]: " 
    //<< M[1][2]
    //<< "\n M[2][0]: " 
    //<< M[2][0]
    //<< "\n M[2][1]: " 
    //<< M[2][1]
    //<< "\n M[2][2]: " 
    //<< M[2][2]
    //<< G4endl;

 //for (G4int wi = 0; wi < 3; wi = wi + 1) {
    //for (G4int wj = 0; wj < 3; wj = wj + 1) {
        //G4cout 
    	   //<< "\n M[i][j]: " 
           //<< M[wi][wj]
           //<< G4endl;
    //}
 //}


 //Calculation of eigenvalues of symmetric 3X3 matrix M
 p1 = (M[0][1])*(M[0][1]) + (M[0][2])*(M[0][2]) + (M[1][2])*(M[1][2]);
 traceM = M[0][0] + M[1][1] + M[2][2]; //The sum of all diagonal values
 Iu[0][0] = Iu[1][1] = Iu[2][2] = 1;
 Iu[0][1] = Iu[1][0] = Iu[0][2] = Iu[2][0] = Iu[1][2] = Iu[2][1] = 0; //Iu is the identity matrix
 if (p1 == 0) { //M is diagonal
    eig1 = M[0][0];
    eig2 = M[1][1];
    eig3 = M[2][2];
    lambda = eig1;
    if (eig2>lambda)   {
     lambda = eig2;
    }
    if (eig3>lambda)   {
     lambda = eig3;
    }
 }
 else {
    q = traceM/3;
    p2 = (M[0][0] - q)*(M[0][0] - q) + (M[1][1] - q)*(M[1][1] - q) + (M[2][2] - q)*(M[2][2] - q) + 2*p1;
    p = sqrt(p2/6);
    for (G4int bi = 0; bi < 3; bi = bi + 1) {
    	for (G4int bj = 0; bj < 3; bj = bj + 1) {
           if (p>0) {
              Ba[bi][bj] = (M[bi][bj] - q*Iu[bi][bj])/p;
           }
           else {
              Ba[bi][bj] = 0;
           }
    	}
    }
    detBa = (Ba[0][0])*((Ba[1][1])*(Ba[2][2]) - (Ba[1][2])*(Ba[2][1])) - (Ba[0][1])*((Ba[1][0])*(Ba[2][2]) - (Ba[1][2])*(Ba[2][0])) + (Ba[0][2])*((Ba[1][0])*(Ba[2][1]) - (Ba[1][1])*(Ba[2][0]));
    r = detBa/2;

    //In exact arithmetic for a symmetric matrix -1 <= r <= 1 but computation error can leave it slightly outside this range.
    if (r <= -1) {
       phi = 3.14159/3;
    }
    else if (r >= 1) {
       phi = 0;
    }
    else {
       phi = acos(r)/3;
    }

    //The eigenvalues satisfy eig3 <= eig2 <= eig1
    eig1 = q + 2*p*cos(phi);
    eig3 = q + 2*p*cos(phi + 2*3.14159/3);
    eig2 = 3*q - eig1 - eig3; //Since trace(M) = eig1 + eig2 + eig3
    lambda = eig1;
 }

 //Calculation of the eigenvector of maximum eigenvalue lambda with Gaussian elimination
 Ga[0][0] = M[0][0] - lambda;
 Ga[0][1] = Ga[1][0] = M[0][1];
 Ga[0][2] = Ga[2][0] = M[0][2];
 Ga[1][1] = M[1][1] - lambda;
 Ga[1][2] = Ga[2][1] = M[1][2];
 Ga[2][2] = M[2][2] - lambda;

 //G4cout 
    //<< "\n Before the Gaussian elimination: "   	
    //<< "\n Ga[0][0]: " 
    //<< Ga[0][0]
    //<< "\n Ga[0][1]: " 
    //<< Ga[0][1]
    //<< "\n Ga[0][2]: " 
    //<< Ga[0][2]
    //<< "\n Ga[1][0]: " 
    //<< Ga[1][0]
    //<< "\n Ga[1][1]: " 
    //<< Ga[1][1]
    //<< "\n Ga[1][2]: " 
    //<< Ga[1][2]
    //<< "\n Ga[2][0]: " 
    //<< Ga[2][0]
    //<< "\n Ga[2][1]: " 
    //<< Ga[2][1]
    //<< "\n Ga[2][2]: " 
    //<< Ga[2][2]
    //<< G4endl;
 
 G4int h = 0; //Initialization of the pivot row
 G4int k = 0; //Initialization of the pivot column
 G4int i_max;
 G4double maxGa, checkGa, cGa, Ga0, Ga1, Ga2, fa;
 while ((h<=2)&&(k<=2)) {
    //Find the k-th pivot
    i_max = h; //argmax(ii = h ... 2, abs(G[ii,k]))
    maxGa = abs(Ga[h][k]);
    for (G4int ii = h; ii < 3; ii = ii + 1) {
       checkGa = abs(Ga[ii][k]);
       if (checkGa > maxGa) {
	  maxGa = checkGa;
          i_max = ii;
       }
    }

    cGa = Ga[i_max][k];
    if (cGa == 0) { //No pivot in this column, pass to next column
       k = k+1;
    }
    else {
       //Swap rows h, i_max
       Ga0 = Ga[h][0];
       Ga1 = Ga[h][1];
       Ga2 = Ga[h][2];
       Ga[h][0] = Ga[i_max][0];
       Ga[h][1] = Ga[i_max][1];
       Ga[h][2] = Ga[i_max][2];
       Ga[i_max][0] = Ga0;
       Ga[i_max][1] = Ga1;
       Ga[i_max][2] = Ga2;

       //Do for all rows below pivot:
       for (G4int ij = h+1; ij < 3; ij = ij + 1) {
          fa = Ga[ij][k]/Ga[h][k];
          //Fill with zeros the lowest part of pivot column:
          Ga[ij][k] = 0;
          //Do for all remaining elements in current row:
          for (G4int jj = k+1; jj < 3; jj = jj + 1) {
             Ga[ij][jj] = Ga[ij][jj] - fa*Ga[h][jj];
          }
       }

       //Increase pivot row and column
       h = h+1;
       k = k+1;
    }
 }

 //G4cout 
    //<< "\n After the Gaussian elimination: "   	
    //<< "\n Ga[0][0]: " 
    //<< Ga[0][0]
    //<< "\n Ga[0][1]: " 
    //<< Ga[0][1]
    //<< "\n Ga[0][2]: " 
    //<< Ga[0][2]
    //<< "\n Ga[1][0]: " 
    //<< Ga[1][0]
    //<< "\n Ga[1][1]: " 
    //<< Ga[1][1]
    //<< "\n Ga[1][2]: " 
    //<< Ga[1][2]
    //<< "\n Ga[2][0]: " 
    //<< Ga[2][0]
    //<< "\n Ga[2][1]: " 
    //<< Ga[2][1]
    //<< "\n Ga[2][2]: " 
    //<< Ga[2][2]
    //<< G4endl;

 //Matrix3f A;
 //A << M[0][0], M[0][1], M[0][2], M[1][0], M[1][1], M[1][2], M[2][0], M[2][1], M[2][2];

 //SelfAdjointEigenSolver <Matrix3f> es(A);
 //lambda = es.eigenvalues().[0];
 //lambdaVector = es.eigenvectors().col(0);
 //if (es.eigenvalues().[1]>lambda)   {
     //lambda = es.eigenvalues().[1];
     //lambdaVector = es.eigenvectors().col(1);
 //}
 //if (es.eigenvalues().[2]>lambda)   {
     //lambda = es.eigenvalues().[2];
     //lambdaVector = es.eigenvectors().col(2);
 //}

 //lambdaVector definition
 Ga11 = Ga[1][1];
 Ga00 = Ga[0][0];
 lambdaVectorZ = 1;
 if (Ga11==0) {
    lambdaVectorY = 0;
 }
 else {
    lambdaVectorY = -(Ga[1][2])*lambdaVectorZ/Ga11;
 }
 if (Ga00==0) {
    lambdaVectorX = 0;
 }
 else {
    lambdaVectorX = -(lambdaVectorY*Ga[0][1] + lambdaVectorZ*Ga[0][2])/Ga00;
 }
 
 //G4cout 
    //<< "\n lambdaVectorX: " 
    //<< lambdaVectorX
    //<< "\n lambdaVectorY: " 
    //<< lambdaVectorY
    //<< "\n lambdaVectorZ: " 
    //<< lambdaVectorZ
    //<< G4endl;

 //The fitted AJ line is Rav + t*lambdaVector, t real

 tA = (A.z() - Rav.z())/lambdaVectorZ;
 Af = G4ThreeVector((Rav.x() + tA*lambdaVectorX), (Rav.y() + tA*lambdaVectorY), A.z());
 AtAx = A.x() - Af.x();
 AtAy = A.y() - Af.y();
 if (AtAx == AtAx)  {
   analysisManager->FillH1(265, AtAx);
 }
 if (AtAy == AtAy)  {
   analysisManager->FillH1(266, AtAy);
 }

 Aztr = A.z() - poszBPIX12;
 A1 = G4ThreeVector(((A.x())*cos(YangleBPIX) + Aztr*sin(YangleBPIX)), (A.y()), (-(A.x())*sin(YangleBPIX) + Aztr*cos(YangleBPIX)));
 Aloc = G4ThreeVector((A1.x()), ((A1.y())*cos(XangleBPIX) - (A1.z())*sin(XangleBPIX)), ((A1.y())*sin(XangleBPIX) + (A1.z())*cos(XangleBPIX)));
 Afztr = Af.z() - poszBPIX12;
 Af1 = G4ThreeVector(((Af.x())*cos(YangleBPIX) + Afztr*sin(YangleBPIX)), (Af.y()), (-(Af.x())*sin(YangleBPIX) + Afztr*cos(YangleBPIX)));
 Afloc = G4ThreeVector((Af1.x()), ((Af1.y())*cos(XangleBPIX) - (Af1.z())*sin(XangleBPIX)), ((Af1.y())*sin(XangleBPIX) + (Af1.z())*cos(XangleBPIX)));
 AtAxloc = Aloc.x() - Afloc.x();
 AtAyloc = Aloc.y() - Afloc.y();
 if ((AtAxloc == AtAxloc) && (AtAxloc != 0)) {
   analysisManager->FillH1(285, AtAxloc);
 }
 if ((AtAyloc == AtAyloc) && (AtAyloc != 0))  {
   analysisManager->FillH1(286, AtAyloc);
 }

 tB = (B.z() - Rav.z())/lambdaVectorZ;
 Bf = G4ThreeVector((Rav.x() + tB*lambdaVectorX), (Rav.y() + tB*lambdaVectorY), B.z());
 BtBx = B.x() - Bf.x();
 BtBy = B.y() - Bf.y();
 if (BtBx == BtBx)  {
   analysisManager->FillH1(267, BtBx);
 }
 if (BtBy == BtBy)  {
   analysisManager->FillH1(268, BtBy);
 }

 Bztr = B.z() - poszBPIX34;
 B1 = G4ThreeVector(((B.x())*cos(YangleBPIX) + Bztr*sin(YangleBPIX)), (B.y()), (-(B.x())*sin(YangleBPIX) + Bztr*cos(YangleBPIX)));
 Bloc = G4ThreeVector((B1.x()), ((B1.y())*cos(XangleBPIX) - (B1.z())*sin(XangleBPIX)), ((B1.y())*sin(XangleBPIX) + (B1.z())*cos(XangleBPIX)));
 Bfztr = Bf.z() - poszBPIX34;
 Bf1 = G4ThreeVector(((Bf.x())*cos(YangleBPIX) + Bfztr*sin(YangleBPIX)), (Bf.y()), (-(Bf.x())*sin(YangleBPIX) + Bfztr*cos(YangleBPIX)));
 Bfloc = G4ThreeVector((Bf1.x()), ((Bf1.y())*cos(XangleBPIX) - (Bf1.z())*sin(XangleBPIX)), ((Bf1.y())*sin(XangleBPIX) + (Bf1.z())*cos(XangleBPIX)));
 BtBxloc = Bloc.x() - Bfloc.x();
 BtByloc = Bloc.y() - Bfloc.y();
 if ((BtBxloc == BtBxloc) && (BtBxloc != 0)) {
   analysisManager->FillH1(287, BtBxloc);
 }
 if ((BtByloc == BtByloc) && (BtByloc != 0))  {
   analysisManager->FillH1(288, BtByloc);
 }

 tC = (C.z() - Rav.z())/lambdaVectorZ;
 Cf = G4ThreeVector((Rav.x() + tC*lambdaVectorX), (Rav.y() + tC*lambdaVectorY), C.z());
 CtCx = C.x() - Cf.x();
 CtCy = C.y() - Cf.y();

 if (CtCx == CtCx)  {
   analysisManager->FillH1(269, CtCx);
 }
 if (CtCy == CtCy)  {
   analysisManager->FillH1(270, CtCy);
 }

 Cztr = C.z() - poszBPIX56;
 C1 = G4ThreeVector(((C.x())*cos(YangleBPIX) + Cztr*sin(YangleBPIX)), (C.y()), (-(C.x())*sin(YangleBPIX) + Cztr*cos(YangleBPIX)));
 Cloc = G4ThreeVector((C1.x()), ((C1.y())*cos(XangleBPIX) - (C1.z())*sin(XangleBPIX)), ((C1.y())*sin(XangleBPIX) + (C1.z())*cos(XangleBPIX)));
 Cfztr = Cf.z() - poszBPIX56;
 Cf1 = G4ThreeVector(((Cf.x())*cos(YangleBPIX) + Cfztr*sin(YangleBPIX)), (Cf.y()), (-(Cf.x())*sin(YangleBPIX) + Cfztr*cos(YangleBPIX)));
 Cfloc = G4ThreeVector((Cf1.x()), ((Cf1.y())*cos(XangleBPIX) - (Cf1.z())*sin(XangleBPIX)), ((Cf1.y())*sin(XangleBPIX) + (Cf1.z())*cos(XangleBPIX)));
 CtCxloc = Cloc.x() - Cfloc.x();
 CtCyloc = Cloc.y() - Cfloc.y();
 if ((CtCxloc == CtCxloc) && (CtCxloc != 0)) {
   analysisManager->FillH1(289, CtCxloc);
 }
 if ((CtCyloc == CtCyloc) && (CtCyloc != 0))  {
   analysisManager->FillH1(290, CtCyloc);
 }

 tD = (D.z() - Rav.z())/lambdaVectorZ;
 Df = G4ThreeVector((Rav.x() + tD*lambdaVectorX), (Rav.y() + tD*lambdaVectorY), D.z());
 DtDx = D.x() - Df.x();
 DtDy = D.y() - Df.y();
 if (DtDx == DtDx)  {
   analysisManager->FillH1(271, DtDx);
 }
 if (DtDy == DtDy)  {
   analysisManager->FillH1(272, DtDy);
 }

 Dztr = D.z() - poszBPIX78;
 D1 = G4ThreeVector(((D.x())*cos(YangleBPIX) + Dztr*sin(YangleBPIX)), (D.y()), (-(D.x())*sin(YangleBPIX) + Dztr*cos(YangleBPIX)));
 Dloc = G4ThreeVector((D1.x()), ((D1.y())*cos(XangleBPIX) - (D1.z())*sin(XangleBPIX)), ((D1.y())*sin(XangleBPIX) + (D1.z())*cos(XangleBPIX)));
 Dfztr = Df.z() - poszBPIX78;
 Df1 = G4ThreeVector(((Df.x())*cos(YangleBPIX) + Dfztr*sin(YangleBPIX)), (Df.y()), (-(Df.x())*sin(YangleBPIX) + Dfztr*cos(YangleBPIX)));
 Dfloc = G4ThreeVector((Df1.x()), ((Df1.y())*cos(XangleBPIX) - (Df1.z())*sin(XangleBPIX)), ((Df1.y())*sin(XangleBPIX) + (Df1.z())*cos(XangleBPIX)));
 DtDxloc = Dloc.x() - Dfloc.x();
 DtDyloc = Dloc.y() - Dfloc.y();
 if ((DtDxloc == DtDxloc) && (DtDxloc != 0)) {
   analysisManager->FillH1(291, DtDxloc);
 }
 if ((DtDyloc == DtDyloc) && (DtDyloc != 0))  {
   analysisManager->FillH1(292, DtDyloc);
 }

 tE = (E.z() - Rav.z())/lambdaVectorZ;
 Et = G4ThreeVector((Rav.x() + tE*lambdaVectorX), (Rav.y() + tE*lambdaVectorY), E.z());
 EtEx = E.x() - Et.x();
 EtEy = E.y() - Et.y();
 if (EtEx == EtEx)  {
   analysisManager->FillH1(119, EtEx);
 }
 if (EtEy == EtEy)  {
   analysisManager->FillH1(120, EtEy);
 }

 Eloc = G4ThreeVector((E.x()), ((E.y())*cos(XangleDUT) - (E.z())*sin(XangleDUT)), ((E.y())*sin(XangleDUT) + (E.z())*cos(XangleDUT)));
 Etloc = G4ThreeVector((Et.x()), ((Et.y())*cos(XangleDUT) - (Et.z())*sin(XangleDUT)), ((Et.y())*sin(XangleDUT) + (Et.z())*cos(XangleDUT)));
 EtExloc = Eloc.x() - Etloc.x();
 EtEyloc = Eloc.y() - Etloc.y();
 if ((EtExloc == EtExloc) && (EtExloc != 0)) {
   analysisManager->FillH1(281, EtExloc);
 }
 if ((EtEyloc == EtEyloc) && (EtEyloc != 0))  {
   analysisManager->FillH1(282, EtEyloc);
 }


 tF = (F.z() - Rav.z())/lambdaVectorZ;
 Ft = G4ThreeVector((Rav.x() + tF*lambdaVectorX), (Rav.y() + tF*lambdaVectorY), F.z());
 FtFx = F.x() - Ft.x();
 FtFy = F.y() - Ft.y();
 if (FtFx == FtFx)  {
   analysisManager->FillH1(121, FtFx);
 }
 if (FtFy == FtFy)  {
   analysisManager->FillH1(122, FtFy);
 }

 Floc = G4ThreeVector((F.x()), ((F.y())*cos(XangleDUT) - (F.z())*sin(XangleDUT)), ((F.y())*sin(XangleDUT) + (F.z())*cos(XangleDUT)));
 Ftloc = G4ThreeVector((Ft.x()), ((Ft.y())*cos(XangleDUT) - (Ft.z())*sin(XangleDUT)), ((Ft.y())*sin(XangleDUT) + (Ft.z())*cos(XangleDUT)));
 FtFxloc = Floc.x() - Ftloc.x();
 FtFyloc = Floc.y() - Ftloc.y();
 if ((FtFxloc == FtFxloc) && (FtFxloc != 0)) {
   analysisManager->FillH1(283, FtFxloc);
 }
 if ((FtFyloc == FtFyloc) && (FtFyloc != 0))  {
   analysisManager->FillH1(284, FtFyloc);
 }

 tG = (G.z() - Rav.z())/lambdaVectorZ;
 Gf = G4ThreeVector((Rav.x() + tG*lambdaVectorX), (Rav.y() + tG*lambdaVectorY), G.z());
 GtGx = G.x() - Gf.x();
 GtGy = G.y() - Gf.y();

 if (GtGx == GtGx)  {
   analysisManager->FillH1(273, GtGx);
 }
 if (GtGy == GtGy)  {
   analysisManager->FillH1(274, GtGy);
 }

 Gztr = G.z() - poszBPIX910;
 G1 = G4ThreeVector(((G.x())*cos(YangleBPIX) + Gztr*sin(YangleBPIX)), (G.y()), (-(G.x())*sin(YangleBPIX) + Gztr*cos(YangleBPIX)));
 Gloc = G4ThreeVector((G1.x()), ((G1.y())*cos(XangleBPIX) - (G1.z())*sin(XangleBPIX)), ((G1.y())*sin(XangleBPIX) + (G1.z())*cos(XangleBPIX)));
 Gfztr = Gf.z() - poszBPIX910;
 Gf1 = G4ThreeVector(((Gf.x())*cos(YangleBPIX) + Gfztr*sin(YangleBPIX)), (Gf.y()), (-(Gf.x())*sin(YangleBPIX) + Gfztr*cos(YangleBPIX)));
 Gfloc = G4ThreeVector((Gf1.x()), ((Gf1.y())*cos(XangleBPIX) - (Gf1.z())*sin(XangleBPIX)), ((Gf1.y())*sin(XangleBPIX) + (Gf1.z())*cos(XangleBPIX)));
 GtGxloc = Gloc.x() - Gfloc.x();
 GtGyloc = Gloc.y() - Gfloc.y();
 if ((GtGxloc == GtGxloc) && (GtGxloc != 0)) {
   analysisManager->FillH1(293, GtGxloc);
 }
 if ((GtGyloc == GtGyloc) && (GtGyloc != 0))  {
   analysisManager->FillH1(294, GtGyloc);
 }

 tH = (H.z() - Rav.z())/lambdaVectorZ;
 Hf = G4ThreeVector((Rav.x() + tH*lambdaVectorX), (Rav.y() + tH*lambdaVectorY), H.z());
 HtHx = H.x() - Hf.x();
 HtHy = H.y() - Hf.y();
 if (HtHx == HtHx)  {
   analysisManager->FillH1(275, HtHx);
 }
 if (HtHy == HtHy)  {
   analysisManager->FillH1(276, HtHy);
 }

 Hztr = H.z() - poszBPIX1112;
 H1 = G4ThreeVector(((H.x())*cos(YangleBPIX) + Hztr*sin(YangleBPIX)), (H.y()), (-(H.x())*sin(YangleBPIX) + Hztr*cos(YangleBPIX)));
 Hloc = G4ThreeVector((H1.x()), ((H1.y())*cos(XangleBPIX) - (H1.z())*sin(XangleBPIX)), ((H1.y())*sin(XangleBPIX) + (H1.z())*cos(XangleBPIX)));
 Hfztr = Hf.z() - poszBPIX1112;
 Hf1 = G4ThreeVector(((Hf.x())*cos(YangleBPIX) + Hfztr*sin(YangleBPIX)), (Hf.y()), (-(Hf.x())*sin(YangleBPIX) + Hfztr*cos(YangleBPIX)));
 Hfloc = G4ThreeVector((Hf1.x()), ((Hf1.y())*cos(XangleBPIX) - (Hf1.z())*sin(XangleBPIX)), ((Hf1.y())*sin(XangleBPIX) + (Hf1.z())*cos(XangleBPIX)));
 HtHxloc = Hloc.x() - Hfloc.x();
 HtHyloc = Hloc.y() - Hfloc.y();
 if ((HtHxloc == HtHxloc) && (HtHxloc != 0)) {
   analysisManager->FillH1(295, HtHxloc);
 }
 if ((HtHyloc == HtHyloc) && (HtHyloc != 0))  {
   analysisManager->FillH1(296, HtHyloc);
 }

 tI = (I.z() - Rav.z())/lambdaVectorZ;
 If = G4ThreeVector((Rav.x() + tI*lambdaVectorX), (Rav.y() + tI*lambdaVectorY), I.z());
 ItIx = I.x() - If.x();
 ItIy = I.y() - If.y();
 if (ItIx == ItIx)  {
   analysisManager->FillH1(277, ItIx);
 }
 if (ItIy == ItIy)  {
   analysisManager->FillH1(278, ItIy);
 }

 Iztr = I.z() - poszBPIX1314;
 I1 = G4ThreeVector(((I.x())*cos(YangleBPIX) + Iztr*sin(YangleBPIX)), (I.y()), (-(I.x())*sin(YangleBPIX) + Iztr*cos(YangleBPIX)));
 Iloc = G4ThreeVector((I1.x()), ((I1.y())*cos(XangleBPIX) - (I1.z())*sin(XangleBPIX)), ((I1.y())*sin(XangleBPIX) + (I1.z())*cos(XangleBPIX)));
 Ifztr = If.z() - poszBPIX1314;
 If1 = G4ThreeVector(((If.x())*cos(YangleBPIX) + Ifztr*sin(YangleBPIX)), (If.y()), (-(If.x())*sin(YangleBPIX) + Ifztr*cos(YangleBPIX)));
 Ifloc = G4ThreeVector((If1.x()), ((If1.y())*cos(XangleBPIX) - (If1.z())*sin(XangleBPIX)), ((If1.y())*sin(XangleBPIX) + (If1.z())*cos(XangleBPIX)));
 ItIxloc = Iloc.x() - Ifloc.x();
 ItIyloc = Iloc.y() - Ifloc.y();
 if ((ItIxloc == ItIxloc) && (ItIxloc != 0)) {
   analysisManager->FillH1(297, ItIxloc);
 }
 if ((ItIyloc == ItIyloc) && (ItIyloc != 0))  {
   analysisManager->FillH1(298, ItIyloc);
 }

 tJ = (J.z() - Rav.z())/lambdaVectorZ;
 Jf = G4ThreeVector((Rav.x() + tJ*lambdaVectorX), (Rav.y() + tJ*lambdaVectorY), J.z());
 JtJx = J.x() - Jf.x();
 JtJy = J.y() - Jf.y();
 if (JtJx == JtJx)  {
   analysisManager->FillH1(279, JtJx);
 }
 if (JtJy == JtJy)  {
   analysisManager->FillH1(280, JtJy);
 }

 Jztr = J.z() - poszBPIX1516;
 J1 = G4ThreeVector(((J.x())*cos(YangleBPIX) + Jztr*sin(YangleBPIX)), (J.y()), (-(J.x())*sin(YangleBPIX) + Jztr*cos(YangleBPIX)));
 Jloc = G4ThreeVector((J1.x()), ((J1.y())*cos(XangleBPIX) - (J1.z())*sin(XangleBPIX)), ((J1.y())*sin(XangleBPIX) + (J1.z())*cos(XangleBPIX)));
 Jfztr = Jf.z() - poszBPIX1516;
 Jf1 = G4ThreeVector(((Jf.x())*cos(YangleBPIX) + Jfztr*sin(YangleBPIX)), (Jf.y()), (-(Jf.x())*sin(YangleBPIX) + Jfztr*cos(YangleBPIX)));
 Jfloc = G4ThreeVector((Jf1.x()), ((Jf1.y())*cos(XangleBPIX) - (Jf1.z())*sin(XangleBPIX)), ((Jf1.y())*sin(XangleBPIX) + (Jf1.z())*cos(XangleBPIX)));
 JtJxloc = Jloc.x() - Jfloc.x();
 JtJyloc = Jloc.y() - Jfloc.y();
 if ((JtJxloc == JtJxloc) && (JtJxloc != 0)) {
   analysisManager->FillH1(299, JtJxloc);
 }
 if ((JtJyloc == JtJyloc) && (JtJyloc != 0))  {
   analysisManager->FillH1(300, JtJyloc);
 }

 fDiff1x = E.x() - EReal.x();
 fDiff1y = E.y() - EReal.y();
 fDiff2x = F.x() - FReal.x();
 fDiff2y = F.y() - FReal.y();

 if (fHitSensor1 == 1)   {
    /*G4cout
	<< "\n E.x(): " 
        << E.x()
	<< "\n EReal.x(): " 
        << EReal.x()
	<< "\n E.y(): " 
        << E.y()
	<< "\n EReal.y(): " 
        << EReal.y()
	<< "\n E.z(): " 
        << E.z()
	<< "\n EReal.z(): " 
        << EReal.z()
	<< G4endl; */
    analysisManager->FillH1(191, E.x());
    analysisManager->FillH1(192, E.y());
    analysisManager->FillH1(193, E.z());
    analysisManager->FillH1(213, fDiff1x);
    analysisManager->FillH1(214, fDiff1y);

    analysisManager->FillH2(9, E.x(), E.y());
    analysisManager->FillH2(19, -Eloc.x(), Eloc.y());
 }
 if (fHitSensor2 == 1)   {
    /*G4cout
	<< "\n F.x(): " 
        << F.x()
	<< "\n FReal.x(): " 
        << FReal.x()
	<< "\n F.y(): " 
        << F.y()
	<< "\n FReal.y(): " 
        << FReal.y()
	<< "\n F.z(): " 
        << F.z()
	<< "\n FReal.z(): " 
        << FReal.z()
	<< G4endl; */
    analysisManager->FillH1(194, F.x());
    analysisManager->FillH1(195, F.y());
    analysisManager->FillH1(196, F.z());
    analysisManager->FillH1(215, fDiff2x);
    analysisManager->FillH1(216, fDiff2y);

    analysisManager->FillH2(10, F.x(), F.y());
    analysisManager->FillH2(20, -Floc.x(), Floc.y());
 }


    /*analysisManager->FillH1(217, A.x());
    analysisManager->FillH1(218, A.y());
    analysisManager->FillH1(219, A.z());

    analysisManager->FillH1(223, B.x());
    analysisManager->FillH1(224, B.y());
    analysisManager->FillH1(225, B.z());

    analysisManager->FillH1(229, C.x());
    analysisManager->FillH1(230, C.y());
    analysisManager->FillH1(231, C.z());

    analysisManager->FillH1(235, D.x());
    analysisManager->FillH1(236, D.y());
    analysisManager->FillH1(237, D.z());

    analysisManager->FillH1(241, G.x());
    analysisManager->FillH1(242, G.y());
    analysisManager->FillH1(243, G.z());

    analysisManager->FillH1(247, H.x());
    analysisManager->FillH1(248, H.y());
    analysisManager->FillH1(249, H.z());

    analysisManager->FillH1(253, I.x());
    analysisManager->FillH1(254, I.y());
    analysisManager->FillH1(255, I.z());

    analysisManager->FillH1(259, J.x());
    analysisManager->FillH1(260, J.y());
    analysisManager->FillH1(261, J.z());*/

 if (fHitPixelDet[1] == 1)   {
    /*G4cout
	<< "\n A.x(): " 
        << A.x()
	<< "\n AReal.x(): " 
        << AReal.x()
	<< "\n A.y(): " 
        << A.y()
	<< "\n AReal.y(): " 
        << AReal.y()
	<< "\n A.z(): " 
        << A.z()
	<< "\n AReal.z(): " 
        << AReal.z()
	<< G4endl; */
    analysisManager->FillH1(217, A.x());
    analysisManager->FillH1(218, A.y());
    analysisManager->FillH1(219, A.z());

    analysisManager->FillH2(1, A.x(), A.y());
    analysisManager->FillH2(11, -Aloc.x(), Aloc.y());
 }
 if (fHitPixelDet[2] == 1)   {
    /*G4cout
	<< "\n A.x(): " 
        << A.x()
	<< "\n AReal.x(): " 
        << AReal.x()
	<< "\n A.y(): " 
        << A.y()
	<< "\n AReal.y(): " 
        << AReal.y()
	<< "\n A.z(): " 
        << A.z()
	<< "\n AReal.z(): " 
        << AReal.z()
	<< G4endl; */
    analysisManager->FillH1(220, A.x());
    analysisManager->FillH1(221, A.y());
    analysisManager->FillH1(222, A.z()); 

    analysisManager->FillH2(1, A.x(), A.y());
    analysisManager->FillH2(11, -Aloc.x(), Aloc.y());
 }
 if (fHitPixelDet[3] == 1)   {
    /*G4cout
	<< "\n B.x(): " 
        << B.x()
	<< "\n BReal.x(): " 
        << BReal.x()
	<< "\n B.y(): " 
        << B.y()
	<< "\n BReal.y(): " 
        << BReal.y()
	<< "\n B.z(): " 
        << B.z()
	<< "\n BReal.z(): " 
        << BReal.z()
	<< G4endl; */
    analysisManager->FillH1(223, B.x());
    analysisManager->FillH1(224, B.y());
    analysisManager->FillH1(225, B.z()); 

    analysisManager->FillH2(2, B.x(), B.y());
    analysisManager->FillH2(12, -Bloc.x(), Bloc.y());
 }
 if (fHitPixelDet[4] == 1)   {
    /*G4cout
	<< "\n B.x(): " 
        << B.x()
	<< "\n BReal.x(): " 
        << BReal.x()
	<< "\n B.y(): " 
        << B.y()
	<< "\n BReal.y(): " 
        << BReal.y()
	<< "\n B.z(): " 
        << B.z()
	<< "\n BReal.z(): " 
        << BReal.z()
	<< G4endl; */
    analysisManager->FillH1(226, B.x());
    analysisManager->FillH1(227, B.y());
    analysisManager->FillH1(228, B.z()); 

    analysisManager->FillH2(2, B.x(), B.y());
    analysisManager->FillH2(12, -Bloc.x(), Bloc.y());
 }
 if (fHitPixelDet[5] == 1)   {
    /*G4cout
	<< "\n C.x(): " 
        << C.x()
	<< "\n CReal.x(): " 
        << CReal.x()
	<< "\n C.y(): " 
        << C.y()
	<< "\n CReal.y(): " 
        << CReal.y()
	<< "\n C.z(): " 
        << C.z()
	<< "\n CReal.z(): " 
        << CReal.z()
	<< G4endl; */
    analysisManager->FillH1(229, C.x());
    analysisManager->FillH1(230, C.y());
    analysisManager->FillH1(231, C.z()); 

    analysisManager->FillH2(3, C.x(), C.y());
    analysisManager->FillH2(13, -Cloc.x(), Cloc.y());
 }
 if (fHitPixelDet[6] == 1)   {
    /*G4cout
	<< "\n C.x(): " 
        << C.x()
	<< "\n CReal.x(): " 
        << CReal.x()
	<< "\n C.y(): " 
        << C.y()
	<< "\n CReal.y(): " 
        << CReal.y()
	<< "\n C.z(): " 
        << C.z()
	<< "\n CReal.z(): " 
        << CReal.z()
	<< G4endl; */
    analysisManager->FillH1(232, C.x());
    analysisManager->FillH1(233, C.y());
    analysisManager->FillH1(234, C.z()); 

    analysisManager->FillH2(3, C.x(), C.y());
    analysisManager->FillH2(13, -Cloc.x(), Cloc.y());
 }
 if (fHitPixelDet[7] == 1)   {
    /*G4cout
	<< "\n D.x(): " 
        << D.x()
	<< "\n DReal.x(): " 
        << DReal.x()
	<< "\n D.y(): " 
        << D.y()
	<< "\n DReal.y(): " 
        << DReal.y()
	<< "\n D.z(): " 
        << D.z()
	<< "\n DReal.z(): " 
        << DReal.z()
	<< G4endl; */
    analysisManager->FillH1(235, D.x());
    analysisManager->FillH1(236, D.y());
    analysisManager->FillH1(237, D.z()); 

    analysisManager->FillH2(4, D.x(), D.y());
    analysisManager->FillH2(14, -Dloc.x(), Dloc.y());
 }
 if (fHitPixelDet[8] == 1)   {
    /*G4cout
	<< "\n D.x(): " 
        << D.x()
	<< "\n DReal.x(): " 
        << DReal.x()
	<< "\n D.y(): " 
        << D.y()
	<< "\n DReal.y(): " 
        << DReal.y()
	<< "\n D.z(): " 
        << D.z()
	<< "\n DReal.z(): " 
        << DReal.z()
	<< G4endl; */
    analysisManager->FillH1(238, D.x());
    analysisManager->FillH1(239, D.y());
    analysisManager->FillH1(240, D.z()); 

    analysisManager->FillH2(4, D.x(), D.y());
    analysisManager->FillH2(14, -Dloc.x(), Dloc.y());
 }

 if (fHitPixelDet[9] == 1)   {
    /*G4cout
	<< "\n G.x(): " 
        << G.x()
	<< "\n GReal.x(): " 
        << GReal.x()
	<< "\n G.y(): " 
        << G.y()
	<< "\n GReal.y(): " 
        << GReal.y()
	<< "\n G.z(): " 
        << G.z()
	<< "\n GReal.z(): " 
        << GReal.z()
	<< G4endl; */
    analysisManager->FillH1(241, G.x());
    analysisManager->FillH1(242, G.y());
    analysisManager->FillH1(243, G.z()); 

    analysisManager->FillH2(5, G.x(), G.y());
    analysisManager->FillH2(15, -Gloc.x(), Gloc.y());
 }
 if (fHitPixelDet[10] == 1)   {
    /*G4cout
	<< "\n G.x(): " 
        << G.x()
	<< "\n GReal.x(): " 
        << GReal.x()
	<< "\n G.y(): " 
        << G.y()
	<< "\n GReal.y(): " 
        << GReal.y()
	<< "\n G.z(): " 
        << G.z()
	<< "\n GReal.z(): " 
        << GReal.z()
	<< G4endl; */
    analysisManager->FillH1(244, G.x());
    analysisManager->FillH1(245, G.y());
    analysisManager->FillH1(246, G.z()); 

    analysisManager->FillH2(5, G.x(), G.y());
    analysisManager->FillH2(15, -Gloc.x(), Gloc.y());
 }
 if (fHitPixelDet[11] == 1)   {
    /*G4cout
	<< "\n H.x(): " 
        << H.x()
	<< "\n HReal.x(): " 
        << HReal.x()
	<< "\n H.y(): " 
        << H.y()
	<< "\n HReal.y(): " 
        << HReal.y()
	<< "\n H.z(): " 
        << H.z()
	<< "\n HReal.z(): " 
        << HReal.z()
	<< G4endl; */
    analysisManager->FillH1(247, H.x());
    analysisManager->FillH1(248, H.y());
    analysisManager->FillH1(249, H.z()); 

    analysisManager->FillH2(6, H.x(), H.y());
    analysisManager->FillH2(16, -Hloc.x(), Hloc.y());
 }
 if (fHitPixelDet[12] == 1)   {
    /*G4cout
	<< "\n H.x(): " 
        << H.x()
	<< "\n HReal.x(): " 
        << HReal.x()
	<< "\n H.y(): " 
        << H.y()
	<< "\n HReal.y(): " 
        << HReal.y()
	<< "\n H.z(): " 
        << H.z()
	<< "\n HReal.z(): " 
        << HReal.z()
	<< G4endl; */
    analysisManager->FillH1(250, H.x());
    analysisManager->FillH1(251, H.y());
    analysisManager->FillH1(252, H.z()); 

    analysisManager->FillH2(6, H.x(), H.y());
    analysisManager->FillH2(16, -Hloc.x(), Hloc.y());
 }
 if (fHitPixelDet[13] == 1)   {
    /*G4cout
	<< "\n I.x(): " 
        << I.x()
	<< "\n IReal.x(): " 
        << IReal.x()
	<< "\n I.y(): " 
        << I.y()
	<< "\n IReal.y(): " 
        << IReal.y()
	<< "\n I.z(): " 
        << I.z()
	<< "\n IReal.z(): " 
        << IReal.z()
	<< G4endl; */
    analysisManager->FillH1(253, I.x());
    analysisManager->FillH1(254, I.y());
    analysisManager->FillH1(255, I.z()); 

    analysisManager->FillH2(7, I.x(), I.y());
    analysisManager->FillH2(17, -Iloc.x(), Iloc.y());
 }
 if (fHitPixelDet[14] == 1)   {
    /*G4cout
	<< "\n I.x(): " 
        << I.x()
	<< "\n IReal.x(): " 
        << IReal.x()
	<< "\n I.y(): " 
        << I.y()
	<< "\n IReal.y(): " 
        << IReal.y()
	<< "\n I.z(): " 
        << I.z()
	<< "\n IReal.z(): " 
        << IReal.z()
	<< G4endl; */
    analysisManager->FillH1(256, I.x());
    analysisManager->FillH1(257, I.y());
    analysisManager->FillH1(258, I.z()); 

    analysisManager->FillH2(7, I.x(), I.y());
    analysisManager->FillH2(17, -Iloc.x(), Iloc.y());
 }
 if (fHitPixelDet[15] == 1)   {
    /*G4cout
	<< "\n J.x(): " 
        << J.x()
	<< "\n JReal.x(): " 
        << JReal.x()
	<< "\n J.y(): " 
        << J.y()
	<< "\n JReal.y(): " 
        << JReal.y()
	<< "\n J.z(): " 
        << J.z()
	<< "\n JReal.z(): " 
        << JReal.z()
	<< G4endl; */
    analysisManager->FillH1(259, J.x());
    analysisManager->FillH1(260, J.y());
    analysisManager->FillH1(261, J.z());

    analysisManager->FillH2(8, J.x(), J.y());
    analysisManager->FillH2(18, -Jloc.x(), Jloc.y()); 
 }
 if (fHitPixelDet[16] == 1)   {
    /*G4cout
	<< "\n J.x(): " 
        << J.x()
	<< "\n JReal.x(): " 
        << JReal.x()
	<< "\n J.y(): " 
        << J.y()
	<< "\n JReal.y(): " 
        << JReal.y()
	<< "\n J.z(): " 
        << J.z()
	<< "\n JReal.z(): " 
        << JReal.z()
	<< G4endl; */
    analysisManager->FillH1(262, J.x());
    analysisManager->FillH1(263, J.y());
    analysisManager->FillH1(264, J.z()); 

    analysisManager->FillH2(8, J.x(), J.y());
    analysisManager->FillH2(18, -Jloc.x(), Jloc.y());
 }


 // fill ntuple  
 //for (int i = 0; i < 10; i = i + 1) {
   //analysisManager->FillNtupleDColumn(i, fEnergyStrip1[i]);
   //analysisManager->FillNtupleDColumn(i+1016, fEnergyStrip2[i]);
 //}
 //analysisManager->AddNtupleRow(); 

 //G4cout 
     //<< "\n Last track length of secondary created in detector calculated = " 
     //<< G4BestUnit(fTrack2, "Length")
     //<< "\n End of Event. Secondary track length from detector secondaries = " 
     //<< G4BestUnit(fSecondaryDetTrackLength, "Length")
     //<< G4endl;

 //Visualize event if there is a track longer than 1 cm
 //if (fTrack2 > 1.0*cm)  {
     //G4cout 
         //<< "\n fTrack2 = " 
         //<< G4BestUnit(fTrack2, "Length")
         //<< G4endl;
     //G4EventManager* evMan = G4EventManager::GetEventManager();
     //evMan->KeepTheCurrentEvent();
 //}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

