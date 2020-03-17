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
// $Id: SteppingAction.cc 67268 2013-02-13 11:38:40Z ihrivnac $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"
#include "StackingAction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4VTrajectory.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4Positron.hh"
#include "G4Electron.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4Gamma.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <cmath>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* DA, RunAction* RA, EventAction* EA, StackingAction* SA)
:G4UserSteppingAction(), fDetectorconstruction(DA), fRunaction(RA), fEventaction(EA), fStackingaction(SA)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
 //get volume of the current step
  G4VPhysicalVolume* volume 
  = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  //step size
  G4double stepSize = step->GetStepLength();  

  //collect energy step by step
  G4double edep = step->GetTotalEnergyDeposit();

  //track position coordinates
  G4StepPoint* pre = step->GetPreStepPoint();
  G4StepPoint* post = step->GetPostStepPoint();

  G4double xp_after,yp_after,zp_after,x_after,y_after,z_after,x,y,z,xPix1,yPix1,zPix1,xPix,yPix,zPix,zPixTr,xPixM;
  //G4double xp, yp, zp;
  xp_after = pre->GetPosition().x();
  yp_after = pre->GetPosition().y();
  zp_after = pre->GetPosition().z();
  x_after = post->GetPosition().x();
  y_after = post->GetPosition().y();
  z_after = post->GetPosition().z();

  //geometry parameters
  G4double Det1SizeZ, Det2SizeZ, Dist, Strip1Depth, Strip1Length, Strip2Depth, Strip2Length, StripDist, StripWidth, BPIXSizeX, BPIXSizeY, BPIXSizeZ, pixelDepth, pixelX, pixelY, ROChor, ROCvert,  XangleDUT, XangleBPIX, YangleBPIX, ElField1, ElField2, poszBPIX12, poszBPIX34, poszBPIX56, poszBPIX78, poszBPIX910, poszBPIX1112, poszBPIX1314, poszBPIX1516;
  //G4double posEndArm1, posBeginningArm2;
  G4int div, div1, div2, divp1, divp2, dirx, diry, dirxM, mod, ROC, row, col;
  G4int totalNbStrips, stripNo;
  Det1SizeZ=fDetectorconstruction->GetSize1();
  Det2SizeZ=fDetectorconstruction->GetSize2();
  Dist=fDetectorconstruction->GetDist();
  Strip1Depth=fDetectorconstruction->GetStrip1Depth();
  Strip1Length=fDetectorconstruction->GetStrip1Length();
  Strip2Depth=fDetectorconstruction->GetStrip2Depth();
  Strip2Length=fDetectorconstruction->GetStrip2Length();
  StripDist=fDetectorconstruction->GetStripDist();
  StripWidth=fDetectorconstruction->GetStripWidth();
  //posEndArm1=-(fDetectorconstruction->Getpos_EndArm1Abs());
  //posBeginningArm2=fDetectorconstruction->Getpos_BeginningArm2Abs();
  poszBPIX12=fDetectorconstruction->GetZPosBPIX12();
  poszBPIX34=fDetectorconstruction->GetZPosBPIX34();
  poszBPIX56=fDetectorconstruction->GetZPosBPIX56();
  poszBPIX78=fDetectorconstruction->GetZPosBPIX78();
  poszBPIX910=fDetectorconstruction->GetZPosBPIX910();
  poszBPIX1112=fDetectorconstruction->GetZPosBPIX1112();
  poszBPIX1314=fDetectorconstruction->GetZPosBPIX1314();
  poszBPIX1516=fDetectorconstruction->GetZPosBPIX1516();
  BPIXSizeX=fDetectorconstruction->GetSizeXBPIX();
  BPIXSizeY=fDetectorconstruction->GetSizeYBPIX();
  BPIXSizeZ=fDetectorconstruction->GetSizeZBPIX();
  pixelX=fDetectorconstruction->GetPixelPitchX();
  pixelY=fDetectorconstruction->GetPixelPitchY();
  pixelDepth=fDetectorconstruction->GetPixelDepth();
  XangleDUT=fDetectorconstruction->GetDUTangleX();
  XangleBPIX=fDetectorconstruction->GetBPIXangleX();
  YangleBPIX=fDetectorconstruction->GetBPIXangleY();
  //ElField1=fDetectorconstruction->GetElField1();
  //ElField2=fDetectorconstruction->GetElField2();
  totalNbStrips = 1016;
  ROChor = 80*pixelX + pixelX;
  ROCvert = 52*pixelY + 2*pixelY;

  //momentum direction when entering BPIX12 and exiting BPIX1516
  G4double momDir1x, momDir1y, momDir1z, momDir2x, momDir2y, momDir2z;

  //position if the DUT wasn't rotated
  //xp = xp_after;
  //yp = yp_after*cos(XangleDUT) - zp_after*sin(XangleDUT);
  //zp = yp_after*sin(XangleDUT) + zp_after*cos(XangleDUT);
  x = x_after;
  y = y_after*cos(XangleDUT) - z_after*sin(XangleDUT);
  z = y_after*sin(XangleDUT) + z_after*cos(XangleDUT);

  //position if the BPIX modules weren't rotated
  zPixTr = 0.0;
  xPix1 = 0.0;
  yPix1 = 0.0;
  zPix1 = 0.0;
  xPix = 0.0;
  yPix = 0.0;
  zPix = 0.0;

  //secret number
  G4int iSecret, jSecret;  

  //get track length, track ID and track's vertex position of the current step
  G4Track* track = step->GetTrack();
  G4double length = track->GetTrackLength();
  G4ThreeVector primVert = track->GetVertexPosition();
  //G4int trackID = track->GetTrackID();

  //first, middle and last point of primary inside each detector
  G4ThreeVector A1, A, A2, B1, B, B2, C1, C, C2, D1, D, D2, E1, E, E2, F1, F, F2, G1, G, G2, H1, H, H2, I1, I, I2, J1, J, J2;

  //particle definition
  const G4ParticleDefinition* particleDefinition = track->GetParticleDefinition();

  //track length of primary particle
  if (track->GetTrackID() == 1)  {
    fRunaction->AddTrackLength(stepSize);
    G4AnalysisManager::Instance()->FillH1(17,stepSize);

    fEventaction->AddPrimaryTrackLength(stepSize);

    if (volume == fDetectorconstruction->GetBPIX12())   {
	if (pre->GetStepStatus() == fGeomBoundary)   {
	    A1 = pre->GetPosition();
	    momDir1x = track->GetMomentumDirection().x();
	    momDir1y = track->GetMomentumDirection().y();
	    momDir1z = track->GetMomentumDirection().z();

            fEventaction->AddMomentumDirection1(momDir1x, momDir1y, momDir1z);
    	    fEventaction->AddAReal(A1);

            /*G4cout 
               << "\n Primary has entered BPIX 12 at: " 
               << G4BestUnit(A1, "Length")
               //<< "\n with momentum direction: " 
               //<< momDir1x
               //<< " "
               //<< momDir1y
               //<< " "
               //<< momDir1z
               << G4endl;*/
        }
	if (post->GetStepStatus() == fGeomBoundary)   {
	    A2 = post->GetPosition();

            /*G4cout 
               << "\n Primary has exited BPIX 12 at: " angle
               << G4BestUnit(A2, "Length")
               << G4endl; */
        }
    }

    if (volume == fDetectorconstruction->GetBPIX34())   {
	if (pre->GetStepStatus() == fGeomBoundary)   {
	    B1 = pre->GetPosition();

    	    fEventaction->AddBReal(B1);

            /*G4cout 
               << "\n Primary has entered BPIX 34 at: " 
               << G4BestUnit(B1, "Length")
               << G4endl; */
        }
	if (post->GetStepStatus() == fGeomBoundary)   {
	    B2 = post->GetPosition();

            /*G4cout 
               << "\n Primary has exited BPIX 34 at: " 
               << G4BestUnit(B2, "Length")
               << G4endl; */
        }
    }

    if (volume == fDetectorconstruction->GetBPIX56())   {
	if (pre->GetStepStatus() == fGeomBoundary)   {
	    C1 = pre->GetPosition();

    	    fEventaction->AddCReal(C1);

            /*G4cout 
               << "\n Primary has entered BPIX 56 at: " 
               << G4BestUnit(C1, "Length")
               << G4endl; */
        }
	if (post->GetStepStatus() == fGeomBoundary)   {
	    C2 = post->GetPosition();


            /*G4cout 
               << "\n Primary has exited BPIX 56 at: " 
               << G4BestUnit(C2, "Length")
               << G4endl; */
        }
    }

    if (volume == fDetectorconstruction->GetBPIX78())   {
	if (pre->GetStepStatus() == fGeomBoundary)   {
	    D1 = pre->GetPosition();

    	    fEventaction->AddDReal(D1);

            /*G4cout 
               << "\n Primary has entered BPIX 78 at: " 
               << G4BestUnit(C1, "Length")
               << G4endl; */
        }
	if (post->GetStepStatus() == fGeomBoundary)   {
	    D2 = post->GetPosition();

            /*G4cout 
               << "\n Primary has exited BPIX 78 at: " 
               << G4BestUnit(D2, "Length")
               << G4endl; */
        }
    }

    if (volume == fDetectorconstruction->GetDet1())   {
	if (pre->GetStepStatus() == fGeomBoundary)   {
	    E1 = pre->GetPosition();

    	    fEventaction->AddEReal(E1);

            /*G4cout 
               << "\n Primary has entered Strip Sensor 1 at: " 
               << G4BestUnit(E1, "Length")
               << G4endl; */
        }
	if (post->GetStepStatus() == fGeomBoundary)   {
	    E2 = post->GetPosition();

            /*G4cout 
               << "\n Primary has exited Strip Sensor 1 at: " 
               << G4BestUnit(E2, "Length")
               << G4endl; */
        }
    }

    if (volume == fDetectorconstruction->GetDet2())   {
	if (pre->GetStepStatus() == fGeomBoundary)   {
	    F1 = pre->GetPosition();

    	    fEventaction->AddFReal(F1);

            /*G4cout 
               << "\n Primary has entered Strip Sensor 2 at: " 
               << G4BestUnit(F1, "Length")
               << G4endl; */
        }
	if (post->GetStepStatus() == fGeomBoundary)   {
	    F2 = post->GetPosition();

            /*G4cout 
               << "\n Primary has exited Strip Sensor 2 at: " 
               << G4BestUnit(F2, "Length")
               << G4endl; */
        }
    }

    if (volume == fDetectorconstruction->GetBPIX910())   {
	if (pre->GetStepStatus() == fGeomBoundary)   {
	    G1 = pre->GetPosition();

    	    fEventaction->AddGReal(G1);

            /*G4cout 
               << "\n Primary has entered BPIX 910 at: " 
               << G4BestUnit(G1, "Length")
               << G4endl; */
        }
	if (post->GetStepStatus() == fGeomBoundary)   {
	    G2 = post->GetPosition();

            /*G4cout 
               << "\n Primary has exited BPIX 910 at: " 
               << G4BestUnit(G2, "Length")
               << G4endl; */
        }
    }

    if (volume == fDetectorconstruction->GetBPIX1112())   {
	if (pre->GetStepStatus() == fGeomBoundary)   {
	    H1 = pre->GetPosition();

    	    fEventaction->AddHReal(H1);

            /*G4cout 
               << "\n Primary has entered BPIX 1112 at: " 
               << G4BestUnit(H1, "Length")
               << G4endl; */
        }
	if (post->GetStepStatus() == fGeomBoundary)   {
	    H2 = post->GetPosition();

            /*G4cout 
               << "\n Primary has exited BPIX 1112 at: " 
               << G4BestUnit(H2, "Length")
               << G4endl; */
        }
    }

    if (volume == fDetectorconstruction->GetBPIX1314())   {
	if (pre->GetStepStatus() == fGeomBoundary)   {
	    I1 = pre->GetPosition();

    	    fEventaction->AddIReal(I1);

            /*G4cout 
               << "\n Primary has entered BPIX 1314 at: " 
               << G4BestUnit(I1, "Length")
               << G4endl; */
        }
	if (post->GetStepStatus() == fGeomBoundary)   {
	    I2 = post->GetPosition();

            /*G4cout 
               << "\n Primary has exited BPIX 1314 at: " 
               << G4BestUnit(I2, "Length")
               << G4endl; */
        }
    }

    if (volume == fDetectorconstruction->GetBPIX1516())   {
	if (pre->GetStepStatus() == fGeomBoundary)   {
	    J1 = pre->GetPosition();

    	    fEventaction->AddJReal(J1);

            /*G4cout 
               << "\n Primary has entered BPIX1516 at: " 
               << G4BestUnit(J1, "Length")
               << G4endl; */
        }
	if (post->GetStepStatus() == fGeomBoundary)   {
	    J2 = post->GetPosition();
	    momDir2x = track->GetMomentumDirection().x();
	    momDir2y = track->GetMomentumDirection().y();
	    momDir2z = track->GetMomentumDirection().z();

            fEventaction->AddMomentumDirection2(momDir2x, momDir2y, momDir2z);

            /*G4cout 
               << "\n Primary has exited BPIX1516 at: " 
               << G4BestUnit(J2, "Length")
               //<< "\n with momentum direction: " 
               //<< momDir2x
               //<< " "
               //<< momDir2y
               //<< " "
               //<< momDir2z
               << G4endl;*/
        }
    }
  }

  //track length of secondaries and tertiaries calculation
  if (track->GetParentID() == 1)  {
    fRunaction->AddSecTrackLength(stepSize);
    G4AnalysisManager::Instance()->FillH1(18,stepSize);

    fEventaction->AddSecondaryTrackLength(stepSize);

    if ((track->GetLogicalVolumeAtVertex() == fDetectorconstruction->GetLDet1()) || (track->GetLogicalVolumeAtVertex() == fDetectorconstruction->GetLDet2()))  {
         fRunaction->AddSecDetTrackLength(stepSize);
         G4AnalysisManager::Instance()->FillH1(20,stepSize);

         fEventaction->AddSecondaryDetTrackLength(stepSize);
         fEventaction->TrackCheck(length);

         //G4cout 
            //<< "\n Secondary produced at: " 
            //<< G4BestUnit(primVert, "Length")
            //<< "\n with track ID: " 
            //<< trackID
            //<< G4endl;

         //G4cout 
            //<< "\n Secondary produced inside a detector has a track length: " 
            //<< G4BestUnit(length, "Length")
            //<< G4endl;

         //G4cout 
            //<< "\n Secondary produced inside a detector's current position: " 
            //<< G4BestUnit(x, "Length")
            //<< G4BestUnit(y, "Length")
            //<< G4BestUnit(z, "Length")
            //<< G4endl;
    }

    if (track->GetLogicalVolumeAtVertex() == fDetectorconstruction->GetLDet1())  {
         fRunaction->AddSecDetTrackLength(stepSize);
         G4AnalysisManager::Instance()->FillH1(20,stepSize);    

         fEventaction->AddSecondaryDetTrackLength(stepSize);
         fEventaction->TrackCheck(length);

         //G4cout 
            //<< "\n Secondary produced at: " 
            //<< G4BestUnit(primVert, "Length")
            //<< "\n with track ID: " 
            //<< trackID
            //<< G4endl;

         //G4cout 
            //<< "\n Secondary produced inside Detector 1 has a track length: " 
            //<< G4BestUnit(length, "Length")
            //<< G4endl;

         //G4cout 
            //<< "\n Secondary produced inside a detector's current position: " 
            //<< G4BestUnit(x, "Length")
            //<< G4BestUnit(y, "Length")
            //<< G4BestUnit(z, "Length")
            //<< G4endl;
    }

    if (track->GetLogicalVolumeAtVertex() == fDetectorconstruction->GetLDet2())  {
         fRunaction->AddSecDetTrackLength(stepSize);
         G4AnalysisManager::Instance()->FillH1(20,stepSize);

         fEventaction->AddSecondaryDetTrackLength(stepSize);
         fEventaction->TrackCheck(length);

         //G4cout 
            //<< "\n Secondary produced at: " 
            //<< G4BestUnit(primVert, "Length")xAnglebpix
            //<< "\n with track ID: " 
            //<< trackID
            //<< G4endl;

         //G4cout 
            //<< "\n Secondary produced inside Detector 2 has a track length: " 
            //<< G4BestUnit(length, "Length")
            //<< G4endl;

         //G4cout 
            //<< "\n Secondary produced inside a detector's current position: " 
            //<< G4BestUnit(x, "Length")
            //<< G4BestUnit(y, "Length")
            //<< G4BestUnit(z, "Length")
            //<< G4endl;
    }

    if ((track->GetLogicalVolumeAtVertex() == fDetectorconstruction->GetLDet1()) && (volume == fDetectorconstruction->GetDet2()))  {
         //G4cout 
            //<< "\n Secondary produced inside Detector 1 has reached Detector 2. TrackID = " 
            //<< trackID
            //<< G4endl;

         if(particleDefinition == G4Positron::Definition())   {
              //G4cout 
                 //<< "\n It was a positron" 
                 //<< G4endl;
         }

         if(particleDefinition == G4Electron::Definition())   {
              //G4cout 
                 //<< "\n It was an electron" 
                 //<< G4endl;
         }

         if(particleDefinition == G4PionPlus::Definition())   {
              //G4cout 
                 //<< "\n It was a pi+" 
                 //<< G4endl;
         }

         if(particleDefinition == G4PionMinus::Definition())   {
              //G4cout 
                 //<< "\n It was a pi-" 
                 //<< G4endl;
         }

         if(particleDefinition == G4MuonPlus::Definition())   {
              //G4cout 
                 //<< "\n It was a mu+" 
                 //<< G4endl;
         }

         if(particleDefinition == G4MuonMinus::Definition())   {
              //G4cout 
                 //<< "\n It was a mu-" 
                 //<< G4endl;
         }

         if(particleDefinition == G4Gamma::Definition())   {
              //G4cout 
                 //<< "\n It was a photon" 
                 //<< G4endl;
         }


   	 if ((z >= Dist/2) && (z <= (Dist/2 + Strip2Depth))&& (x >= (-Strip2Length/2)) && (x <= Strip2Length/2)) {
     	  div = y/(StripWidth+StripDist);
	  if (y == 0)   {
	     iSecret = rand() % 99;
	     if (iSecret < 50)   {
		  stripNo = totalNbStrips/2 + div;
	     }
	     if (iSecret >= 50)   {
	  	  stripNo = totalNbStrips/2 + 1 + div;
	     }
	  }
     	  if (y > 0)    {
	     stripNo = totalNbStrips/2 + 1 + div;
     	  }
     	  if (y < 0)    {
             stripNo = totalNbStrips/2 + div;
     	  }
     	  if ((stripNo <= totalNbStrips) && (stripNo > 0))   {
     	     //G4cout
	         //<< "\n The secondary has passed through strip No "
	         //<< stripNo
                 //<< G4endl;
          }
         }
    }

    if ((track->GetLogicalVolumeAtVertex() == fDetectorconstruction->GetLDet2()) && (volume == fDetectorconstruction->GetDet1()))  {
         //G4cout 
            //<< "\n Secondary produced inside Detector 2 has reached Detector 1. TrackID = " 
            //<< trackID
            //<< G4endl;

         if(particleDefinition == G4Positron::Definition())   {
              //G4cout 
                 //<< "\n It was a positron" 
                 //<< G4endl;
         }

         if(particleDefinition == G4Electron::Definition())   {
              //G4cout 
                 //<< "\n It was an electron" 
                 //<< G4endl;
         }

         if(particleDefinition == G4PionPlus::Definition())   {
              //G4cout 
                 //<< "\n It was a pi+" 
                 //<< G4endl;
         }

         if(particleDefinition == G4PionMinus::Definition())   {
              //G4cout 
                 //<< "\n It was a pi-" 
                 //<< G4endl;
         }

         if(particleDefinition == G4MuonPlus::Definition())   {
              //G4cout 
                 //<< "\n It was a mu+" 
                 //<< G4endl;
         }

         if(particleDefinition == G4MuonMinus::Definition())   {
              //G4cout 
                 //<< "\n It was a mu-" 
                 //<< G4endl;
         }

         if(particleDefinition == G4Gamma::Definition())   {
              //G4cout 
                 //<< "\n It was a photon" 
                 //<< G4endl;
         }

    	 if ((z >= -Dist/2 - Det1SizeZ) && (z <= -Dist/2 - Det1SizeZ + Strip1Depth) && (x >= (-Strip1Length/2)) && (x <= Strip1Length/2)) {
     	  div = y/(StripWidth+StripDist);
	  if (y == 0)   {
	     iSecret = rand() % 99;
	     if (iSecret < 50)   {
		  stripNo = totalNbStrips/2 + div;
	     }
	     if (iSecret >= 50)   {
	       	  stripNo = totalNbStrips/2 + 1 + div;
	     }
	  }
     	  if (y > 0)    {
	      stripNo = totalNbStrips/2 + 1 + div;
     	  }
     	  if (y < 0)    {
              stripNo = totalNbStrips/2 + div;
     	  }
     	  if ((stripNo <= totalNbStrips) && (stripNo > 0))  {
     	      //G4cout
	   	  //<< "\n The secondary has passed through strip No "
	   	  //<< stripNo
           	  //<< G4endl;
     	  }
    	 }
    }
  }

  if (track->GetParentID() > 1)  {
    fRunaction->AddTertTrackLength(stepSize);
    G4AnalysisManager::Instance()->FillH1(19,stepSize);

    fEventaction->AddTertiaryTrackLength(stepSize);
  }

 //continuous energy deposit per event  
 if (volume == fDetectorconstruction->GetBPIX12()) {

   fEventaction->AddEnergyDepositL1(edep);

   mod = -1;
   ROC = -1;
   row = -1;
   col = -1;
   dirx = -1;
   dirxM = -1;
   divp1 = -1;

   zPixTr = z_after - poszBPIX12;
   xPix1 = x_after*cos(YangleBPIX) + zPixTr*sin(YangleBPIX);
   yPix1 = y_after;
   zPix1 = -x_after*sin(YangleBPIX) + zPixTr*cos(YangleBPIX);
   xPix = xPix1;
   yPix = yPix1*cos(XangleBPIX) - zPix1*sin(XangleBPIX);
   zPix = yPix1*sin(XangleBPIX) + zPix1*cos(XangleBPIX);

   if ((zPix <= (BPIXSizeZ/2))  &&  (zPix >= (BPIXSizeZ/2 - pixelDepth))) {
     fEventaction->AddEnergyDepositL1Active(edep);
     /*G4cout
       << "\n zPix12 = "
       << G4BestUnit(zPix, "Length")
       << G4endl; */
     if (xPix == 0)   {
	iSecret = rand() % 99;
     	if (iSecret < 50)   {
	   dirx = 0;
     	}
     	if (iSecret >= 50)   {
   	   dirx = 1;
     	}
     }
     if (xPix > 0)    {
	dirx = 0;
     }
     if (xPix < 0)    {
        dirx = 1;
     }
     if (dirx == 0)  {
        mod = 1; //module ID
	xPixM = xPix - BPIXSizeX/2;
	/*G4cout
	   << "\n xPixM = "
	   << G4BestUnit(xPixM, "Length")
	   << G4endl; */
        if (xPixM == 0)   {
	  iSecret = rand() % 99;
     	  if (iSecret < 50)   {
	     dirxM = 0;
     	  }
     	  if (iSecret >= 50)   {
   	     dirxM = 1;
     	  }
        }
     	if (xPixM > 0)    {
	  dirxM = 0;
     	}
     	if (xPixM < 0)    {
          dirxM = 1;
        }
     }
     if (dirx == 1)  {
        mod = 2; //module ID
	xPixM = xPix + BPIXSizeX/2;
	/*G4cout
	   << "\n xPixM = "
	   << G4BestUnit(xPixM, "Length")
	   << G4endl; */
        if (xPixM == 0)   {
	  iSecret = rand() % 99;
     	  if (iSecret < 50)   {
	     dirxM = 0;
     	  }
     	  if (iSecret >= 50)   {
   	     dirxM = 1;
     	  }
        }
     	if (xPixM > 0)    {
	  dirxM = 0;
     	}
     	if (xPixM < 0)    {
          dirxM = 1;
        }
     }
     div1 = xPix/ROChor;
     div2 = yPix/ROCvert;
     divp1 = xPixM/pixelX;
     ROC = 0;
     if (yPix == 0)   {
	iSecret = rand() % 99;
     	if (iSecret < 50)   {
	   diry = 0;
     	}
     	if (iSecret >= 50)   {
   	   diry = 1;
     	}
     }
     if (yPix > 0)    {
	diry = 0;
     }
     if (yPix < 0)    {
        diry = 1;
     }
     if ((dirxM == 0) && (diry == 0))   {
         if (div2 == 3)  {
             ROC = 1;
         }
         if (div2 == 2)  {
             ROC = 2;
         }
         if (div2 == 1)  {
             ROC = 3;
         }
         if (div2 == 0)  {
             ROC = 4;
         }
     }
     if ((dirxM == 0) && (diry == 1))   {
         if (div2 == 0)  {
             ROC = 5;
         }
         if (div2 == -1)  {
             ROC = 6;
         }  
         if (div2 == -2)  {
             ROC = 7;
         }
         if (div2 == -3)  {
             ROC = 8;
         }
     }
     if ((dirxM == 1) && (diry == 0))   {
         if (div2 == 3)  {
             ROC = 9;
         }
         if (div2 == 2)  {
             ROC = 10;
         }
         if (div2 == 1)  {
             ROC = 11;
         }
         if (div2 == 0)  {
             ROC = 12;
         }
     }
     if ((dirxM == 1) && (diry == 1))   {
         if (div2 == 0)  {
             ROC = 13;
         }
         if (div2 == -1)  {
             ROC = 14;
         }
         if (div2 == -2)  {
             ROC = 15;
         }
         if (div2 == -3)  {
             ROC = 16;
         }
     }
     if (dirxM == 0)  {
	 if ((divp1 == 0) || (divp1 == 1))  {
	    col = 80;
         }
         if ((divp1 > 1) && (divp1 <= 80))  {
	    col = 80 - divp1 + 1;
         } 
     }
     if (dirxM == 1)  {
	 if ((divp1 == -79) || (divp1 == -80))  {
	    col = 80;
         }
         if ((divp1 > -79) && (divp1 <= 0))  {
	    col = -divp1 + 1;
         } 
     }

     if ((ROC == 1) || (ROC == 9))  {
     	divp2 = (yPix-3*ROCvert)/pixelY;
     }
     if ((ROC == 2) || (ROC == 10))  {
     	divp2 = (yPix-2*ROCvert)/pixelY;
     }
     if ((ROC == 3) || (ROC == 11))  {
     	divp2 = (yPix-1*ROCvert)/pixelY;
     }
     if ((ROC == 4) || (ROC == 12))  {
     	divp2 = (yPix-0*ROCvert)/pixelY;
     }

     if ((ROC == 5) || (ROC == 13))  {
     	divp2 = (yPix+0*ROCvert)/pixelY;
     }
     if ((ROC == 6) || (ROC == 14))  {
     	divp2 = (yPix+1*ROCvert)/pixelY;
     }
     if ((ROC == 7) || (ROC == 15))  {
     	divp2 = (yPix+2*ROCvert)/pixelY;
     }
     if ((ROC == 8) || (ROC == 16))  {
     	divp2 = (yPix+3*ROCvert)/pixelY;
     }

     if ((ROC == 1) || (ROC == 2) || (ROC == 3) || (ROC == 4) || (ROC == 9) || (ROC == 10) || (ROC == 11) || (ROC == 12))  {
	if ((divp2 == 0) || (divp2 == 1))  {
	   row = 52;
	}
	if ((divp2 == 52) || (divp2 == 53))  {
	   row = 1;
	}
	if ((divp2 < 52) && (divp2 > 1))  {
	   row = 52 - divp2 + 1;
	}
     }
     if ((ROC == 5) || (ROC == 6) || (ROC == 7) || (ROC == 8) || (ROC == 13) || (ROC == 14) || (ROC == 15) || (ROC == 16))  {
	if ((divp2 == 0) || (divp2 == -1))  {
	   row = 1;
	}
	if ((divp2 == -52) || (divp2 == -53))  {
	   row = 52;
	}
	if ((divp2 > -52) && (divp2 < -1))  {
	   row = -divp2;
	}
     }

     if ((mod>=1) && (mod<17) && (ROC>=1) && (ROC<17) && (row>=1) && (row<53) && (col>=1) && (col<81))   {
	   fEventaction->AddEnergyPixel(edep,mod,ROC,row,col);
	   /*G4cout
	      << "\n Enery added to pixel = "
	      << G4BestUnit(edep, "Energy")
	      << "\n module = "
	      << mod
	      << G4endl;*/ 
     }   
   }

   /*G4cout
     << "\n edep = "
     << G4BestUnit(edep, "Energy")
     << "\n module = "
     << mod
     << "\n ROC = "
     << ROC
     << "\n row = "
     << row
     << "\n col = "
     << col
     << "\n diry = "
     << diry
     << "\n dirxM = "
     << dirxM
     << "\n divp1 = "
     << divp1
     << G4endl;*/ 
 } 

 if (volume == fDetectorconstruction->GetBPIX34()) {
   fEventaction->AddEnergyDepositL2(edep);

   mod = -1;
   ROC = -1;
   row = -1;
   col = -1;
   dirx = -1;
   dirxM = -1;
   divp1 = -1;

   zPixTr = z_after - poszBPIX34;
   xPix1 = x_after*cos(YangleBPIX) + zPixTr*sin(YangleBPIX);
   yPix1 = y_after;
   zPix1 = -x_after*sin(YangleBPIX) + zPixTr*cos(YangleBPIX);
   xPix = xPix1;
   yPix = yPix1*cos(XangleBPIX) - zPix1*sin(XangleBPIX);
   zPix = yPix1*sin(XangleBPIX) + zPix1*cos(XangleBPIX);

   if ((zPix <= (BPIXSizeZ/2))  &&  (zPix >= (BPIXSizeZ/2 - pixelDepth))) {
     /*G4cout
       << "\n zPix34 = "
       << G4BestUnit(zPix, "Length")
       << G4endl; */
     if (xPix == 0)   {
	iSecret = rand() % 99;
     	if (iSecret < 50)   {
	   dirx = 0;
     	}
     	if (iSecret >= 50)   {
   	   dirx = 1;
     	}
     }
     if (xPix > 0)    {
	dirx = 0;
     }
     if (xPix < 0)    {
        dirx = 1;
     }
     if (dirx == 0)  {
        mod = 3; //module ID
	xPixM = xPix - BPIXSizeX/2;
	/*G4cout
	   << "\n xPixM = "
	   << G4BestUnit(xPixM, "Length")
	   << G4endl; */
        if (xPixM == 0)   {
	  iSecret = rand() % 99;
     	  if (iSecret < 50)   {
	     dirxM = 0;
     	  }
     	  if (iSecret >= 50)   {
   	     dirxM = 1;
     	  }
        }
     	if (xPixM > 0)    {
	  dirxM = 0;
     	}
     	if (xPixM < 0)    {
          dirxM = 1;
        }
     }
     if (dirx == 1)  {
        mod = 4; //module ID
	xPixM = xPix + BPIXSizeX/2;
	/*G4cout
	   << "\n xPixM = "
	   << G4BestUnit(xPixM, "Length")
	   << G4endl; */
        if (xPixM == 0)   {
	  iSecret = rand() % 99;
     	  if (iSecret < 50)   {
	     dirxM = 0;
     	  }
     	  if (iSecret >= 50)   {
   	     dirxM = 1;
     	  }
        }
     	if (xPixM > 0)    {
	  dirxM = 0;
     	}
     	if (xPixM < 0)    {
          dirxM = 1;
        }
     }
     div1 = xPix/ROChor;
     div2 = yPix/ROCvert;
     divp1 = xPixM/pixelX;
     ROC = 0;
     if (yPix == 0)   {
	iSecret = rand() % 99;
     	if (iSecret < 50)   {
	   diry = 0;
     	}
     	if (iSecret >= 50)   {
   	   diry = 1;
     	}
     }
     if (yPix > 0)    {
	diry = 0;
     }
     if (yPix < 0)    {
        diry = 1;
     }
     if ((dirxM == 0) && (diry == 0))   {
         if (div2 == 3)  {
             ROC = 1;
         }
         if (div2 == 2)  {
             ROC = 2;
         }
         if (div2 == 1)  {
             ROC = 3;
         }
         if (div2 == 0)  {
             ROC = 4;
         }
     }
     if ((dirxM == 0) && (diry == 1))   {
         if (div2 == 0)  {
             ROC = 5;
         }
         if (div2 == -1)  {
             ROC = 6;
         }  
         if (div2 == -2)  {
             ROC = 7;
         }
         if (div2 == -3)  {
             ROC = 8;
         }
     }
     if ((dirxM == 1) && (diry == 0))   {
         if (div2 == 3)  {
             ROC = 9;
         }
         if (div2 == 2)  {
             ROC = 10;
         }
         if (div2 == 1)  {
             ROC = 11;
         }
         if (div2 == 0)  {
             ROC = 12;
         }
     }
     if ((dirxM == 1) && (diry == 1))   {
         if (div2 == 0)  {
             ROC = 13;
         }
         if (div2 == -1)  {
             ROC = 14;
         }
         if (div2 == -2)  {
             ROC = 15;
         }
         if (div2 == -3)  {
             ROC = 16;
         }
     }
     if (dirxM == 0)  {
	 if ((divp1 == 0) || (divp1 == 1))  {
	    col = 80;
         }
         if ((divp1 > 1) && (divp1 <= 80))  {
	    col = 80 - divp1 + 1;
         } 
     }
     if (dirxM == 1)  {
	 if ((divp1 == -79) || (divp1 == -80))  {
	    col = 80;
         }
         if ((divp1 > -79) && (divp1 <= 0))  {
	    col = -divp1 + 1;
         } 
     }

     if ((ROC == 1) || (ROC == 9))  {
     	divp2 = (yPix-3*ROCvert)/pixelY;
     }
     if ((ROC == 2) || (ROC == 10))  {
     	divp2 = (yPix-2*ROCvert)/pixelY;
     }
     if ((ROC == 3) || (ROC == 11))  {
     	divp2 = (yPix-1*ROCvert)/pixelY;
     }
     if ((ROC == 4) || (ROC == 12))  {
     	divp2 = (yPix-0*ROCvert)/pixelY;
     }

     if ((ROC == 5) || (ROC == 13))  {
     	divp2 = (yPix+0*ROCvert)/pixelY;
     }
     if ((ROC == 6) || (ROC == 14))  {
     	divp2 = (yPix+1*ROCvert)/pixelY;
     }
     if ((ROC == 7) || (ROC == 15))  {
     	divp2 = (yPix+2*ROCvert)/pixelY;
     }
     if ((ROC == 8) || (ROC == 16))  {
     	divp2 = (yPix+3*ROCvert)/pixelY;
     }

     if ((ROC == 1) || (ROC == 2) || (ROC == 3) || (ROC == 4) || (ROC == 9) || (ROC == 10) || (ROC == 11) || (ROC == 12))  {
	if ((divp2 == 0) || (divp2 == 1))  {
	   row = 52;
	}
	if ((divp2 == 52) || (divp2 == 53))  {
	   row = 1;
	}
	if ((divp2 < 52) && (divp2 > 1))  {
	   row = 52 - divp2 + 1;
	}
     }
     if ((ROC == 5) || (ROC == 6) || (ROC == 7) || (ROC == 8) || (ROC == 13) || (ROC == 14) || (ROC == 15) || (ROC == 16))  {
	if ((divp2 == 0) || (divp2 == -1))  {
	   row = 1;
	}
	if ((divp2 == -52) || (divp2 == -53))  {
	   row = 52;
	}
	if ((divp2 > -52) && (divp2 < -1))  {
	   row = -divp2;
	}
     }

     if ((mod>=1) && (mod<17) && (ROC>=1) && (ROC<17) && (row>=1) && (row<53) && (col>=1) && (col<81))   {
	   fEventaction->AddEnergyPixel(edep,mod,ROC,row,col);
	   /*G4cout
	      << "\n Enery added to pixel = "
	      << G4BestUnit(edep, "Energy")
	      << "\n module = "
	      << mod
	      << G4endl;*/ 
     }   
   }

   /*G4cout
     << "\n edep = "
     << G4BestUnit(edep, "Energy")
     << "\n module = "
     << mod
     << "\n ROC = "
     << ROC
     << "\n row = "
     << row
     << "\n col = "
     << col
     << "\n diry = "
     << diry
     << "\n dirxM = "
     << dirxM
     << "\n divp1 = "
     << divp1
     << G4endl;*/ 
 }

 if (volume == fDetectorconstruction->GetBPIX56()) {
   fEventaction->AddEnergyDepositL3(edep);

   mod = -1;
   ROC = -1;
   row = -1;
   col = -1;
   dirx = -1;
   dirxM = -1;
   divp1 = -1;

   zPixTr = z_after - poszBPIX56;
   xPix1 = x_after*cos(YangleBPIX) + zPixTr*sin(YangleBPIX);
   yPix1 = y_after;
   zPix1 = -x_after*sin(YangleBPIX) + zPixTr*cos(YangleBPIX);
   xPix = xPix1;
   yPix = yPix1*cos(XangleBPIX) - zPix1*sin(XangleBPIX);
   zPix = yPix1*sin(XangleBPIX) + zPix1*cos(XangleBPIX);

   if ((zPix <= (BPIXSizeZ/2))  &&  (zPix >= (BPIXSizeZ/2 - pixelDepth))) {
     /*G4cout
       << "\n zPix56 = "
       << G4BestUnit(zPix, "Length")
       << G4endl; */
     if (xPix == 0)   {
	iSecret = rand() % 99;
     	if (iSecret < 50)   {
	   dirx = 0;
     	}
     	if (iSecret >= 50)   {
   	   dirx = 1;
     	}
     }
     if (xPix > 0)    {
	dirx = 0;
     }
     if (xPix < 0)    {
        dirx = 1;
     }
     if (dirx == 0)  {
        mod = 5; //module ID
	xPixM = xPix - BPIXSizeX/2;
	/*G4cout
	   << "\n xPixM = "
	   << G4BestUnit(xPixM, "Length")
	   << G4endl; */
        if (xPixM == 0)   {
	  iSecret = rand() % 99;
     	  if (iSecret < 50)   {
	     dirxM = 0;
     	  }
     	  if (iSecret >= 50)   {
   	     dirxM = 1;
     	  }
        }
     	if (xPixM > 0)    {
	  dirxM = 0;
     	}
     	if (xPixM < 0)    {
          dirxM = 1;
        }
     }
     if (dirx == 1)  {
        mod = 6; //module ID
	xPixM = xPix + BPIXSizeX/2;
	/*G4cout
	   << "\n xPixM = "
	   << G4BestUnit(xPixM, "Length")
	   << G4endl; */
        if (xPixM == 0)   {
	  iSecret = rand() % 99;
     	  if (iSecret < 50)   {
	     dirxM = 0;
     	  }
     	  if (iSecret >= 50)   {
   	     dirxM = 1;
     	  }
        }
     	if (xPixM > 0)    {
	  dirxM = 0;
     	}
     	if (xPixM < 0)    {
          dirxM = 1;
        }
     }
     div1 = xPix/ROChor;
     div2 = yPix/ROCvert;
     divp1 = xPixM/pixelX;
     ROC = 0;
     if (yPix == 0)   {
	iSecret = rand() % 99;
     	if (iSecret < 50)   {
	   diry = 0;
     	}
     	if (iSecret >= 50)   {
   	   diry = 1;
     	}
     }
     if (yPix > 0)    {
	diry = 0;
     }
     if (yPix < 0)    {
        diry = 1;
     }
     if ((dirxM == 0) && (diry == 0))   {
         if (div2 == 3)  {
             ROC = 1;
         }
         if (div2 == 2)  {
             ROC = 2;
         }
         if (div2 == 1)  {
             ROC = 3;
         }
         if (div2 == 0)  {
             ROC = 4;
         }
     }
     if ((dirxM == 0) && (diry == 1))   {
         if (div2 == 0)  {
             ROC = 5;
         }
         if (div2 == -1)  {
             ROC = 6;
         }  
         if (div2 == -2)  {
             ROC = 7;
         }
         if (div2 == -3)  {
             ROC = 8;
         }
     }
     if ((dirxM == 1) && (diry == 0))   {
         if (div2 == 3)  {
             ROC = 9;
         }
         if (div2 == 2)  {
             ROC = 10;
         }
         if (div2 == 1)  {
             ROC = 11;
         }
         if (div2 == 0)  {
             ROC = 12;
         }
     }
     if ((dirxM == 1) && (diry == 1))   {
         if (div2 == 0)  {
             ROC = 13;
         }
         if (div2 == -1)  {
             ROC = 14;
         }
         if (div2 == -2)  {
             ROC = 15;
         }
         if (div2 == -3)  {
             ROC = 16;
         }
     }
     if (dirxM == 0)  {
	 if ((divp1 == 0) || (divp1 == 1))  {
	    col = 80;
         }
         if ((divp1 > 1) && (divp1 <= 80))  {
	    col = 80 - divp1 + 1;
         } 
     }
     if (dirxM == 1)  {
	 if ((divp1 == -79) || (divp1 == -80))  {
	    col = 80;
         }
         if ((divp1 > -79) && (divp1 <= 0))  {
	    col = -divp1 + 1;
         } 
     }

     if ((ROC == 1) || (ROC == 9))  {
     	divp2 = (yPix-3*ROCvert)/pixelY;
     }
     if ((ROC == 2) || (ROC == 10))  {
     	divp2 = (yPix-2*ROCvert)/pixelY;
     }
     if ((ROC == 3) || (ROC == 11))  {
     	divp2 = (yPix-1*ROCvert)/pixelY;
     }
     if ((ROC == 4) || (ROC == 12))  {
     	divp2 = (yPix-0*ROCvert)/pixelY;
     }

     if ((ROC == 5) || (ROC == 13))  {
     	divp2 = (yPix+0*ROCvert)/pixelY;
     }
     if ((ROC == 6) || (ROC == 14))  {
     	divp2 = (yPix+1*ROCvert)/pixelY;
     }
     if ((ROC == 7) || (ROC == 15))  {
     	divp2 = (yPix+2*ROCvert)/pixelY;
     }
     if ((ROC == 8) || (ROC == 16))  {
     	divp2 = (yPix+3*ROCvert)/pixelY;
     }

     if ((ROC == 1) || (ROC == 2) || (ROC == 3) || (ROC == 4) || (ROC == 9) || (ROC == 10) || (ROC == 11) || (ROC == 12))  {
	if ((divp2 == 0) || (divp2 == 1))  {
	   row = 52;
	}
	if ((divp2 == 52) || (divp2 == 53))  {
	   row = 1;
	}
	if ((divp2 < 52) && (divp2 > 1))  {
	   row = 52 - divp2 + 1;
	}
     }
     if ((ROC == 5) || (ROC == 6) || (ROC == 7) || (ROC == 8) || (ROC == 13) || (ROC == 14) || (ROC == 15) || (ROC == 16))  {
	if ((divp2 == 0) || (divp2 == -1))  {
	   row = 1;
	}
	if ((divp2 == -52) || (divp2 == -53))  {
	   row = 52;
	}
	if ((divp2 > -52) && (divp2 < -1))  {
	   row = -divp2;
	}
     }

     if ((mod>=1) && (mod<17) && (ROC>=1) && (ROC<17) && (row>=1) && (row<53) && (col>=1) && (col<81))   {
	   fEventaction->AddEnergyPixel(edep,mod,ROC,row,col);
	   /*G4cout
	      << "\n Enery added to pixel = "
	      << G4BestUnit(edep, "Energy")
	      << "\n module = "
	      << mod
	      << G4endl;*/ 
     }   
   }

   /*G4cout
     << "\n edep = "
     << G4BestUnit(edep, "Energy")
     << "\n module = "
     << mod
     << "\n ROC = "
     << ROC
     << "\n row = "
     << row
     << "\n col = "
     << col
     << "\n diry = "
     << diry
     << "\n dirxM = "
     << dirxM
     << "\n divp1 = "
     << divp1
     << G4endl;*/ 
 }

 if (volume == fDetectorconstruction->GetBPIX78()) {
   fEventaction->AddEnergyDepositL4(edep);

   mod = -1;
   ROC = -1;
   row = -1;
   col = -1;
   dirx = -1;
   dirxM = -1;
   divp1 = -1;

   zPixTr = z_after - poszBPIX78;
   xPix1 = x_after*cos(YangleBPIX) + zPixTr*sin(YangleBPIX);
   yPix1 = y_after;
   zPix1 = -x_after*sin(YangleBPIX) + zPixTr*cos(YangleBPIX);
   xPix = xPix1;
   yPix = yPix1*cos(XangleBPIX) - zPix1*sin(XangleBPIX);
   zPix = yPix1*sin(XangleBPIX) + zPix1*cos(XangleBPIX);

   if ((zPix <= (BPIXSizeZ/2))  &&  (zPix >= (BPIXSizeZ/2 - pixelDepth))) {
     /*G4cout
       << "\n zPix78 = "
       << G4BestUnit(zPix, "Length")
       << G4endl; */
     if (xPix == 0)   {
	iSecret = rand() % 99;
     	if (iSecret < 50)   {
	   dirx = 0;
     	}
     	if (iSecret >= 50)   {
   	   dirx = 1;
     	}
     }
     if (xPix > 0)    {
	dirx = 0;
     }
     if (xPix < 0)    {
        dirx = 1;
     }
     if (dirx == 0)  {
        mod = 7; //module ID
	xPixM = xPix - BPIXSizeX/2;
	/*G4cout
	   << "\n xPixM = "
	   << G4BestUnit(xPixM, "Length")
	   << G4endl; */
        if (xPixM == 0)   {
	  iSecret = rand() % 99;
     	  if (iSecret < 50)   {
	     dirxM = 0;
     	  }
     	  if (iSecret >= 50)   {
   	     dirxM = 1;
     	  }
        }
     	if (xPixM > 0)    {
	  dirxM = 0;
     	}
     	if (xPixM < 0)    {
          dirxM = 1;
        }
     }
     if (dirx == 1)  {
        mod = 8; //module ID
	xPixM = xPix + BPIXSizeX/2;
	/*G4cout
	   << "\n xPixM = "
	   << G4BestUnit(xPixM, "Length")
	   << G4endl; */
        if (xPixM == 0)   {
	  iSecret = rand() % 99;
     	  if (iSecret < 50)   {
	     dirxM = 0;
     	  }
     	  if (iSecret >= 50)   {
   	     dirxM = 1;
     	  }
        }
     	if (xPixM > 0)    {
	  dirxM = 0;
     	}
     	if (xPixM < 0)    {
          dirxM = 1;
        }
     }
     div1 = xPix/ROChor;
     div2 = yPix/ROCvert;
     divp1 = xPixM/pixelX;
     ROC = 0;
     if (yPix == 0)   {
	iSecret = rand() % 99;
     	if (iSecret < 50)   {
	   diry = 0;
     	}
     	if (iSecret >= 50)   {
   	   diry = 1;
     	}
     }
     if (yPix > 0)    {
	diry = 0;
     }
     if (yPix < 0)    {
        diry = 1;
     }
     if ((dirxM == 0) && (diry == 0))   {
         if (div2 == 3)  {
             ROC = 1;
         }
         if (div2 == 2)  {
             ROC = 2;
         }
         if (div2 == 1)  {
             ROC = 3;
         }
         if (div2 == 0)  {
             ROC = 4;
         }
     }
     if ((dirxM == 0) && (diry == 1))   {
         if (div2 == 0)  {
             ROC = 5;
         }
         if (div2 == -1)  {
             ROC = 6;
         }  
         if (div2 == -2)  {
             ROC = 7;
         }
         if (div2 == -3)  {
             ROC = 8;
         }
     }
     if ((dirxM == 1) && (diry == 0))   {
         if (div2 == 3)  {
             ROC = 9;
         }
         if (div2 == 2)  {
             ROC = 10;
         }
         if (div2 == 1)  {
             ROC = 11;
         }
         if (div2 == 0)  {
             ROC = 12;
         }
     }
     if ((dirxM == 1) && (diry == 1))   {
         if (div2 == 0)  {
             ROC = 13;
         }
         if (div2 == -1)  {
             ROC = 14;
         }
         if (div2 == -2)  {
             ROC = 15;
         }
         if (div2 == -3)  {
             ROC = 16;
         }
     }
     if (dirxM == 0)  {
	 if ((divp1 == 0) || (divp1 == 1))  {
	    col = 80;
         }
         if ((divp1 > 1) && (divp1 <= 80))  {
	    col = 80 - divp1 + 1;
         } 
     }
     if (dirxM == 1)  {
	 if ((divp1 == -79) || (divp1 == -80))  {
	    col = 80;
         }
         if ((divp1 > -79) && (divp1 <= 0))  {
	    col = -divp1 + 1;
         } 
     }

     if ((ROC == 1) || (ROC == 9))  {
     	divp2 = (yPix-3*ROCvert)/pixelY;
     }
     if ((ROC == 2) || (ROC == 10))  {
     	divp2 = (yPix-2*ROCvert)/pixelY;
     }
     if ((ROC == 3) || (ROC == 11))  {
     	divp2 = (yPix-1*ROCvert)/pixelY;
     }
     if ((ROC == 4) || (ROC == 12))  {
     	divp2 = (yPix-0*ROCvert)/pixelY;
     }

     if ((ROC == 5) || (ROC == 13))  {
     	divp2 = (yPix+0*ROCvert)/pixelY;
     }
     if ((ROC == 6) || (ROC == 14))  {
     	divp2 = (yPix+1*ROCvert)/pixelY;
     }
     if ((ROC == 7) || (ROC == 15))  {
     	divp2 = (yPix+2*ROCvert)/pixelY;
     }
     if ((ROC == 8) || (ROC == 16))  {
     	divp2 = (yPix+3*ROCvert)/pixelY;
     }

     if ((ROC == 1) || (ROC == 2) || (ROC == 3) || (ROC == 4) || (ROC == 9) || (ROC == 10) || (ROC == 11) || (ROC == 12))  {
	if ((divp2 == 0) || (divp2 == 1))  {
	   row = 52;
	}
	if ((divp2 == 52) || (divp2 == 53))  {
	   row = 1;
	}
	if ((divp2 < 52) && (divp2 > 1))  {
	   row = 52 - divp2 + 1;
	}
     }
     if ((ROC == 5) || (ROC == 6) || (ROC == 7) || (ROC == 8) || (ROC == 13) || (ROC == 14) || (ROC == 15) || (ROC == 16))  {
	if ((divp2 == 0) || (divp2 == -1))  {
	   row = 1;
	}
	if ((divp2 == -52) || (divp2 == -53))  {
	   row = 52;
	}
	if ((divp2 > -52) && (divp2 < -1))  {
	   row = -divp2;
	}
     }

     if ((mod>=1) && (mod<17) && (ROC>=1) && (ROC<17) && (row>=1) && (row<53) && (col>=1) && (col<81))   {
	   fEventaction->AddEnergyPixel(edep,mod,ROC,row,col);
	   /*G4cout
	      << "\n Enery added to pixel = "
	      << G4BestUnit(edep, "Energy")
	      << "\n module = "
	      << mod
	      << G4endl;*/ 
     }   
   }

   /*G4cout
     << "\n edep = "
     << G4BestUnit(edep, "Energy")
     << "\n module = "
     << mod
     << "\n ROC = "
     << ROC
     << "\n row = "
     << row
     << "\n col = "
     << col
     << "\n diry = "
     << diry
     << "\n dirxM = "
     << dirxM
     << "\n divp1 = "
     << divp1
     << G4endl;*/ 
 }

 if (volume == fDetectorconstruction->GetBPIX910()) {
   fEventaction->AddEnergyDepositL5(edep);

   mod = -1;
   ROC = -1;
   row = -1;
   col = -1;
   dirx = -1;
   dirxM = -1;
   divp1 = -1;

   zPixTr = z_after - poszBPIX910;
   xPix1 = x_after*cos(YangleBPIX) + zPixTr*sin(YangleBPIX);
   yPix1 = y_after;
   zPix1 = -x_after*sin(YangleBPIX) + zPixTr*cos(YangleBPIX);
   xPix = xPix1;
   yPix = yPix1*cos(XangleBPIX) - zPix1*sin(XangleBPIX);
   zPix = yPix1*sin(XangleBPIX) + zPix1*cos(XangleBPIX);

   if ((zPix <= (BPIXSizeZ/2))  &&  (zPix >= (BPIXSizeZ/2 - pixelDepth))) {
     /*G4cout
       << "\n zPix910 = "
       << G4BestUnit(zPix, "Length")
       << G4endl; */
     if (xPix == 0)   {
	iSecret = rand() % 99;
     	if (iSecret < 50)   {
	   dirx = 0;
     	}
     	if (iSecret >= 50)   {
   	   dirx = 1;
     	}
     }
     if (xPix > 0)    {
	dirx = 0;
     }
     if (xPix < 0)    {
        dirx = 1;
     }
     if (dirx == 0)  {
        mod = 9; //module ID
	xPixM = xPix - BPIXSizeX/2;
	/*G4cout
	   << "\n xPixM = "
	   << G4BestUnit(xPixM, "Length")
	   << G4endl; */
        if (xPixM == 0)   {
	  iSecret = rand() % 99;
     	  if (iSecret < 50)   {
	     dirxM = 0;
     	  }
     	  if (iSecret >= 50)   {
   	     dirxM = 1;
     	  }
        }
     	if (xPixM > 0)    {
	  dirxM = 0;
     	}
     	if (xPixM < 0)    {
          dirxM = 1;
        }
     }
     if (dirx == 1)  {
        mod = 10; //module ID
	xPixM = xPix + BPIXSizeX/2;
	/*G4cout
	   << "\n xPixM = "
	   << G4BestUnit(xPixM, "Length")
	   << G4endl; */
        if (xPixM == 0)   {
	  iSecret = rand() % 99;
     	  if (iSecret < 50)   {
	     dirxM = 0;
     	  }
     	  if (iSecret >= 50)   {
   	     dirxM = 1;
     	  }
        }
     	if (xPixM > 0)    {
	  dirxM = 0;
     	}
     	if (xPixM < 0)    {
          dirxM = 1;
        }
     }
     div1 = xPix/ROChor;
     div2 = yPix/ROCvert;
     divp1 = xPixM/pixelX;
     ROC = 0;
     if (yPix == 0)   {
	iSecret = rand() % 99;
     	if (iSecret < 50)   {
	   diry = 0;
     	}
     	if (iSecret >= 50)   {
   	   diry = 1;
     	}
     }
     if (yPix > 0)    {
	diry = 0;
     }
     if (yPix < 0)    {
        diry = 1;
     }
     if ((dirxM == 0) && (diry == 0))   {
         if (div2 == 3)  {
             ROC = 1;
         }
         if (div2 == 2)  {
             ROC = 2;
         }
         if (div2 == 1)  {
             ROC = 3;
         }
         if (div2 == 0)  {
             ROC = 4;
         }
     }
     if ((dirxM == 0) && (diry == 1))   {
         if (div2 == 0)  {
             ROC = 5;
         }
         if (div2 == -1)  {
             ROC = 6;
         }  
         if (div2 == -2)  {
             ROC = 7;
         }
         if (div2 == -3)  {
             ROC = 8;
         }
     }
     if ((dirxM == 1) && (diry == 0))   {
         if (div2 == 3)  {
             ROC = 9;
         }
         if (div2 == 2)  {
             ROC = 10;
         }
         if (div2 == 1)  {
             ROC = 11;
         }
         if (div2 == 0)  {
             ROC = 12;
         }
     }
     if ((dirxM == 1) && (diry == 1))   {
         if (div2 == 0)  {
             ROC = 13;
         }
         if (div2 == -1)  {
             ROC = 14;
         }
         if (div2 == -2)  {
             ROC = 15;
         }
         if (div2 == -3)  {
             ROC = 16;
         }
     }
     if (dirxM == 0)  {
	 if ((divp1 == 0) || (divp1 == 1))  {
	    col = 80;
         }
         if ((divp1 > 1) && (divp1 <= 80))  {
	    col = 80 - divp1 + 1;
         } 
     }
     if (dirxM == 1)  {
	 if ((divp1 == -79) || (divp1 == -80))  {
	    col = 80;
         }
         if ((divp1 > -79) && (divp1 <= 0))  {
	    col = -divp1 + 1;
         } 
     }

     if ((ROC == 1) || (ROC == 9))  {
     	divp2 = (yPix-3*ROCvert)/pixelY;
     }
     if ((ROC == 2) || (ROC == 10))  {
     	divp2 = (yPix-2*ROCvert)/pixelY;
     }
     if ((ROC == 3) || (ROC == 11))  {
     	divp2 = (yPix-1*ROCvert)/pixelY;
     }
     if ((ROC == 4) || (ROC == 12))  {
     	divp2 = (yPix-0*ROCvert)/pixelY;
     }

     if ((ROC == 5) || (ROC == 13))  {
     	divp2 = (yPix+0*ROCvert)/pixelY;
     }
     if ((ROC == 6) || (ROC == 14))  {
     	divp2 = (yPix+1*ROCvert)/pixelY;
     }
     if ((ROC == 7) || (ROC == 15))  {
     	divp2 = (yPix+2*ROCvert)/pixelY;
     }
     if ((ROC == 8) || (ROC == 16))  {
     	divp2 = (yPix+3*ROCvert)/pixelY;
     }

     if ((ROC == 1) || (ROC == 2) || (ROC == 3) || (ROC == 4) || (ROC == 9) || (ROC == 10) || (ROC == 11) || (ROC == 12))  {
	if ((divp2 == 0) || (divp2 == 1))  {
	   row = 52;
	}
	if ((divp2 == 52) || (divp2 == 53))  {
	   row = 1;
	}
	if ((divp2 < 52) && (divp2 > 1))  {
	   row = 52 - divp2 + 1;
	}
     }
     if ((ROC == 5) || (ROC == 6) || (ROC == 7) || (ROC == 8) || (ROC == 13) || (ROC == 14) || (ROC == 15) || (ROC == 16))  {
	if ((divp2 == 0) || (divp2 == -1))  {
	   row = 1;
	}
	if ((divp2 == -52) || (divp2 == -53))  {
	   row = 52;
	}
	if ((divp2 > -52) && (divp2 < -1))  {
	   row = -divp2;
	}
     }

     if ((mod>=1) && (mod<17) && (ROC>=1) && (ROC<17) && (row>=1) && (row<53) && (col>=1) && (col<81))   {
	   fEventaction->AddEnergyPixel(edep,mod,ROC,row,col);
	   /*G4cout
	      << "\n Enery added to pixel = "
	      << G4BestUnit(edep, "Energy")
	      << "\n module = "
	      << mod
	      << G4endl;*/ 
     }   
   }

   /*G4cout
     << "\n edep = "
     << G4BestUnit(edep, "Energy")
     << "\n module = "
     << mod
     << "\n ROC = "
     << ROC
     << "\n row = "
     << row
     << "\n col = "
     << col
     << "\n diry = "
     << diry
     << "\n dirxM = "
     << dirxM
     << "\n divp1 = "
     << divp1
     << G4endl;*/ 
 }

 if (volume == fDetectorconstruction->GetBPIX1112()) {
   fEventaction->AddEnergyDepositL6(edep);

   mod = -1;
   ROC = -1;
   row = -1;
   col = -1;
   dirx = -1;
   dirxM = -1;
   divp1 = -1;

   zPixTr = z_after - poszBPIX1112;
   xPix1 = x_after*cos(YangleBPIX) + zPixTr*sin(YangleBPIX);
   yPix1 = y_after;
   zPix1 = -x_after*sin(YangleBPIX) + zPixTr*cos(YangleBPIX);
   xPix = xPix1;
   yPix = yPix1*cos(XangleBPIX) - zPix1*sin(XangleBPIX);
   zPix = yPix1*sin(XangleBPIX) + zPix1*cos(XangleBPIX);

   if ((zPix <= (BPIXSizeZ/2))  &&  (zPix >= (BPIXSizeZ/2 - pixelDepth))) {
     /*G4cout
       << "\n zPix1112 = "
       << G4BestUnit(zPix, "Length")
       << G4endl; */
     if (xPix == 0)   {
	iSecret = rand() % 99;
     	if (iSecret < 50)   {
	   dirx = 0;
     	}
     	if (iSecret >= 50)   {
   	   dirx = 1;
     	}
     }
     if (xPix > 0)    {
	dirx = 0;
     }
     if (xPix < 0)    {
        dirx = 1;
     }
     if (dirx == 0)  {
        mod = 11; //module ID
	xPixM = xPix - BPIXSizeX/2;
	/*G4cout
	   << "\n xPixM = "
	   << G4BestUnit(xPixM, "Length")
	   << G4endl; */
        if (xPixM == 0)   {
	  iSecret = rand() % 99;
     	  if (iSecret < 50)   {
	     dirxM = 0;
     	  }
     	  if (iSecret >= 50)   {
   	     dirxM = 1;
     	  }
        }
     	if (xPixM > 0)    {
	  dirxM = 0;
     	}
     	if (xPixM < 0)    {
          dirxM = 1;
        }
     }
     if (dirx == 1)  {
        mod = 12; //module ID
	xPixM = xPix + BPIXSizeX/2;
	/*G4cout
	   << "\n xPixM = "
	   << G4BestUnit(xPixM, "Length")
	   << G4endl; */
        if (xPixM == 0)   {
	  iSecret = rand() % 99;
     	  if (iSecret < 50)   {
	     dirxM = 0;
     	  }
     	  if (iSecret >= 50)   {
   	     dirxM = 1;
     	  }
        }
     	if (xPixM > 0)    {
	  dirxM = 0;
     	}
     	if (xPixM < 0)    {
          dirxM = 1;
        }
     }
     div1 = xPix/ROChor;
     div2 = yPix/ROCvert;
     divp1 = xPixM/pixelX;
     ROC = 0;
     if (yPix == 0)   {
	iSecret = rand() % 99;
     	if (iSecret < 50)   {
	   diry = 0;
     	}
     	if (iSecret >= 50)   {
   	   diry = 1;
     	}
     }
     if (yPix > 0)    {
	diry = 0;
     }
     if (yPix < 0)    {
        diry = 1;
     }
     if ((dirxM == 0) && (diry == 0))   {
         if (div2 == 3)  {
             ROC = 1;
         }
         if (div2 == 2)  {
             ROC = 2;
         }
         if (div2 == 1)  {
             ROC = 3;
         }
         if (div2 == 0)  {
             ROC = 4;
         }
     }
     if ((dirxM == 0) && (diry == 1))   {
         if (div2 == 0)  {
             ROC = 5;
         }
         if (div2 == -1)  {
             ROC = 6;
         }  
         if (div2 == -2)  {
             ROC = 7;
         }
         if (div2 == -3)  {
             ROC = 8;
         }
     }
     if ((dirxM == 1) && (diry == 0))   {
         if (div2 == 3)  {
             ROC = 9;
         }
         if (div2 == 2)  {
             ROC = 10;
         }
         if (div2 == 1)  {
             ROC = 11;
         }
         if (div2 == 0)  {
             ROC = 12;
         }
     }
     if ((dirxM == 1) && (diry == 1))   {
         if (div2 == 0)  {
             ROC = 13;
         }
         if (div2 == -1)  {
             ROC = 14;
         }
         if (div2 == -2)  {
             ROC = 15;
         }
         if (div2 == -3)  {
             ROC = 16;
         }
     }
     if (dirxM == 0)  {
	 if ((divp1 == 0) || (divp1 == 1))  {
	    col = 80;
         }
         if ((divp1 > 1) && (divp1 <= 80))  {
	    col = 80 - divp1 + 1;
         } 
     }
     if (dirxM == 1)  {
	 if ((divp1 == -79) || (divp1 == -80))  {
	    col = 80;
         }
         if ((divp1 > -79) && (divp1 <= 0))  {
	    col = -divp1 + 1;
         } 
     }

     if ((ROC == 1) || (ROC == 9))  {
     	divp2 = (yPix-3*ROCvert)/pixelY;
     }
     if ((ROC == 2) || (ROC == 10))  {
     	divp2 = (yPix-2*ROCvert)/pixelY;
     }
     if ((ROC == 3) || (ROC == 11))  {
     	divp2 = (yPix-1*ROCvert)/pixelY;
     }
     if ((ROC == 4) || (ROC == 12))  {
     	divp2 = (yPix-0*ROCvert)/pixelY;
     }

     if ((ROC == 5) || (ROC == 13))  {
     	divp2 = (yPix+0*ROCvert)/pixelY;
     }
     if ((ROC == 6) || (ROC == 14))  {
     	divp2 = (yPix+1*ROCvert)/pixelY;
     }
     if ((ROC == 7) || (ROC == 15))  {
     	divp2 = (yPix+2*ROCvert)/pixelY;
     }
     if ((ROC == 8) || (ROC == 16))  {
     	divp2 = (yPix+3*ROCvert)/pixelY;
     }

     if ((ROC == 1) || (ROC == 2) || (ROC == 3) || (ROC == 4) || (ROC == 9) || (ROC == 10) || (ROC == 11) || (ROC == 12))  {
	if ((divp2 == 0) || (divp2 == 1))  {
	   row = 52;
	}
	if ((divp2 == 52) || (divp2 == 53))  {
	   row = 1;
	}
	if ((divp2 < 52) && (divp2 > 1))  {
	   row = 52 - divp2 + 1;
	}
     }
     if ((ROC == 5) || (ROC == 6) || (ROC == 7) || (ROC == 8) || (ROC == 13) || (ROC == 14) || (ROC == 15) || (ROC == 16))  {
	if ((divp2 == 0) || (divp2 == -1))  {
	   row = 1;
	}
	if ((divp2 == -52) || (divp2 == -53))  {
	   row = 52;
	}
	if ((divp2 > -52) && (divp2 < -1))  {
	   row = -divp2;
	}
     }

     if ((mod>=1) && (mod<17) && (ROC>=1) && (ROC<17) && (row>=1) && (row<53) && (col>=1) && (col<81))   {
	   fEventaction->AddEnergyPixel(edep,mod,ROC,row,col);
	   /*G4cout
	      << "\n Enery added to pixel = "
	      << G4BestUnit(edep, "Energy")
	      << "\n module = "
	      << mod
	      << G4endl;*/ 
     }   
   }

   /*G4cout
     << "\n edep = "
     << G4BestUnit(edep, "Energy")
     << "\n module = "
     << mod
     << "\n ROC = "
     << ROC
     << "\n row = "
     << row
     << "\n col = "
     << col
     << "\n diry = "
     << diry
     << "\n dirxM = "
     << dirxM
     << "\n divp1 = "
     << divp1
     << G4endl;*/ 
 }

 if (volume == fDetectorconstruction->GetBPIX1314()) {
   fEventaction->AddEnergyDepositL7(edep);

   mod = -1;
   ROC = -1;
   row = -1;
   col = -1;
   dirx = -1;
   dirxM = -1;
   divp1 = -1;

   zPixTr = z_after - poszBPIX1314;
   xPix1 = x_after*cos(YangleBPIX) + zPixTr*sin(YangleBPIX);
   yPix1 = y_after;
   zPix1 = -x_after*sin(YangleBPIX) + zPixTr*cos(YangleBPIX);
   xPix = xPix1;
   yPix = yPix1*cos(XangleBPIX) - zPix1*sin(XangleBPIX);
   zPix = yPix1*sin(XangleBPIX) + zPix1*cos(XangleBPIX);

   if ((zPix <= (BPIXSizeZ/2))  &&  (zPix >= (BPIXSizeZ/2 - pixelDepth))) {
     /*G4cout
       << "\n zPix1314 = "
       << G4BestUnit(zPix, "Length")
       << G4endl; */
     if (xPix == 0)   {
	iSecret = rand() % 99;
     	if (iSecret < 50)   {
	   dirx = 0;
     	}
     	if (iSecret >= 50)   {
   	   dirx = 1;
     	}
     }
     if (xPix > 0)    {
	dirx = 0;
     }
     if (xPix < 0)    {
        dirx = 1;
     }
     if (dirx == 0)  {
        mod = 13; //module ID
	xPixM = xPix - BPIXSizeX/2;
	/*G4cout
	   << "\n xPixM = "
	   << G4BestUnit(xPixM, "Length")
	   << G4endl; */
        if (xPixM == 0)   {
	  iSecret = rand() % 99;
     	  if (iSecret < 50)   {
	     dirxM = 0;
     	  }
     	  if (iSecret >= 50)   {
   	     dirxM = 1;
     	  }
        }
     	if (xPixM > 0)    {
	  dirxM = 0;
     	}
     	if (xPixM < 0)    {
          dirxM = 1;
        }
     }
     if (dirx == 1)  {
        mod = 14; //module ID
	xPixM = xPix + BPIXSizeX/2;
	/*G4cout
	   << "\n xPixM = "
	   << G4BestUnit(xPixM, "Length")
	   << G4endl; */
        if (xPixM == 0)   {
	  iSecret = rand() % 99;
     	  if (iSecret < 50)   {
	     dirxM = 0;
     	  }
     	  if (iSecret >= 50)   {
   	     dirxM = 1;
     	  }
        }
     	if (xPixM > 0)    {
	  dirxM = 0;
     	}
     	if (xPixM < 0)    {
          dirxM = 1;
        }
     }
     div1 = xPix/ROChor;
     div2 = yPix/ROCvert;
     divp1 = xPixM/pixelX;
     ROC = 0;
     if (yPix == 0)   {
	iSecret = rand() % 99;
     	if (iSecret < 50)   {
	   diry = 0;
     	}
     	if (iSecret >= 50)   {
   	   diry = 1;
     	}
     }
     if (yPix > 0)    {
	diry = 0;
     }
     if (yPix < 0)    {
        diry = 1;
     }
     if ((dirxM == 0) && (diry == 0))   {
         if (div2 == 3)  {
             ROC = 1;
         }
         if (div2 == 2)  {
             ROC = 2;
         }
         if (div2 == 1)  {
             ROC = 3;
         }
         if (div2 == 0)  {
             ROC = 4;
         }
     }
     if ((dirxM == 0) && (diry == 1))   {
         if (div2 == 0)  {
             ROC = 5;
         }
         if (div2 == -1)  {
             ROC = 6;
         }  
         if (div2 == -2)  {
             ROC = 7;
         }
         if (div2 == -3)  {
             ROC = 8;
         }
     }
     if ((dirxM == 1) && (diry == 0))   {
         if (div2 == 3)  {
             ROC = 9;
         }
         if (div2 == 2)  {
             ROC = 10;
         }
         if (div2 == 1)  {
             ROC = 11;
         }
         if (div2 == 0)  {
             ROC = 12;
         }
     }
     if ((dirxM == 1) && (diry == 1))   {
         if (div2 == 0)  {
             ROC = 13;
         }
         if (div2 == -1)  {
             ROC = 14;
         }
         if (div2 == -2)  {
             ROC = 15;
         }
         if (div2 == -3)  {
             ROC = 16;
         }
     }
     if (dirxM == 0)  {
	 if ((divp1 == 0) || (divp1 == 1))  {
	    col = 80;
         }
         if ((divp1 > 1) && (divp1 <= 80))  {
	    col = 80 - divp1 + 1;
         } 
     }
     if (dirxM == 1)  {
	 if ((divp1 == -79) || (divp1 == -80))  {
	    col = 80;
         }
         if ((divp1 > -79) && (divp1 <= 0))  {
	    col = -divp1 + 1;
         } 
     }

     if ((ROC == 1) || (ROC == 9))  {
     	divp2 = (yPix-3*ROCvert)/pixelY;
     }
     if ((ROC == 2) || (ROC == 10))  {
     	divp2 = (yPix-2*ROCvert)/pixelY;
     }
     if ((ROC == 3) || (ROC == 11))  {
     	divp2 = (yPix-1*ROCvert)/pixelY;
     }
     if ((ROC == 4) || (ROC == 12))  {
     	divp2 = (yPix-0*ROCvert)/pixelY;
     }

     if ((ROC == 5) || (ROC == 13))  {
     	divp2 = (yPix+0*ROCvert)/pixelY;
     }
     if ((ROC == 6) || (ROC == 14))  {
     	divp2 = (yPix+1*ROCvert)/pixelY;
     }
     if ((ROC == 7) || (ROC == 15))  {
     	divp2 = (yPix+2*ROCvert)/pixelY;
     }
     if ((ROC == 8) || (ROC == 16))  {
     	divp2 = (yPix+3*ROCvert)/pixelY;
     }

     if ((ROC == 1) || (ROC == 2) || (ROC == 3) || (ROC == 4) || (ROC == 9) || (ROC == 10) || (ROC == 11) || (ROC == 12))  {
	if ((divp2 == 0) || (divp2 == 1))  {
	   row = 52;
	}
	if ((divp2 == 52) || (divp2 == 53))  {
	   row = 1;
	}
	if ((divp2 < 52) && (divp2 > 1))  {
	   row = 52 - divp2 + 1;
	}
     }
     if ((ROC == 5) || (ROC == 6) || (ROC == 7) || (ROC == 8) || (ROC == 13) || (ROC == 14) || (ROC == 15) || (ROC == 16))  {
	if ((divp2 == 0) || (divp2 == -1))  {
	   row = 1;
	}
	if ((divp2 == -52) || (divp2 == -53))  {
	   row = 52;
	}
	if ((divp2 > -52) && (divp2 < -1))  {
	   row = -divp2;
	}
     }

     if ((mod>=1) && (mod<17) && (ROC>=1) && (ROC<17) && (row>=1) && (row<53) && (col>=1) && (col<81))   {
	   fEventaction->AddEnergyPixel(edep,mod,ROC,row,col);
	   /*G4cout
	      << "\n Enery added to pixel = "
	      << G4BestUnit(edep, "Energy")
	      << "\n module = "
	      << mod
	      << G4endl;*/ 
     }   
   }

   /*G4cout
     << "\n edep = "
     << G4BestUnit(edep, "Energy")
     << "\n module = "
     << mod
     << "\n ROC = "
     << ROC
     << "\n row = "
     << row
     << "\n col = "
     << col
     << "\n diry = "
     << diry
     << "\n dirxM = "
     << dirxM
     << "\n divp1 = "
     << divp1
     << G4endl;*/ 
 }

 if (volume == fDetectorconstruction->GetBPIX1516()) {
   fEventaction->AddEnergyDepositL8(edep);

   mod = -1;
   ROC = -1;
   row = -1;
   col = -1;
   dirx = -1;
   dirxM = -1;
   divp1 = -1;

   zPixTr = z_after - poszBPIX1516;
   xPix1 = x_after*cos(YangleBPIX) + zPixTr*sin(YangleBPIX);
   yPix1 = y_after;
   zPix1 = -x_after*sin(YangleBPIX) + zPixTr*cos(YangleBPIX);
   xPix = xPix1;
   yPix = yPix1*cos(XangleBPIX) - zPix1*sin(XangleBPIX);
   zPix = yPix1*sin(XangleBPIX) + zPix1*cos(XangleBPIX);

   if ((zPix <= (BPIXSizeZ/2))  &&  (zPix >= (BPIXSizeZ/2 - pixelDepth))) {
     /*G4cout
       << "\n zPix1516 = "
       << G4BestUnit(zPix, "Length")
       << G4endl; */
     if (xPix == 0)   {
	iSecret = rand() % 99;
     	if (iSecret < 50)   {
	   dirx = 0;
     	}
     	if (iSecret >= 50)   {
   	   dirx = 1;
     	}
     }
     if (xPix > 0)    {
	dirx = 0;
     }
     if (xPix < 0)    {
        dirx = 1;
     }
     if (dirx == 0)  {
        mod = 15; //module ID
	xPixM = xPix - BPIXSizeX/2;
	/*G4cout
	   << "\n xPixM = "
	   << G4BestUnit(xPixM, "Length")
	   << G4endl; */
        if (xPixM == 0)   {
	  iSecret = rand() % 99;
     	  if (iSecret < 50)   {
	     dirxM = 0;
     	  }
     	  if (iSecret >= 50)   {
   	     dirxM = 1;
     	  }
        }
     	if (xPixM > 0)    {
	  dirxM = 0;
     	}
     	if (xPixM < 0)    {
          dirxM = 1;
        }
     }
     if (dirx == 1)  {
        mod = 16; //module ID
	xPixM = xPix + BPIXSizeX/2;
	/*G4cout
	   << "\n xPixM = "
	   << G4BestUnit(xPixM, "Length")
	   << G4endl; */
        if (xPixM == 0)   {
	  iSecret = rand() % 99;
     	  if (iSecret < 50)   {
	     dirxM = 0;
     	  }
     	  if (iSecret >= 50)   {
   	     dirxM = 1;
     	  }
        }
     	if (xPixM > 0)    {
	  dirxM = 0;
     	}
     	if (xPixM < 0)    {
          dirxM = 1;
        }
     }
     div1 = xPix/ROChor;
     div2 = yPix/ROCvert;
     divp1 = xPixM/pixelX;
     ROC = 0;
     if (yPix == 0)   {
	iSecret = rand() % 99;
     	if (iSecret < 50)   {
	   diry = 0;
     	}
     	if (iSecret >= 50)   {
   	   diry = 1;
     	}
     }
     if (yPix > 0)    {
	diry = 0;
     }
     if (yPix < 0)    {
        diry = 1;
     }
     if ((dirxM == 0) && (diry == 0))   {
         if (div2 == 3)  {
             ROC = 1;
         }
         if (div2 == 2)  {
             ROC = 2;
         }
         if (div2 == 1)  {
             ROC = 3;
         }
         if (div2 == 0)  {
             ROC = 4;
         }
     }
     if ((dirxM == 0) && (diry == 1))   {
         if (div2 == 0)  {
             ROC = 5;
         }
         if (div2 == -1)  {
             ROC = 6;
         }  
         if (div2 == -2)  {
             ROC = 7;
         }
         if (div2 == -3)  {
             ROC = 8;
         }
     }
     if ((dirxM == 1) && (diry == 0))   {
         if (div2 == 3)  {
             ROC = 9;
         }
         if (div2 == 2)  {
             ROC = 10;
         }
         if (div2 == 1)  {
             ROC = 11;
         }
         if (div2 == 0)  {
             ROC = 12;
         }
     }
     if ((dirxM == 1) && (diry == 1))   {
         if (div2 == 0)  {
             ROC = 13;
         }
         if (div2 == -1)  {
             ROC = 14;
         }
         if (div2 == -2)  {
             ROC = 15;
         }
         if (div2 == -3)  {
             ROC = 16;
         }
     }
     if (dirxM == 0)  {
	 if ((divp1 == 0) || (divp1 == 1))  {
	    col = 80;
         }
         if ((divp1 > 1) && (divp1 <= 80))  {
	    col = 80 - divp1 + 1;
         } 
     }
     if (dirxM == 1)  {
	 if ((divp1 == -79) || (divp1 == -80))  {
	    col = 80;
         }
         if ((divp1 > -79) && (divp1 <= 0))  {
	    col = -divp1 + 1;
         } 
     }

     if ((ROC == 1) || (ROC == 9))  {
     	divp2 = (yPix-3*ROCvert)/pixelY;
     }
     if ((ROC == 2) || (ROC == 10))  {
     	divp2 = (yPix-2*ROCvert)/pixelY;
     }
     if ((ROC == 3) || (ROC == 11))  {
     	divp2 = (yPix-1*ROCvert)/pixelY;
     }
     if ((ROC == 4) || (ROC == 12))  {
     	divp2 = (yPix-0*ROCvert)/pixelY;
     }

     if ((ROC == 5) || (ROC == 13))  {
     	divp2 = (yPix+0*ROCvert)/pixelY;
     }
     if ((ROC == 6) || (ROC == 14))  {
     	divp2 = (yPix+1*ROCvert)/pixelY;
     }
     if ((ROC == 7) || (ROC == 15))  {
     	divp2 = (yPix+2*ROCvert)/pixelY;
     }
     if ((ROC == 8) || (ROC == 16))  {
     	divp2 = (yPix+3*ROCvert)/pixelY;
     }

     if ((ROC == 1) || (ROC == 2) || (ROC == 3) || (ROC == 4) || (ROC == 9) || (ROC == 10) || (ROC == 11) || (ROC == 12))  {
	if ((divp2 == 0) || (divp2 == 1))  {
	   row = 52;
	}
	if ((divp2 == 52) || (divp2 == 53))  {
	   row = 1;
	}
	if ((divp2 < 52) && (divp2 > 1))  {
	   row = 52 - divp2 + 1;
	}
     }
     if ((ROC == 5) || (ROC == 6) || (ROC == 7) || (ROC == 8) || (ROC == 13) || (ROC == 14) || (ROC == 15) || (ROC == 16))  {
	if ((divp2 == 0) || (divp2 == -1))  {
	   row = 1;
	}
	if ((divp2 == -52) || (divp2 == -53))  {
	   row = 52;
	}
	if ((divp2 > -52) && (divp2 < -1))  {
	   row = -divp2;
	}
     }

     if ((mod>=1) && (mod<17) && (ROC>=1) && (ROC<17) && (row>=1) && (row<53) && (col>=1) && (col<81))   {
	   fEventaction->AddEnergyPixel(edep,mod,ROC,row,col);
	   /*G4cout
	      << "\n Enery added to pixel = "
	      << G4BestUnit(edep, "Energy")
	      << "\n module = "
	      << mod
	      << G4endl;*/ 
     }   
   }

   /*G4cout
     << "\n edep = "
     << G4BestUnit(edep, "Energy")
     << "\n module = "
     << mod
     << "\n ROC = "
     << ROC
     << "\n row = "
     << row
     << "\n col = "
     << col
     << "\n diry = "
     << diry
     << "\n dirxM = "
     << dirxM
     << "\n divp1 = "
     << divp1
     << G4endl;*/ 
 }



 if (volume == fDetectorconstruction->GetDet1()) {
   fEventaction->AddEnergyDeposit1 (edep);

   if ((z >= -Dist/2 - Det1SizeZ) && (z <= -Dist/2 - Det1SizeZ + Strip1Depth) && (x >= (-Strip1Length/2)) && (x <= Strip1Length/2)) {
     div = y/(StripWidth+StripDist);
     if (y == 0)   {
	iSecret = rand() % 99;
     	if (iSecret < 50)   {
	   stripNo = totalNbStrips/2 + div;
     	}
     	if (iSecret >= 50)   {
   	   stripNo = totalNbStrips/2 + 1 + div;
     	}
     }
     if (y > 0)    {
	stripNo = totalNbStrips/2 + 1 + div;
     }
     if (y < 0)    {
        stripNo = totalNbStrips/2 + div;
     }
     if ((stripNo <= totalNbStrips) && (stripNo > 0))  {
     	//G4cout
	   //<< "\n Continuous energy deposition through strip No "
           //<< stripNo
           //<< "\n x = "
	   //<< x
           //<< "\n div = "
	   //<< div
           //<< G4endl;
        if (x > 0)  {
	   fEventaction->AddEnergyStrip1a(edep,stripNo);
        }
        if (x < 0)  {
	   fEventaction->AddEnergyStrip1b(edep,stripNo);
        }
        if (x == 0)  {
	    jSecret = rand() % 99;
	    if (jSecret < 50)   {
	       fEventaction->AddEnergyStrip1a(edep,stripNo);
	    }
	    if (jSecret >= 50)   {
	       fEventaction->AddEnergyStrip1b(edep,stripNo);
	    }
        }
     }
   }
 }

 if (volume == fDetectorconstruction->GetDet2()) {
   fEventaction->AddEnergyDeposit2 (edep);

   if ((z >= Dist/2) && (z <= (Dist/2 + Strip2Depth))&& (x >= (-Strip2Length/2)) && (x <= Strip2Length/2)) {
     div = y/(StripWidth+StripDist);
     if (y == 0)   {
	iSecret = rand() % 99;
     	if (iSecret < 50)   {
	   stripNo = totalNbStrips/2 + div;
     	}
     	if (iSecret >= 50)   {
   	   stripNo = totalNbStrips/2 + 1 + div;  //momentum direction when entering Pixel Detector 1 and exiting Pixel Detector 2
     	}
     }
     if (y > 0)    {
	stripNo = totalNbStrips/2 + 1 + div;
     }
     if (y < 0)    {
        stripNo = totalNbStrips/2 + div;
     }
     if ((stripNo <= totalNbStrips)&& (stripNo > 0))   {
     	//G4cout
	   //<< "\n Continuous energy deposition through strip No "
           //<< stripNo
           //<< "\n x = "
	   //<< x
           //<< "\n div = "
	   //<< div
           //<< G4endl;
        if (x > 0)  {
	   fEventaction->AddEnergyStrip2a(edep,stripNo);
        }
        if (x < 0)  {
	   fEventaction->AddEnergyStrip2b(edep,stripNo);
        }
        if (x == 0)  {
	    jSecret = rand() % 99;
	    if (jSecret < 50)   {
	       fEventaction->AddEnergyStrip2a(edep,stripNo);
	    }
	    if (jSecret >= 50)   {
	       fEventaction->AddEnergyStrip2b(edep,stripNo);
	    }
        }
     }
   }
 }
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

