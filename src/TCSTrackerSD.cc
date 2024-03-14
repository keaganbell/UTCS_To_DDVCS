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

#include "TCSTrackerSD.hh"
#include "TCSEventAction.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include <string>

#include "TCSTrackInformation.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSTrackerSD::TCSTrackerSD(const G4String& name,
			   const G4String& hitsCollectionName) 
 : G4VSensitiveDetector(name),
   fHitsCollection(0), lastID(-1)
{
  collectionName.insert(hitsCollectionName);

  G4cout << "TCSTrackerXSD::TCSTrackerXSD: constructed" << G4endl;
  G4cout << "   name = " << name << G4endl;
  G4cout << "   hitsCollectionName = " << hitsCollectionName << G4endl;
  //  getchar();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSTrackerSD::~TCSTrackerSD() 
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSTrackerSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection

  fHitsCollection 
  = new TCSTrackerHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce

  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 

  //  G4cout << "TCSTrackerSD::Initialize: initialized" << G4endl;
  // G4cout << "   SensitiveDetectorName = " << SensitiveDetectorName << G4endl;
  //  G4cout << "   collectionName = " << collectionName[0] << G4endl;
  //  getchar();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TCSTrackerSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  TCSTrackInformation* info =
    (TCSTrackInformation*)(step->GetTrack()->GetUserInformation());

  string vname =
    step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName();

  // Particle in the tracker's drift gas volume.
  ///  if (vname == "Drift_PV") {
  if (vname.substr(0,8) == "Drift_PV") {

    G4int layer = -1;

    //    string WorldName =
      //step->GetPreStepPoint()->GetTouchableHandle()->GetVolume(2)->GetName();
    //if (WorldName.substr(0,7)== "tracker" && WorldName.substr(8,5) == "World")
    //layer = atoi(WorldName.substr(7,1).c_str()) - 1;

    string detName =
      step->GetPreStepPoint()->GetTouchableHandle()->GetVolume(1)->GetName();
    if (detName.substr(0,7)== "Tracker")
      layer = atoi(detName.substr(7,1).c_str()) - 1;
    else
      G4cout << "*** TCSTrackerSD::ProcessHits: wrong detName = "
	     << detName << " ***" << G4endl;

    //    G4ThreeVector trackPos = step->GetTrack()->GetPosition();
    G4ThreeVector prestepPos = step->GetPreStepPoint()->GetPosition();
    G4ThreeVector poststepPos = step->GetPostStepPoint()->GetPosition();

    //G4ThreeVector loctrackPos =step->GetPreStepPoint()->GetTouchableHandle()->
    //      GetHistory()->GetTopTransform().TransformPoint(trackPos);
    G4ThreeVector locprePos = step->GetPreStepPoint()->GetTouchableHandle()->
      GetHistory()->GetTopTransform().TransformPoint(prestepPos);
    G4ThreeVector locpostPos = step->GetPreStepPoint()->GetTouchableHandle()->
      GetHistory()->GetTopTransform().TransformPoint(poststepPos);

    G4double totedep = step->GetTotalEnergyDeposit();
    G4double nonionedep = step->GetNonIonizingEnergyDeposit();

    G4double length = step->GetStepLength();

    G4int trackID = step->GetTrack()->GetTrackID();
    //    G4int stepNumber = step->GetTrack()->GetCurrentStepNumber();

    G4int origPID = info->GetOriginalParticle()->GetPDGEncoding();

    G4int PID = step->GetTrack()->GetDefinition()->GetPDGEncoding();

    //    G4ThreeVector mom = step->GetTrack()->GetMomentum();
    //    G4double P = sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]);
    //    TCSTrackerHit* hit = new TCSTrackerHit(pid, P, localPos, layer);
    //    fHitsCollection->insert( hit );

    G4int quarter = TCSEventAction::GetQuarter(prestepPos.getX(),
					       prestepPos.getY());

    TCSTrackerHit* hit = new TCSTrackerHit(
				   (locprePos.getX()+locpostPos.getX())/2,
				   (locprePos.getY()+locpostPos.getY())/2,
				   totedep - nonionedep, length,
				   quarter, layer,
				   PID, origPID, trackID);

    fHitsCollection->insert( hit );

    /*
    G4cout << "TCSTrackerSD::ProcessHits: vol. name = " << vname << G4endl;
    G4cout << "  pid = " << PID << "  origpid = " << origPID << G4endl;
    G4cout << "  locprePos   = " << locprePos.getX() << " " 
	                         << locprePos.getY() << " "
	                         << locprePos.getZ() << " " << G4endl;
    G4cout << "  locpostPos  = " << locpostPos.getX() << " " 
	                         << locpostPos.getY() << " "
	                         << locpostPos.getZ() << " " << G4endl;
    G4cout << "  Edep total  = " << totedep << G4endl;
    G4cout << "  Edep nonion = " << nonionedep << G4endl;
    G4cout << "  length      = " << length << G4endl;
    */
    /*
    G4cout << "TCSTrackerSD::ProcessHits: vol. name = " << vname << G4endl;
    G4cout << "  pid = " << pid << "  P = " << P << G4endl;
    G4cout << "  prestepPos   = " << prestepPos.getX() << " " 
	                          << prestepPos.getY() << " "
	                          << prestepPos.getZ() << " " << G4endl;
    G4cout << "  locprePos   = " << locprePos.getX() << " " 
	                         << locprePos.getY() << " "
	                         << locprePos.getZ() << " " << G4endl;
    G4cout << "  locpostPos  = " << locpostPos.getX() << " " 
	                         << locpostPos.getY() << " "
	                         << locpostPos.getZ() << " " << G4endl;
    G4cout << "  loctrackPos = " << loctrackPos.getX() << " " 
	                         << loctrackPos.getY() << " "
	                         << loctrackPos.getZ() << " " << G4endl;
    G4cout << "  Edep total  = " << totedep << G4endl;
    G4cout << "  Edep nonion = " << nonionedep << G4endl;
    G4cout << "  track ID    = " << trackID << G4endl;
    G4cout << "  step number = " << stepNumber << G4endl;
    G4cout << "  origPID     = " << origPID << G4endl;
    G4cout << "  volume 1    = " <<
      step->GetPreStepPoint()->GetTouchableHandle()->GetVolume(1)->GetName()
	   << G4endl;
    G4cout << "  volume 2    = " <<
      step->GetPreStepPoint()->GetTouchableHandle()->GetVolume(2)->GetName()
	   << G4endl;
    G4cout << "  WorldName   = " << WorldName << G4endl;
    G4cout << "  layer       = " << layer << G4endl;
    G4cout << "  quarter     = " << quarter << G4endl;
    getchar();
    */

    ////    return true;
  }

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSTrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  //  if ( verboseLevel>1 ) { 
  //  G4int nofHits = fHitsCollection->entries();
  //   G4cout << G4endl
  //          << "-------->Hits Collection: in this event they are " << nofHits 
  //          << " hits in the tracker chambers: " << G4endl;
  //  for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->delete();
  //}

  //  G4cout << "TCSTrackerXSD::EndOfEvent: end of event" << G4endl;
  //  getchar();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
