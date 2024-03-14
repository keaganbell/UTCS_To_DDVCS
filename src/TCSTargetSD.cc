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

#include "TCSTargetSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "TCSTrackInformation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSTargetSD::TCSTargetSD(const G4String& name,
			 const G4String& hitsCollectionName) 
 : G4VSensitiveDetector(name),
   fHitsCollection(NULL)
{
  collectionName.insert(hitsCollectionName);

  G4cout << "TCSTargetSD::TCSTargetSD: constructed" << G4endl;
  G4cout << "   name = " << name << G4endl;
  G4cout << "   hitsCollectionName = " << hitsCollectionName << G4endl;
  //  getchar();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSTargetSD::~TCSTargetSD() 
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSTargetSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  //  G4cout <<"TCSTargetSD::Initialize: creating hit collection.."<< G4endl;
  //  getchar();

  fHitsCollection 
    = new TCSTargetHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce

  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 

  //  G4cout << "TCSTargetSD::Initialize: initialized" << G4endl;
  //  getchar();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TCSTargetSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  //  G4cout << "TCSTargetSD::ProcessHits: entered" << G4endl;
  
  // Particle leaving volume, save kinetic energy for flux calc-s.
  if (step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==
      "TargetAssembly_PV" &&
      step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==
      "ChamberAssembly_PV" &&
      step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
    G4int pid = step->GetTrack()->GetDefinition()->GetPDGEncoding();
    G4ThreeVector pos = step->GetTrack()->GetPosition();
    G4double ekin = step->GetTrack()->GetKineticEnergy();
    TCSTargetHit* hit = new TCSTargetHit(pid, ekin, pos, 1);
    fHitsCollection->insert( hit );
    ////    return true;
  }

  //  TCSTrackInformation* info =
  //    (TCSTrackInformation*)(step->GetTrack()->GetUserInformation());
  //  G4cout << " Original Track ID " << info->GetOriginalTrackID() << G4endl;
  //  G4cout << " Original particle "
  //	 << info->GetOriginalParticle()->GetPDGEncoding() << G4endl;
  //  getchar();

  //  G4cout << "TCSTargetSD::ProcessHits:" << G4endl;
  //  G4cout << " PreStepPoint volume:" <<
  //    step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName()
  //  	 << G4endl;
  //  G4cout << " PostStepPoint volume:" <<
  //    step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName()
  //  	 << G4endl;
  //  getchar();

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSTargetSD::EndOfEvent(G4HCofThisEvent*)
{
  //  if ( verboseLevel>1 ) { 
  //    G4int nofHits = fHitsCollection->entries();
  //    G4cout << G4endl
  //	   << "-------->Hits Collection: in this event there are " << nofHits 
  //	   << " hits in the tracker chambers: " << G4endl;
  //for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->delete();
  //  }

  //  G4cout << "TCSTargetSD::EndOfEvent: end of event" << G4endl;
  //  getchar();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
