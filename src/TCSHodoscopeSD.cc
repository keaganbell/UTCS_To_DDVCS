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

#include "TCSHodoscopeSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "TCSTrackInformation.hh"
#include "TCSHodoscopeConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSHodoscopeSD::TCSHodoscopeSD(const G4String& name,
				   const G4String& hitsCollectionName) 
 : G4VSensitiveDetector(name),
   fHitsCollection(NULL)
   ////   fHitsCollection(0), lastID(-1)
{
  collectionName.insert(hitsCollectionName);

  G4cout << "TCSHodoscopeSD::TCSHodoscopeSD: constructed" << G4endl;
  G4cout << "   name = " << name << G4endl;
  G4cout << "   hitsCollectionName = " << hitsCollectionName << G4endl;
  //  getchar();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSHodoscopeSD::~TCSHodoscopeSD() 
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSHodoscopeSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  //  G4cout << "TCSHodoscopeSD::Initialize: creating hit collection. "
  //  << G4endl;
  //  getchar();

  fHitsCollection 
  = new TCSHodoscopeHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce

  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 

  //  G4cout << "TCSHodoscopeSD::Initialize: initialized" << G4endl;
  //  getchar();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TCSHodoscopeSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  //  G4cout << "TCSHodoscopeSD::ProcessHits: entered" << G4endl;
  
  //  G4int pid = step->GetTrack()->GetDefinition()->GetPDGEncoding();
  TCSTrackInformation* info =
    (TCSTrackInformation*)(step->GetTrack()->GetUserInformation());
  G4int pid = info->GetOriginalParticle()->GetPDGEncoding();

  G4ThreeVector pos = step->GetTrack()->GetPosition();

  // Particle in the crystals, save energy deposit for dose calc-s.

  if (step->GetTrack()->GetVolume()->GetName() == "hModule_PV") {

    G4double edep = step->GetTotalEnergyDeposit();
    if (edep > 0.) {

      TCSHodoscopeConstruction cc;
      int ncol = cc.GetNCOL();
      int nrow = cc.GetNROW();
      ///int nmod=step->GetPostStepPoint()->GetTouchableHandle()->GetCopyNumber(0);
      int nmod=step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(0);
      int row = nmod/ncol;
      int col = nmod - row*ncol;

      if (row < 0 || row >= nrow)
	G4cout << "*** TCSHodoscopeSD::ProcessHits: wrong row = " << row
	       << G4endl;
      if (col < 0 || col >= ncol)
	G4cout << "*** TCSHodoscopeSD::ProcessHits: wrong col = " << col
	       << G4endl;

      TCSHodoscopeHit* hit_d = new TCSHodoscopeHit(col,row,pid,edep,pos, 0);
      fHitsCollection->insert( hit_d );
      //  hit_d->Print();
      //  getchar();
    }

    return true;
  }

  // Particle entering calorimeter, save kinetic energy for flux calc-s.
  /*
  if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==
      "Module_phys" &&
      step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==
      "Hodoscope_PV" &&
      step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {

    TCSHodoscopeConstruction cc;
    int ncol = cc.GetNCOL();
    int nrow = cc.GetNROW();
    int nmod=step->GetPostStepPoint()->GetTouchableHandle()->GetCopyNumber(0);
    int row = nmod/ncol;
    int col = nmod - row*ncol;

    if (row < 0 || row >= nrow)
      G4cout << "*** TCSHodoscopeSD::ProcessHits: wrong row = " << row
	     << G4endl;
    if (col < 0 || col >= ncol)
      G4cout << "*** TCSHodoscopeSD::ProcessHits: wrong col = " << col
	     << G4endl;

    G4double ekin = step->GetTrack()->GetKineticEnergy();
    TCSHodoscopeHit* hit_f=new TCSHodoscopeHit(col,row,pid, ekin, pos, 1);
    fHitsCollection->insert( hit_f );
    return true;
  }
  */

  /*
  if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==
  "Module_phys" &&
  step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
    G4cout << "At Module_phys boundary, PreStep Vol. = " << step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() << G4endl;
    return true;
  }
  */

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSHodoscopeSD::EndOfEvent(G4HCofThisEvent*)
{
  //  G4cout << "TCSHodoscopeSD::EndOfEvent:" << G4endl;
  //  G4cout << "  Number of entries in hit collection = "
  //	 << fHitsCollection->entries() << G4endl;

  ////  for ( G4int i=0; i<fHitsCollection->entries(); i++ )
  ////    (*fHitsCollection)[i]->delete();

  //  getchar();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
