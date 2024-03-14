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

#include "TCSSteppingAction.hh"
#include "TCSEventAction.hh"
#include "TCSDetectorConstruction.hh"

#include "globals.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcommand.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSSteppingAction::TCSSteppingAction(TCSEventAction* eventAction)
: G4UserSteppingAction(), fEventAction(eventAction), fScoringVolume(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSSteppingAction::~TCSSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSSteppingAction::UserSteppingAction(const G4Step* step)
{
  // horrible hack: skip if particle is not a neutron
  // (attempt to get more neutrons for the event generator test!)
  /// if (step->GetTrack()->GetDefinition()->GetPDGEncoding()!=2112)     return;
  ////

  /*
  G4cout << "TCSSteppingAction::UserSteppingAction" << G4endl;
  G4cout << " PreStepPoint volume:" <<
    step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName()
	 << G4endl;
  G4cout << " PostStepPoint volume:" <<
    step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName()
	 << G4endl;
  if (step->GetPreStepPoint()->GetStepStatus() == fGeomBoundary)
    G4cout << " PreStepPoint status: on Boundary !" << G4endl;
  if (step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary)
    G4cout << " PostStepPoint status: on Boundary !" << G4endl;
  getchar();
  */

  if (step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==
      "TargetAssembly_PV" &&
      step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==
      "ChamberAssembly_PV" &&
      step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
    G4int pid = step->GetTrack()->GetDefinition()->GetPDGEncoding();
    switch (pid) {
    case  11 : fEventAction->AddNTar_eneg(1);
      break;
    case -11 : fEventAction->AddNTar_epos(1);
      break;
    default: ;
    }
  }

  G4LogicalVolume* volume =
 step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

  if (!fScoringVolume) { 
    const TCSDetectorConstruction* detectorConstruction
      = static_cast<const TCSDetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detectorConstruction->GetScoringVolume();   
  }

  // Check if we are in scoring volume.
  if (volume != fScoringVolume) return;

  // Collect energy deposited in this step.
  G4double edepStep = step->GetTotalEnergyDeposit();
  fEventAction->AddEdep(edepStep);  

  //G4int ncol = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();
  //G4int nrow =step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1);
  //G4int nh = step->GetPreStepPoint()->GetTouchableHandle()->GetHistoryDepth();
  //  cout << "TCSSteppingAction::UserSteppingAction: ncol = " << ncol
  //       << "  nrow = " << nrow << "  Edep = " << edepStep << G4endl;
  //  cout << "  history depth = " << nh << G4endl;
  //  getchar();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
