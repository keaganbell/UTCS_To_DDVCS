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
// $Id: B3StackingAction.cc 66536 2012-12-19 14:32:36Z ihrivnac $
// 
/// \file B3StackingAction.cc
/// \brief Implementation of the B3StackingAction class

#include "TCSStackingAction.hh"

#include "G4Track.hh"
//#include "G4NeutrinoE.hh"
#include <fstream>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSStackingAction::TCSStackingAction()
{
  ifstream file("stacking_control.txt"); // Open the file for reading.

  string line;
  istringstream iss;

  getline(file, line);  iss.str(line);
  iss >> fPrimaryOnly;
  //  getline(file, line);  iss.str(line);
  //  iss >> fAcceptanceOnly;

  file.close();

  G4cout << "TCSStackingAction::TCSStackingAction:" << G4endl;

  G4cout << "  fPrimaryOnly = " << fPrimaryOnly;
  if (fPrimaryOnly)
    G4cout << " ==> only primary particle tracking!" << G4endl;
  else
    G4cout << G4endl;

  //  G4cout << "  fAcceptanceOnly = " << fAcceptanceOnly;
  //  if (fAcceptanceOnly)
  //    G4cout << " ==> tracks killed in the OVC body!" << G4endl;
  //  else
  //    G4cout << G4endl;

 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSStackingAction::~TCSStackingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
TCSStackingAction::ClassifyNewTrack(const G4Track* track)
{
  if (fPrimaryOnly) {
    if (track->GetParentID() != 0) return fKill;   //kill secondaries
  }

  //  if (fAcceptanceOnly) {
  //if (step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==
  //	"OVCCan_PV") return fKill;
  //  }

  //kill secondary neutrino
  //  if (track->GetDefinition() == G4NeutrinoE::NeutrinoE()) return fKill;
  //  else return fUrgent;

  return fUrgent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
