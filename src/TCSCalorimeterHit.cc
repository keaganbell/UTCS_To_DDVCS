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

#include "TCSCalorimeterHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

G4ThreadLocal G4Allocator<TCSCalorimeterHit>* TCSCalorimeterHitAllocator=0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSCalorimeterHit::TCSCalorimeterHit() : G4VHit(),
	 fCol(-1), fRow(-1), fPID(0), fEnergy(999999.), fPos(G4ThreeVector()),
	 fBoundaryFlag(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSCalorimeterHit::~TCSCalorimeterHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSCalorimeterHit::TCSCalorimeterHit(const TCSCalorimeterHit& right) : G4VHit()
{
  fCol    = right.fCol;
  fRow    = right.fRow;
  fPID    = right.fPID;
  fEnergy   = right.fEnergy;
  fPos    = right.fPos;
  fBoundaryFlag = right.fBoundaryFlag;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSCalorimeterHit::TCSCalorimeterHit(G4int col, G4int row, G4int pid,
	     G4double energy, G4ThreeVector pos, G4int boundary_flag) {
  fCol    = col;
  fRow    = row;
  fPID    = pid;
  fEnergy   = energy;
  fPos    = pos;
  fBoundaryFlag = boundary_flag;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const TCSCalorimeterHit&
TCSCalorimeterHit::operator=(const TCSCalorimeterHit& right)
{
  fCol    = right.fCol;
  fRow    = right.fRow;
  fPID    = right.fPID;
  fEnergy = right.fEnergy;
  fPos    = right.fPos;
  fBoundaryFlag = right.fBoundaryFlag;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int TCSCalorimeterHit::operator==(const TCSCalorimeterHit& right) const
{
  return ( this == &right ) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSCalorimeterHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(fPos);
    circle.SetScreenSize(4.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSCalorimeterHit::Print()
{
  G4cout << "TCSCalorimeterHit: col = " << fCol << "  row = " << fRow
	 << "  particle id = " << fPID << "  energy = " << fEnergy
	 << "  boundary flag = " << fBoundaryFlag << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
