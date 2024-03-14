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

#include "TCSCalorimeterPMTHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

G4ThreadLocal G4Allocator<TCSCalorimeterPMTHit>*
TCSCalorimeterPMTHitAllocator=0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSCalorimeterPMTHit::TCSCalorimeterPMTHit() : G4VHit(),
	       fCol(-1), fRow(-1), fPID(0), fNpe(-1), fPos(G4ThreeVector())
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSCalorimeterPMTHit::~TCSCalorimeterPMTHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSCalorimeterPMTHit::TCSCalorimeterPMTHit(const TCSCalorimeterPMTHit& right) :
  G4VHit()
{
  fCol    = right.fCol;
  fRow    = right.fRow;
  fPID    = right.fPID;
  fNpe    = right.fNpe;
  fPos    = right.fPos;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSCalorimeterPMTHit::TCSCalorimeterPMTHit(G4int col, G4int row, G4int pid,
					   int Npe, G4ThreeVector pos) {
  fCol    = col;
  fRow    = row;
  fPID    = pid;
  fNpe    = Npe;
  fPos    = pos;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const TCSCalorimeterPMTHit&
TCSCalorimeterPMTHit::operator=(const TCSCalorimeterPMTHit& right)
{
  fCol    = right.fCol;
  fRow    = right.fRow;
  fPID    = right.fPID;
  fNpe    = right.fNpe;
  fPos    = right.fPos;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int TCSCalorimeterPMTHit::operator==(const TCSCalorimeterPMTHit& right) const
{
  return ( this == &right ) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSCalorimeterPMTHit::Draw()
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

void TCSCalorimeterPMTHit::Print()
{
  G4cout << "TCSCalorimeterPMTHit: col = " << fCol << "  row = " << fRow
	 << "  particle id = " << fPID << "  Npe = " << fNpe;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
