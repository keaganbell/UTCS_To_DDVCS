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

#ifndef TCSCalorimeterHit_h
#define TCSCalorimeterHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"

class TCSCalorimeterHit : public G4VHit {

public:

  TCSCalorimeterHit();
  TCSCalorimeterHit(const TCSCalorimeterHit&);
  TCSCalorimeterHit(G4int col, G4int row, G4int pid, G4double energy,
		    G4ThreeVector pos, G4int boundary_flag);
  virtual ~TCSCalorimeterHit();
  
  // operators
  const TCSCalorimeterHit& operator=(const TCSCalorimeterHit&);
  G4int operator==(const TCSCalorimeterHit&) const;
  
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  
  // methods from base class
  virtual void Draw();
  virtual void Print();
  
  // Get methods
  G4int GetCol() const { return fCol; };
  G4int GetRow() const { return fRow; };
  G4int GetPID() const { return fPID; };
  G4double GetEnergy() const { return fEnergy; };
  G4ThreeVector GetPos() const { return fPos; };
  G4int GetBoundaryFlag() const { return fBoundaryFlag; };
  
private:

  // fEnergy is kin.energy for particles entering calorimeter,
  // otherwise deposited energy.
  // fBoundaryFlag is true for particles entering calorimeter,
  // otherwise false.
  
  G4int         fCol;
  G4int         fRow;
  G4int         fPID;
  G4double      fEnergy;
  G4ThreeVector fPos;
  G4int        fBoundaryFlag;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<TCSCalorimeterHit> TCSCalorimeterHitsCollection;

extern G4ThreadLocal G4Allocator<TCSCalorimeterHit>* TCSCalorimeterHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* TCSCalorimeterHit::operator new(size_t)
{
  if(!TCSCalorimeterHitAllocator)
    TCSCalorimeterHitAllocator = new G4Allocator<TCSCalorimeterHit>;
  return (void *) TCSCalorimeterHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void TCSCalorimeterHit::operator delete(void *hit)
{
  TCSCalorimeterHitAllocator->FreeSingle((TCSCalorimeterHit*) hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
