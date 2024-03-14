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

#ifndef TCSCalorimeterPMTHit_h
#define TCSCalorimeterPMTHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"

class TCSCalorimeterPMTHit : public G4VHit {

public:

  TCSCalorimeterPMTHit();
  TCSCalorimeterPMTHit(const TCSCalorimeterPMTHit&);
  TCSCalorimeterPMTHit(G4int col, G4int row, G4int pid, int npe,
		       G4ThreeVector pos);
  virtual ~TCSCalorimeterPMTHit();
  
  // operators
  const TCSCalorimeterPMTHit& operator=(const TCSCalorimeterPMTHit&);
  G4int operator==(const TCSCalorimeterPMTHit&) const;
  
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  
  // methods from base class
  virtual void Draw();
  virtual void Print();
  
  // Get methods
  G4int GetCol() const { return fCol; };
  G4int GetRow() const { return fRow; };
  G4int GetPID() const { return fPID; };
  G4double GetNpe() const { return fNpe; };
  G4ThreeVector GetPos() const { return fPos; };
  
private:

  G4int         fCol;
  G4int         fRow;
  G4int         fPID;
  G4double      fNpe;
  G4ThreeVector fPos;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<TCSCalorimeterPMTHit> TCSCalorimeterPMTHitsCollection;

extern G4ThreadLocal G4Allocator<TCSCalorimeterPMTHit>* 
TCSCalorimeterPMTHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* TCSCalorimeterPMTHit::operator new(size_t)
{
  if(!TCSCalorimeterPMTHitAllocator)
    TCSCalorimeterPMTHitAllocator = new G4Allocator<TCSCalorimeterPMTHit>;
  return (void *) TCSCalorimeterPMTHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void TCSCalorimeterPMTHit::operator delete(void *hit)
{
  TCSCalorimeterPMTHitAllocator->FreeSingle((TCSCalorimeterPMTHit*) hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
