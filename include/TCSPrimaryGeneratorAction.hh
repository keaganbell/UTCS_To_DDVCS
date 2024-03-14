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

#ifndef TCSPrimaryGeneratorAction_h
#define TCSPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include <fstream>

//#include "TCSGen.hh"

using namespace std;

class G4ParticleGun;
class G4VPrimaryGenerator;
class G4Event;
//class G4Box;

/// The primary generator action class with particle gun.

class TCSPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {

public:

  TCSPrimaryGeneratorAction();    
  virtual ~TCSPrimaryGeneratorAction();

  // method from the base class
  virtual void GeneratePrimaries(G4Event*);         
  
  // method to access particle gun
  const G4ParticleGun* GetParticleGun() const { return fParticleGun; }

  const G4VPrimaryGenerator* GetPrimaryGenerator() const { return fHEPEvt; }

  enum mode {beam, brem, tcs, hep_data};
  mode GetMode() {return fMode;}

  double GetBremEnergy(double Ee, double Eg_min, double Eg_max);
  G4ThreeVector GenRandomDirection(const G4ThreeVector Axis,
				   const double Theta, const double Phi);

private:

  G4ParticleGun*  fParticleGun;    // pointer a to G4 gun class

  // Beam default parameter.
  string fParticleName;
  G4double fPmin, fPmax;           // momentum range
  G4double fX0, fY0, fZ0;          // beam position
  G4double fDX, fDY, fDZ;          // beam dimensions
  G4double fThetaX, fThetaY;       // beam deflection angles in X and Y [deg]
  G4double fdThetaX;               // beam divegence (+/-dThetaX if dThetaY>=0)
  G4double fdThetaY;               // beam divegence (+/-dThetaY if >=0)
  mode fMode;                      // beam/tcs mode  

  //  TCSGen fTCSGen;
  //  long fTCSEntryNum;
  //  uint fTCSPartNum;
  G4VPrimaryGenerator* fHEPEvt;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
