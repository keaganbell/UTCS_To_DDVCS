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

#ifndef TCSEventAction_h
#define TCSEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "TCSHistoManager.hh"
#include "TCSHodoscopeHit.hh"
#include "TCSTrackerHit.hh"

// Event action class

class G4PrimaryParticle;
class TCSHistoManager;

class TCSEventAction : public G4UserEventAction {

public:
  TCSEventAction(TCSHistoManager*);
  virtual ~TCSEventAction();
    
  virtual void BeginOfEventAction(const G4Event* event);
  virtual void EndOfEventAction(const G4Event* event);
  TCSHistoManager *GetHistoManager() {return fHistoManager;}
  void AddEdep(G4double edep) { fEdep += edep; }
  void AddNTar_epos(int n) { fNTar_epos += n; }
  void AddNTar_eneg(int n) { fNTar_eneg += n; }
  G4int GetEvtNo() {return fEvtNo;}
  
private:

  G4double  fEdep;
  int fNTar_epos;    //number of positrons from gamma conversion in target.
  int fNTar_eneg;    //number of electrons from gamma conversion in target.
  TCSHistoManager *fHistoManager;
  G4int fPrintModulo;
  G4int fTargetCollID;
  G4int fCalorimeterCollID;
  G4int fCalorimeterPMTCollID;
  G4int fHodoscopeCollID;
  G4int fTrackerCollID;
  G4int fEvtNo;

  //  void AddHodoHit(TCSHodoHitsCollection* HC, HodoHitContainer& HodoHitCont);
  void AddTrackerHit(TCSTrackerHitsCollection* HC,
  		     TrackerHitContainer& TrackerHitCont);

  static int GetQuarter(double x, double y);

  friend class TCSTrackerSD;   //to use GetQuarter
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
