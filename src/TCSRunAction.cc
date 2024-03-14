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

#include "TCSRunAction.hh"
#include "TCSPrimaryGeneratorAction.hh"
#include "TCSDetectorConstruction.hh"
#include "TCSRun.hh"
#include "TCSHistoManager.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <ctime>
#include "G4Timer.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSRunAction::TCSRunAction(TCSHistoManager* histo)
  : G4UserRunAction(), fHistoManager(histo), fAutoSeed(false)
//  : G4UserRunAction(), fHistoManager(histo), fAutoSeed(true)

{ 
  // add new units for dose
  // 
  const G4double kilogray  = 1.e+3*gray;
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;

  new G4UnitDefinition("kilogray",  "kiloGy"  , "Dose", kilogray);
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray);        

  fTimer = new G4Timer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSRunAction::~TCSRunAction()
{
  delete fTimer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* TCSRunAction::GenerateRun()
{
  return new TCSRun; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSRunAction::BeginOfRunAction(const G4Run*)
{ 

  long seeds[2] = {0, 0};;
  if (fAutoSeed) {
    // automatic (time-based) random seeds for each run
    G4cout << "*************************************" << G4endl;
    G4cout << "*** AUTOSEED ON: Time-based seeds ***" << G4endl;
    G4cout << "*************************************" << G4endl;
    time_t systime = time(NULL);
    seeds[0] = (long) systime;
    seeds[1] = (long) (systime*G4UniformRand());
  } else {
    G4cout << "******************************************" << G4endl;
    G4cout << "*** AUTOSEED OFF: Seeds from data file ***" << G4endl;
    G4cout << "******************************************" << G4endl;
    ifstream ifs("random_seeds.dat");
    ifs >> seeds[0] >> seeds[1];
    ifs.close();
  }

  G4Random::setTheSeeds(seeds);
  G4Random::showEngineStatus();

  //histograms
  fHistoManager->book(); 

  fTimer->Start();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSRunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  const TCSRun* b1Run = static_cast<const TCSRun*>(run);

  // Compute dose
  /*
  G4double edep  = b1Run->GetEdep();
  G4double edep2 = b1Run->GetEdep2();
  G4double rms = edep2 - edep*edep/nofEvents;
  if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;

  const TCSDetectorConstruction* detectorConstruction
   = static_cast<const TCSDetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double mass = detectorConstruction->GetScoringVolume()->GetMass();
  G4double dose = edep/mass;
  G4double rmsDose = rms/mass;
  */
  
  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.

  /*
  const TCSPrimaryGeneratorAction* generatorAction
   = static_cast<const TCSPrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());

  G4String runCondition;

  if (generatorAction)
  {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }
  */        

  // Print

  cout << "# of e- from gamma conversion in target: " << b1Run->GetNTar_eneg()
       << endl;
  cout << "# of e+ from gamma conversion in target: " << b1Run->GetNTar_epos()
       << endl;
    
  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------";
  }

  /*  
  G4cout
     << G4endl
     << " The run consists of " << nofEvents << " "<< runCondition
     << G4endl
     << " Dose in scoring volumes : " 
     << G4BestUnit(dose,"Dose") << " +- " << G4BestUnit(rmsDose,"Dose")
     << G4endl
     << "    Deposited energy = " << G4BestUnit(edep,"Energy")
     << G4endl
     << "    Mass             = " << G4BestUnit(mass,"Mass")
     << G4endl
     << "------------------------------------------------------------"
     << G4endl
     << G4endl;

  dose *= 1./nofEvents * 1./eplus * microampere * 3600. * second;
  rmsDose *= 1./nofEvents * 1./eplus * microampere * 3600. * second;
  G4cout << " Dose in scoring volume : " 
	 << G4BestUnit(dose,"Dose") << " +- " << G4BestUnit(rmsDose,"Dose")
	 << "/(uA*hour)" << G4endl;
  */

  //save histograms
  //
  fHistoManager->PrintStatistic();
  fHistoManager->save();   

  fTimer->Stop();
  G4cout << "Number of events " << nofEvents << ",  " << *fTimer << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
