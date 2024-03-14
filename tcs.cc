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
//
// Neutral Particle Spectrometer (based on Gabi's Amagnet.cc).
// Roadmap:
// define geometry using gdml
// define magnetic field
// define interpolated mag. field
// define scoring meshes
// define ROOT histograms and trees

#include "TCSDetectorConstruction.hh"
//#include "TCSActionInitialization.hh"
#include "TCSPrimaryGeneratorAction.hh"
#include "TCSRunAction.hh"
#include "TCSEventAction.hh"
#include "TCSTrackingAction.hh"
#include "TCSStackingAction.hh"
#include "TCSSteppingAction.hh"
#include "TCSHistoManager.hh"

// attempt to use the GDML parser
#include <vector>
////#include "G4GDMLParser.hh"

////#ifdef G4MULTITHREADED
////#include "G4MTRunManager.hh"
////#else
#include "G4RunManager.hh"
////#endif

#include "G4LogicalVolumeStore.hh"
#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
////#include "QGSP_BERT.hh"
#include "NPSAddOptics.hh"

#include "G4VModularPhysicsList.hh"
////#include "G4EmStandardPhysics.hh"
//#include "G4EmStandardPhysics_option4.hh"
//#include "G4EmLivermorePhysics.hh"

#include "G4StepLimiterPhysics.hh"
#include "G4ScoringManager.hh"
#include "G4UserLimits.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

int main(int argc,char** argv)
{

  cout << "==> TCS setup: argc = " << argc << endl;
  for (int i=0; i<argc; i++) {
    cout << "==>  argv " << i << " = " << argv[i] << endl;
  }

  // Detect interactive mode (if no arguments) and define UI session

  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
    cout << "==>  UI session defined." << endl;
  }

  // Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  //Initialize Random engine.
  //  G4int myseed = 1234;
  G4int myseed = 0;
  ifstream ifs("random_seeds.dat");
  ifs >> myseed;
  ifs.close();
  G4Random::setTheSeed(myseed);
  G4cout << "===> Random Engine initialized with seed (index) " << myseed
	 << " <===" << G4endl;

  // Construct the default run manager

  G4RunManager* runManager = new G4RunManager;

  TCSDetectorConstruction* setup = new TCSDetectorConstruction;

  runManager->SetUserInitialization(setup);

  //  G4Colour  white   ()              ;  // white
  G4Colour  white   (1.0, 1.0, 1.0) ;  // white
  G4Colour  grey    (0.5, 0.5, 0.5) ;  // grey
  G4Colour  black   (0.0, 0.0, 0.0) ;  // black
  G4Colour  red     (1.0, 0.0, 0.0) ;  // red
  G4Colour  green   (0.0, 1.0, 0.0) ;  // green
  G4Colour  blue    (0.0, 0.0, 1.0) ;  // blue
  G4Colour  cyan    (0.0, 1.0, 1.0) ;  // cyan
  G4Colour  magenta (1.0, 0.0, 1.0) ;  // magenta 
  G4Colour  yellow  (1.0, 1.0, 0.0) ;  // yellow
  G4Colour  copper  (0.5, 0.5, 0.0) ;
  G4Colour  orange  (1.0, 0.5, 0.0) ;

  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  ////  G4VModularPhysicsList* physicsList = new QGSP_BERT;

  physicsList->RegisterPhysics(new NPSAddOptics("Cerenkov & Scintillation") );
  physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  runManager->SetUserInitialization(physicsList);

  // User action initialization
  //  runManager->SetUserInitialization(new TCSActionInitialization());

  // Set a HistoManager.
  // Look for root filename in the input line, otherwise do the default.

  char const *kinFileName="ddvcs_gen.kin_data";
  char const *rootFileName="ddvcs_setup.root";
  if ( argc >= 3 )    kinFileName=argv[2];
  if ( argc >= 4 )    rootFileName=argv[3];

  TCSHistoManager*  histo = new TCSHistoManager(kinFileName,rootFileName);

  // Set user action classes.

  TCSPrimaryGeneratorAction* gen_action = new TCSPrimaryGeneratorAction();
  runManager->SetUserAction(gen_action);

  TCSRunAction* run_action = new TCSRunAction(histo);  
  runManager->SetUserAction(run_action);

  TCSEventAction* event_action = new TCSEventAction(histo);
  runManager->SetUserAction(event_action);

  TCSTrackingAction* tracking_action = new TCSTrackingAction();
  runManager->SetUserAction(tracking_action);

  TCSStackingAction* stacking_action = new TCSStackingAction();
  runManager->SetUserAction(stacking_action);

  TCSSteppingAction* stepping_action = new TCSSteppingAction(event_action);
  runManager->SetUserAction(stepping_action);

  // Visualization manager construction
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
  cout << "==>  VisManager defined & initialized." << endl;
    
  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if ( argc != 1 ) {

    // execute an argument macro file if exist
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);

  }
  else {

    UImanager->ApplyCommand("/control/execute ./init_vis.mac"); 
    cout << "==>  init_vis.mac executed." << endl;

    // start interactive session
    ui->SessionStart();
    cout << "==>  UI session started." << endl;
    delete ui;

  }   //argc!=1

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
