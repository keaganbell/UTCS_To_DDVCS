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

#include "TCSDetectorConstruction.hh"
#include "TCSCalorimeterConstruction.hh"
#include "TCSHodoscopeConstruction.hh"
#include "TCSTargetSD.hh"
#include "TCSCalorimeterSD.hh"
#include "TCSCalorimeterPMTSD.hh"
#include "TCSHodoscopeSD.hh"
#include "TCSTrackerSD.hh"
#include "G4SDManager.hh"

#include "TCSMuonDetectorConstruction.hh"

// **** Magnetic field ******
// New include files - used for magnetic field
#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "SimpleField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"
#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4NystromRK4.hh"
#include "G4SimpleHeum.hh"
#include "G4SimpleRunge.hh"
//#include "G4ChordFinder.hh"
//#include "G4EqMagElectricField.hh"
#include "G4PropagatorInField.hh"
#include "G4UserLimits.hh"

// Shapes...
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4UnionSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4TwoVector.hh"
#include "G4ExtrudedSolid.hh"
//#include "UExtrudedSolid.hh"
#include "G4SystemOfUnits.hh"

// Others...
#include "G4AutoDelete.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSDetectorConstruction::TCSDetectorConstruction()
  : G4VUserDetectorConstruction(), fScoringVolume(0), fField(0), fEquation(0),
    fStepper(0), fChordFinder(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSDetectorConstruction::~TCSDetectorConstruction()
{
  delete fField;
  delete fEquation;
  delete fStepper;
  delete fChordFinder;
  delete fScoringVolume;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* TCSDetectorConstruction::Construct()
{

  G4NistManager* man = G4NistManager::Instance();
  //  man->SetVerbose(1);
  G4Material* Air = man->FindOrBuildMaterial("G4_AIR");

  // Create the hall
  //
  //fXWorld, fYWorld, fZWorld are defined in include/TCSDetectorConstruction.hh
  //worldBox defines outer most cubic boundary : Unique volume '; contains all other volumes; define the global coordinate system with origin at the center of the world volume; position of tracks defines wrt global co ordinate
  G4Box * WorldBox = new G4Box("WorldBox",fXWorld/2, fYWorld/2, fZWorld/2);
  G4LogicalVolume * WorldLV = new G4LogicalVolume(WorldBox, Air,"WorldLV");
  physWorld = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), WorldLV, "World",
				0, false, 0);

  // Construct calorimeter.
  
  TCSCalorimeterConstruction CalorimeterConstruction; //temp off may23
  CalorimeterConstruction.Construct(); //temp off may23
  
  //GetCalorimeter is defined in CalorimeterConstruction.cc // temp off may23
  G4LogicalVolume* Calorimeter_log = CalorimeterConstruction.GetCalorimeter();//temp off may23
  
  //Calo is a structure defined in TCSDetectorConstruction.hh
  //NQuarter is defined inside the Calo strcuture and is const int with value 4
  
  for (int quarter=0; quarter<Calo.NQuarter; quarter++) //temp off may23, 2023
  PositionCalorimeter(Calorimeter_log, quarter); //temp off may23
  
  // Construct hodoscope.

  //TCSHodoscopeConstruction HodoscopeConstruction; //temp off
  //HodoscopeConstruction.Construct(); //temp off
  //G4LogicalVolume* Hodoscope_log = HodoscopeConstruction.GetHodoscope(); //temp off

  //for (int quarter=0; quarter<Calo.NQuarter; quarter++) //temp off
  //PositionHodoscope(Hodoscope_log, quarter); //temp off

  
  // Read trackers from the gdml files.

  /*    //temp off may23, 2023
  G4LogicalVolume* Tracker1_log = GetGDMLVolume(
		   "tcs_gdmls/pointer_referenced/tracker1_ref.gdml",
		   "TrackerAssembly0xe91030");
  //		   "tracker1World0xe915c0");
  for (int quarter=0; quarter<Calo.NQuarter; quarter++)
    PositionTracker(Tracker1_log, quarter, 1);

  G4LogicalVolume* Tracker2_log = GetGDMLVolume(
		   "tcs_gdmls/pointer_referenced/tracker2_ref.gdml",
		   "TrackerAssembly0x1a10030");
		   //		   "tracker2World0x1a105c0");
  for (int quarter=0; quarter<Calo.NQuarter; quarter++)
    PositionTracker(Tracker2_log, quarter, 2);

  G4LogicalVolume* Tracker3_log = GetGDMLVolume(
		   "tcs_gdmls/pointer_referenced/tracker3_ref.gdml",
		   "TrackerAssembly0x1760030");
		   //		   "tracker3World0x17605c0");
  for (int quarter=0; quarter<Calo.NQuarter; quarter++)
    PositionTracker(Tracker3_log, quarter, 3);
  */ //temp off may23, 2023

  ///  G4LogicalVolume* Tracker4_log = GetGDMLVolume(
  ///		   "tcs_gdmls/pointer_referenced/tracker4_ref.gdml",
  ///		   "TrackerAssembly0x1951ec0");
		   //		   "tracker3World0x17605c0");
  ///  for (int quarter=0; quarter<Calo.NQuarter; quarter++)
  ///    PositionTracker(Tracker4_log, quarter, 4);

  
  //Scattering chamber.
  
  fParser.ReadModule("tcs_gdmls/scattering_chamber_hallac_may4_2023.gdml"); //deb see next line
  //fParser.ReadModule("tcs_gdmls/scattering_chamber_hallac_may4_2023_test.gdml");
  //fParser.ReadModule("tcs_gdmls/target_hallac.gdml");
  G4LogicalVolume* chamber_log = fParser.GetVolume("ChamberAssembly"); //deb see next line 
  //G4LogicalVolume* chamber_log = fParser.GetVolume("TargetAssembly");
  G4RotationMatrix*  cham_rot = new G4RotationMatrix();
  //cham_rot->rotateY(90.*degree); //to not rotate the target
  cham_rot->rotateY(0.*degree);
  new G4PVPlacement(cham_rot,
		    G4ThreeVector(0.,0.,0.),
		    chamber_log,
		    "ChamberAssembly_PV",
		    //"Assembly_PV",
		    physWorld->GetLogicalVolume(), //its mother  volume
		    false,
		    0,
                    false);
  // Setup Magnetic Field here!!!
  ConstructField(); // deb commented 

  // Setup sensitive detector!
  ConstructSDandField(); //deb commented 
  
  //--implement the beam pipe--//

  // beam pipe commentout : oct 31
  
  
  fParser.ReadModule("tcs_gdmls/beam_pipe.gdml");
  G4LogicalVolume* BeamPipe_log = fParser.GetVolume("BeamPipeAssembly");
  G4RotationMatrix*  BeamPipe_rot = new G4RotationMatrix();
  BeamPipe_rot->rotateY(0.*degree); 

  double beamp_shift = 1.01*(beam_pipe_pos.rin_cham_len + 
                       beam_pipe_pos.thick_cham_len + 
                      (beam_pipe_pos.len1 +
                       beam_pipe_pos.len2 +
		       beam_pipe_pos.len3 +
		       beam_pipe_pos.len4)/2.);

  new G4PVPlacement(BeamPipe_rot,
		    G4ThreeVector(0.,0.,beamp_shift),
		    BeamPipe_log,
		    "BeamPipeAssembly_PV",
		    physWorld->GetLogicalVolume(), //its mother  volume
		    false,
		    0,
                    false);
  
  
  //---simple magnetic field mimicking SBS magnets---//
  //---need to change later---//

  
  fParser.ReadModule("tcs_gdmls/simple_magnet.gdml");
  G4LogicalVolume* SBSMagnet_log = fParser.GetVolume("MagnetAssembly");
  G4RotationMatrix*  SBSMagnet_rot = new G4RotationMatrix();
  SBSMagnet_rot->rotateY(0.); 
  new G4PVPlacement(SBSMagnet_rot,
		    G4ThreeVector(0.,0.,SBS_mag.Distance),
		    SBSMagnet_log,
		    "MagnetAssembly_PV",
		    physWorld->GetLogicalVolume(), //its mother  volume
		    false,
		    0,
                    false);
  
  
  
  //---create simple magnetic field--//
  G4Box * FieldVolumeBox = new G4Box("FieldVolumeBox", (SBS_mag.mag_vol_x)/2., (SBS_mag.mag_vol_y)/2., (SBS_mag.mag_vol_z)/2.);
   
  G4Tubs* cyl = new G4Tubs("Cylinder",0, 8.2*cm, 1.01*(SBS_mag.mag_vol_z)/2, 0., 360.*deg); 
  //format is ("Cylinder", r_min, r_max, z:+z to -z, phi_min, phi_max)

  G4SubtractionSolid* FieldVolumeBox_sub =
  new G4SubtractionSolid("FieldVolumeBox-Cyl", FieldVolumeBox, cyl);

  cout << "SBS_mag.mag_vol_x = " << SBS_mag.mag_vol_x << endl;

  //G4LogicalVolume * FieldVolume_LV = new G4LogicalVolume(FieldVolumeBox_sub, Air,"FieldVolume_LV");
  G4LogicalVolume * FieldVolume_LV = new G4LogicalVolume(FieldVolumeBox, Air,"FieldVolume_LV"); //oct 31 second blob
  G4RotationMatrix*  FieldVolume_rot = new G4RotationMatrix();
  FieldVolume_rot->rotateY(0.); 
  new G4PVPlacement(FieldVolume_rot,
		    G4ThreeVector(0.,0.,SBS_mag.Distance),
		    FieldVolume_LV,
		    "FieldVolume_PV",
		    physWorld->GetLogicalVolume(), //its mother  volume
		    false,
		    0,
                    false);
  
  // Define Magnetic Field//
  //magnet bore z direction length = 1.2192 m 

  //G4ThreeVector fieldVector((2.4/1.2192)*tesla, 0.0, 0.0);
  //G4ThreeVector fieldVector(2.4*tesla, 0.0, 0.0);
  G4ThreeVector fieldVector(0.0*tesla, 0.0, 0.0);
  G4UniformMagField* SBSmagneticField = new G4UniformMagField(fieldVector);

  G4FieldManager* fieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fieldManager->SetDetectorField(SBSmagneticField);
  fieldManager->CreateChordFinder(SBSmagneticField);

  FieldVolume_LV->SetFieldManager(fieldManager, true);

  
  // G4Box* SBS_mag_field_volume = new G4Box("SBS_mag_field_volume", 469.9*mm, 1.0*cm, 1.0*cm);

  /// Limit tracking steps in rapid canging magnetic field close to target.
  //fParser.GetVolume("TargetMat")->SetUserLimits(new G4UserLimits(1.*mm));
  //fParser.GetVolume("TargetAssembly")->SetUserLimits(new G4UserLimits(1.*mm));
  //fParser.GetVolume("MagPoleAssembly")->SetUserLimits(new G4UserLimits(1.*mm));
  //fParser.GetVolume("LN2ShieldAssembly")->SetUserLimits(new G4UserLimits(1.*mm));
  //chamber_log->SetUserLimits(new G4UserLimits(1.*mm));

  //always return the physical World
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSDetectorConstruction::ConstructField() 
{
  fField = new SimpleField();

  static G4TransportationManager* trMgr= 
    G4TransportationManager::GetTransportationManager(); 

  // A field object is held by a field manager
  // Find the global Field Manager
  G4FieldManager* globalFieldMgr= trMgr->GetFieldManager();

  // Relative accuracy values:
  G4double minEps= 1.0e-5;  //   Minimum & value for largest steps
  G4double maxEps= 1.0e-4;  //   Maximum & value for smallest steps

  globalFieldMgr->SetMinimumEpsilonStep( minEps );
  globalFieldMgr->SetMaximumEpsilonStep( maxEps );
  globalFieldMgr->SetDeltaOneStep( 0.5e-3 * mm );  // 0.5 micrometer

  globalFieldMgr->SetDetectorField(fField);

  globalFieldMgr->CreateChordFinder(fField);

  globalFieldMgr->GetChordFinder()->SetDeltaChord(1e-4*m);

  fEquation = new G4Mag_UsualEqRhs (fField);

  //  fStepper = new G4ClassicalRK4 (fEquation);
  //  fStepper = new G4NystromRK4 (fEquation);
  //  fStepper = new G4HelixImplicitEuler(fEquation);
  fStepper = new G4HelixExplicitEuler(fEquation);

  
  /*
  static G4TransportationManager* trMgr= 
    G4TransportationManager::GetTransportationManager(); 



    //The ChordFinder is an helper class to track particles 
    //in magnetic fields, it sets the accuracy to be used.

    fEquation = new G4Mag_UsualEqRhs (fField);

    fStepper = new G4ClassicalRK4 (fEquation);    //Default choice.
    //    fStepper = new G4ImplicitEuler (fEquation);    //300 ev/min
    //    fStepper = new G4ExplicitEuler (fEquation);    //<300 ev/min
    //    fStepper = new G4SimpleRunge (fEquation);    //300 ev/min
    //    fStepper = new G4DormandPrinceRK56 (fEquation);  //can't link

    //Mag. field
    //    fStepper = new G4HelixImplicitEuler (fEquation);    //does not work
    //    fStepper = new G4HelixExplicitEuler (fEquation);    //slow
    //    fStepper = new G4HelixSimpleRunge (fEquation);    //does not work
    //    fStepper = new G4NystromRK4 (fEquation);    //slow
    //    fStepper = new G4SimpleHeum (fEquation);    //300 ev/min

    ///    fChordFinder = new G4ChordFinder(fField,1e-5*m,fStepper);
    fChordFinder = new G4ChordFinder(fField,1e-4*m,fStepper);
    globalFieldMgr->SetChordFinder(fChordFinder);
    globalFieldMgr->GetChordFinder()->SetDeltaChord(1e-4*m);
    globalFieldMgr->SetDeltaIntersection(1e-4*m);
    globalFieldMgr->SetDeltaOneStep(1e-4*m);
    //    globalFieldMgr->SetEpsilonMax(1e-4);  not present in 10.4
    //    globalFieldMgr->SetEpsilonMin(1e-5);

    G4cout << "Magnetic field has been constructed " << 
      "in TCSDetectorConstruction::ConstructField()" << G4endl;
    fieldIsInitialized = true; 
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void TCSDetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

  //Avoid double initialization.
  static G4ThreadLocal G4bool initialized = false;
  if ( ! initialized ) {

    // Calorimeter SD

    TCSCalorimeterSD* caloSD = new TCSCalorimeterSD("CalorimeterSD",
						   "CalorimeterHitsCollection");
    G4SDManager::GetSDMpointer()->AddNewDetector(caloSD);

    // Register the messenger for deleting
    //  G4AutoDelete::Register(caloSD);
  
    //  SetSensitiveDetector("CalorimeterAssembly", tcsSD, true);
    ////    SetSensitiveDetector("Block", caloSD, true);
    ////    SetSensitiveDetector("caloWorld", caloSD, true);

    SetSensitiveDetector("Block_log", caloSD, true); //temp off may23, 2023
    SetSensitiveDetector("Counter_log", caloSD, true); //temp off may23, 2023
    SetSensitiveDetector("Calorimeter_LV", caloSD, true); //temp off may23, 2023

    TCSCalorimeterPMTSD* caloPMTSD = new TCSCalorimeterPMTSD("CalorimeterPMTSD",
    				     "CalorimeterPMTHitsCollection"); //temp off may23, 2023
    G4SDManager::GetSDMpointer()->AddNewDetector(caloPMTSD); //temp off may23, 2023

    SetSensitiveDetector("Cathode_log", caloPMTSD, true); //temp off may23, 2023

    // Hodoscope SD

    //TCSHodoscopeSD* hodoSD = new TCSHodoscopeSD("HodoscopeSD", //temp off
    //"HodoscopeHitsCollection"); //temp off
    //G4SDManager::GetSDMpointer()->AddNewDetector(hodoSD); //temp off 
    //SetSensitiveDetector("hodo_module_LV", hodoSD, true); //temp off
    //SetSensitiveDetector("Hodoscope_LV", hodoSD, true); //temp off

    // Tracker SD

    // TCSTrackerSD* trackerSD = new TCSTrackerSD("TrackerSD", //temp off may23, 2023
    //					       "TrackerHitsCollection"); //temp off may23, 2023
    // G4SDManager::GetSDMpointer()->AddNewDetector(trackerSD); //temp off may23, 2023
    //    SetSensitiveDetector("Drift", trackerSD, true);
    //SetSensitiveDetector("Drift0xe15800", trackerSD, true);   //tracker 1 //temp off may23, 2023
    //SetSensitiveDetector("Drift0x1994830", trackerSD, true);  //tracker 2 //temp off may23, 2023
    //SetSensitiveDetector("Drift0x16e4830", trackerSD, true);  //tracker 3 //temp off may23, 2023
    ///    SetSensitiveDetector("Drift0x19512f0", trackerSD, true);  //tracker 4
    //    SetSensitiveDetector("tracker1World0xe915c0", trackerSD, true);
    //    SetSensitiveDetector("tracker2World0x1a105c0", trackerSD, true);
    //    SetSensitiveDetector("tracker3World0x17605c0", trackerSD, true);

     TCSTargetSD* targetSD = new TCSTargetSD("TargetSD",
					    "TargetHitsCollection");
     G4SDManager::GetSDMpointer()->AddNewDetector(targetSD);

     SetSensitiveDetector("TargetAssembly", targetSD, true);

    // Create global magnetic field messenger.
    // Uniform magnetic field is then created automatically if
    // the field value is not zero.
    //  G4ThreeVector fieldValue = G4ThreeVector();
    //fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
    //fMagFieldMessenger->SetVerboseLevel(1);
  
    // Register the field messenger for deleting
    //G4AutoDelete::Register(fMagFieldMessenger);

    initialized=true;
  }

}

//==============================================================================

void TCSDetectorConstruction::PositionCalorimeter(
			    G4LogicalVolume* Calorimeter_log, int quarter) {

  // Positioning of the calorimeter quarters.

  ///  double phi, theta_tilt, theta_pos;
    double phi, theta_tilt;

  switch (quarter) {
  case 0:
    phi        = Calo.RotationAngle;
    theta_tilt = -Calo.TiltAngle;
    ///    theta_pos  =  Calo.PositionAngle;
    break;
  case 1:
    phi        = -Calo.RotationAngle;
    theta_tilt = -Calo.TiltAngle;
    ///    theta_pos  =  Calo.PositionAngle;
    break;
  case 2:
    phi        = -Calo.RotationAngle;
    theta_tilt =  Calo.TiltAngle;
    ///    theta_pos  = -Calo.PositionAngle;
    break;
  case 3:
    phi        =  Calo.RotationAngle;
    theta_tilt =  Calo.TiltAngle;
    ///    theta_pos  = -Calo.PositionAngle;
    break;
  default:
    phi        = 0.;
    theta_tilt = 0.;
    ///    theta_pos  = 0.;
  }


  /*
  //This is how it should be.
 // u, v, w are the daughter axes, projected on the mother frame
  G4ThreeVector u = G4ThreeVector(cos(phi), 0.,-sin(phi));
  G4ThreeVector v = G4ThreeVector(-sin(theta_tilt)*sin(phi), cos(theta_tilt),
				  -sin(theta_tilt)*cos(phi));
  G4ThreeVector w = G4ThreeVector(sin(phi), 0., cos(phi));
  G4RotationMatrix rotm = G4RotationMatrix(u, v, w);
  G4cout << "Direct rotation matrix : ";
  rotm.print(G4cout);     
  */

  //This is consistent with gdml coding.
  G4RotationMatrix rotm;
  rotm.rotateY(phi);
  rotm.rotateX(theta_tilt);
  ///  rotm.rotateZ(0.);

  ///G4ThreeVector position=G4ThreeVector(sin(phi)*cos(theta_pos), sin(theta_pos),
  ///				       cos(phi)*cos(theta_pos));
  ///  position *= Calo.Distance;
  G4ThreeVector position(0., 0., Calo.Distance);
  cout <<"Calorimeter Z Distance = " << Calo.Distance << endl;
  position *= rotm;

  G4Transform3D transform = G4Transform3D(rotm, position);

  new G4PVPlacement(transform,                     //position, rotation        
                    Calorimeter_log,               //logical volume
                    "Calorimeter_PV",              //name
		    physWorld->GetLogicalVolume(), //its mother  volume
                    false,                         //no boolean operation
                    quarter);                      //copy number
}
//==============================================================================

void TCSDetectorConstruction::PositionHodoscope(
			    G4LogicalVolume* Hodoscope_log, int quarter) {

  // Positioning of the hodoscope quarters.

  ///  double phi, theta_tilt, theta_pos;
    double phi, theta_tilt;

  switch (quarter) {
  case 0:
    phi        =  Hodo.RotationAngle;
    theta_tilt = -Hodo.TiltAngle;
    ///    theta_pos  =  Hodo.PositionAngle;
    break;
  case 1:
    phi        = -Hodo.RotationAngle;
    theta_tilt = -Hodo.TiltAngle;
    ///    theta_pos  =  Hodo.PositionAngle;
    break;
  case 2:
    phi        = -Hodo.RotationAngle;
    theta_tilt =  Hodo.TiltAngle;
    ///    theta_pos  = -Hodo.PositionAngle;
    break;
  case 3:
    phi        =  Hodo.RotationAngle;
    theta_tilt =  Hodo.TiltAngle;
    ///    theta_pos  = -Hodo.PositionAngle;
    break;
  default:
    phi        = 0.;
    theta_tilt = 0.;
    ///    theta_pos  = 0.;
  }

  //This is consistent with gdml coding.
  G4RotationMatrix rotm;
  rotm.rotateY(phi);
  rotm.rotateX(theta_tilt);
  ///  rotm.rotateZ(0.);

  ///G4ThreeVector position=G4ThreeVector(sin(phi)*cos(theta_pos), sin(theta_pos),
  ///				       cos(phi)*cos(theta_pos));
  ///  position *= Hodo.Distance;
  G4ThreeVector position(0., 0., Hodo.Distance);
  position *= rotm;

  G4Transform3D transform = G4Transform3D(rotm, position);

  new G4PVPlacement(transform,                     //position, rotation        
                    Hodoscope_log,               //logical volume
                    "Hodoscope_PV",              //name
		    physWorld->GetLogicalVolume(), //its mother  volume
                    false,                         //no boolean operation
                    quarter);                      //copy number
}

//==============================================================================


void TCSDetectorConstruction::PositionMuon(
					   G4LogicalVolume* Muon_log, int quarter) { // Do not include int quarter now

  
  double phi, theta_tilt;
  
  cout <<"Position Muon Detector" << endl; // check 
  phi        =  Muon.RotationAngle; // check
  theta_tilt = -Muon.TiltAngle; // check 
  

//This is consistent with gdml coding.
 G4RotationMatrix rotm;
 rotm.rotateY(phi);
 rotm.rotateX(theta_tilt);
  ///  rotm.rotateZ(0.);

  ///G4ThreeVector position=G4ThreeVector(sin(phi)*cos(theta_pos), sin(theta_pos),
  ///				       cos(phi)*cos(theta_pos));
  ///  position *= Hodo.Distance;
  G4ThreeVector position(0., 0., Muon.Distance);
 //G4ThreeVector position(0., 0., 12.m);

 position *= rotm;

  G4Transform3D transform = G4Transform3D(rotm, position);

  new G4PVPlacement(transform,                     //position, rotation        
                    Muon_log,               //logical volume
                    "Muon_PV",              //name
		    physWorld->GetLogicalVolume(), //its mother  volume
                    false,                         //no boolean operation
                    quarter);                      //copy number
}


//==============================================================================

void TCSDetectorConstruction::PositionTracker(G4LogicalVolume* Tracker_log,
					      int quarter, int layer) {

  // Positioning of the tracker quarters.

  ///  double phi, theta_tilt, theta_pos;
    double phi, theta_tilt;

  switch (quarter) {
  case 0:
    phi        =  Tracker.RotationAngle;
    theta_tilt = -Tracker.TiltAngle;
    ///    theta_pos  =  Tracker.PositionAngle;
    break;
  case 1:
    phi        = -Tracker.RotationAngle;
    theta_tilt = -Tracker.TiltAngle;
    ///    theta_pos  =  Tracker.PositionAngle;
    break;
  case 2:
    phi        = -Tracker.RotationAngle;
    theta_tilt =  Tracker.TiltAngle;
    ///    theta_pos  = -Tracker.PositionAngle;
    break;
  case 3:
    phi        =  Tracker.RotationAngle;
    theta_tilt =  Tracker.TiltAngle;
    ///    theta_pos  = -Tracker.PositionAngle;
    break;
  default:
    phi        = 0.;
    theta_tilt = 0.;
    ///    theta_pos  = 0.;
  }

  /*
  //This is how it should be.
 // u, v, w are the daughter axes, projected on the mother frame
  G4ThreeVector u = G4ThreeVector(cos(phi), 0.,-sin(phi));
  G4ThreeVector v = G4ThreeVector(-sin(theta_tilt)*sin(phi), cos(theta_tilt),
				  -sin(theta_tilt)*cos(phi));
  G4ThreeVector w = G4ThreeVector(sin(phi), 0., cos(phi));
  G4RotationMatrix rotm = G4RotationMatrix(u, v, w);
  G4cout << "Direct rotation matrix : ";
  rotm.print(G4cout);     
  */

  //This is consistent with gdml coding.
  G4RotationMatrix rotm;
  rotm.rotateY(phi);
  rotm.rotateX(theta_tilt);
  ///  rotm.rotateZ(0.);

  //G4ThreeVector position=G4ThreeVector(sin(phi)*cos(theta_pos), sin(theta_pos),
  //				       cos(phi)*cos(theta_pos));
  //  position *= Tracker.Distance[layer-1];
  G4ThreeVector position(0., 0., Tracker.Distance[layer-1]);
  position *= rotm;

  G4Transform3D transform = G4Transform3D(rotm, position);

  new G4PVPlacement(transform,                     //position, rotation        
                    Tracker_log,                   //logical volume
                    "Tracker"+to_string(layer),    //name
		    physWorld->GetLogicalVolume(), //its mother  volume
                    false,                         //no boolean operation
                    quarter);                      //copy number
}

//==============================================================================

G4LogicalVolume* TCSDetectorConstruction::GetGDMLVolume(const string file_name,
							const string vol_name) {
  fParser.ReadModule(file_name);
  return fParser.GetVolume(vol_name);
}
