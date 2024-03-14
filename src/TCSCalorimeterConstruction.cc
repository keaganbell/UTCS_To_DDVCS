#include "TCSCalorimeterConstruction.hh"
#include "NPSModuleConstruction.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4SystemOfUnits.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
//#include "G4MultiUnion.hh"
#include "G4SubtractionSolid.hh"

#include <iostream>

using namespace std;

TCSCalorimeterConstruction::TCSCalorimeterConstruction() {
  //  Construct();  call separately
}

TCSCalorimeterConstruction::~TCSCalorimeterConstruction() {;}

G4LogicalVolume* TCSCalorimeterConstruction::GetCalorimeter() {
  return fCalorimeter;
}

void TCSCalorimeterConstruction::Construct() {

  G4NistManager* man = G4NistManager::Instance();
  //  man->SetVerbose(1);

  // Module construction and its sizes.

  NPSModuleConstruction* ModuleConstruction = new NPSModuleConstruction(man);

  double xmod = ModuleConstruction->GetSizeX();
  double ymod = ModuleConstruction->GetSizeY();
  double zmod = ModuleConstruction->GetSizeZ();

  cout << "TCSCalorimeterConstruction::GetCalorimeter: module sizes:" << endl;
  cout << "  xmod = " << xmod/cm << " cm" << endl;
  cout << "  ymod = " << ymod/cm << " cm" << endl;
  cout << "  zmod = " << zmod/cm << " cm" << endl;

  // Calorimeter construction.

  double xstep = xmod + fFrameThick;
  double ystep = ymod + fFrameThick;

  double xcal = xstep*fNCOL + fFrameThick;
  double ycal = ystep*fNROW + fFrameThick;
  double zcal = zmod;
  cout << "TCSCalorimeterConstruction::GetCalorimeter: calo matrix sizes:"
       << endl;
  cout << "  xcal = " << xcal/cm << " cm" << endl;
  cout << "  ycal = " << ycal/cm << " cm" << endl;
  cout << "  zcal = " << zcal/cm << " cm" << endl;

  //Calorimeter's case sizes.
  double xcase = xcal + 3.5*cm;
  double ycase = ycal + 8.5*cm;
  double zcase = zcal + 12.*cm;
  cout << "                              Calorimeter's case sizes:" << endl;
  cout << "  xcase = " << xcase/cm << " cm" << endl;
  cout << "  ycase = " << ycase/cm << " cm" << endl;
  cout << "  zcase = " << zcase/cm << " cm" << endl;

  //  getchar();

  G4Material* Air = man->FindOrBuildMaterial("G4_AIR");

  G4Box* caloBox = new G4Box("caloBox", xcase/2, ycase/2, zcase/2);
  fCalorimeter = new G4LogicalVolume(caloBox, Air, "Calorimeter_LV", 0, 0, 0);

  // Positioning of modules in the calorimeter. Numbering of rows from bottom
  // top, columns from left to right.

  int copy_number = 0;
  const double zpos = 0.;

  double ypos = -ycal/2 + ystep/2;
  for (int irow = 0; irow<fNROW; irow++) {

    ///    double xpos = -xcal/2 + xstep/2;
    double xpos = xcal/2 - xstep/2;
    for (int icol = 0; icol<fNCOL; icol++) {
      G4LogicalVolume* Module_log = ModuleConstruction->GetModule();

      new G4PVPlacement(0,                                //no rotation
			G4ThreeVector(xpos,ypos,zpos),    //at position
			Module_log,                       //its logical volume
			"Module_phys",                    //its name
			fCalorimeter,                     //its mother  volume
			false,                            //no boolean operation
			copy_number,                      //copy number
			0);                               //overlaps checking

      copy_number++;
      ///      xpos += xstep;
      xpos -= xstep;
    }

    ypos += ystep;
  }

  // Support frames, made from carbon composite (carbon fiber + epoxy,
  // front frame) and mu-metal (rear frame).

  G4Element* C  = man->FindOrBuildElement("C");
  G4Element* O  = man->FindOrBuildElement("O");
  G4Element* H  = man->FindOrBuildElement("H");

  //Carbon fiber of average density.
  G4Material* carbon_fiber = new G4Material("Carbon fiber", 1.84*g/cm3, 1);
  carbon_fiber->AddElement(C, 1);

  G4Element* Ni  = man->FindOrBuildElement("Ni");
  G4Element* Mo  = man->FindOrBuildElement("Mo");
  G4Element* Fe  = man->FindOrBuildElement("Fe");

  //mu-metal (ASTM A753 Alloy 4) according to mu-metal.com.
  G4Material* mu_metal = new G4Material("mu-metal", 8.7*g/cm3, 3);
  mu_metal->AddElement(Ni, 80.*perCent);
  mu_metal->AddElement(Mo,  5.*perCent);
  mu_metal->AddElement(Fe, 15.*perCent);

  //Bisphenol A epoxy resin, with typical epoxy density.
  G4Material* epoxy = new G4Material("Bisphenol A", 1.25*g/cm3, 3);
  epoxy->AddElement(C, 17);
  epoxy->AddElement(O,  3);
  epoxy->AddElement(H, 19);

  G4Material* carbon_composite = new G4Material("Carbon composite",1.4*g/cm3,2);
  carbon_composite->AddMaterial(carbon_fiber, 50.*perCent);
  carbon_composite->AddMaterial(epoxy,        50.*perCent);

  /*
  G4Box* cellBox = new G4Box("cellBox", xmod/2, ymod/2, fFrameWidth/2+0.5*mm);
  G4MultiUnion* cellsSolid = new G4MultiUnion("cellsSolid");

  ypos = -ycal/2 + ystep/2;
  for (int irow = 0; irow<fNROW; irow++) {

    double xpos = -xcal/2 + xstep/2;
    for (int icol = 0; icol<fNCOL; icol++) {

      G4RotationMatrix rotm  = G4RotationMatrix();
      G4ThreeVector pos = G4ThreeVector(xpos, ypos, zpos);
      G4Transform3D tr  = G4Transform3D(rotm, pos);

      cellsSolid->AddNode(*cellBox, tr);

      xpos += xstep;
    }

    ypos += ystep;
  }

  cellsSolid->Voxelize();

  G4Box* frameBox = new G4Box("frameBox", xcal/2, ycal/2, fFrameWidth/2);

  G4SubtractionSolid* frameSolid =
    new G4SubtractionSolid("frameSolid", frameBox, cellsSolid);

  G4LogicalVolume* frameLog = new G4LogicalVolume(frameSolid, carbon_composite,
						  "frame", 0, 0, 0);

  new G4PVPlacement(0,                                          //no rotation
		    G4ThreeVector(0.,0.,-zcal/2+fFrameWidth/2), //at position
		    frameLog,                                   //logical volume
		    "front frame",                              //its name
		    fCalorimeter,                               //mother  volume
		    false,                                      //no bool. op.
		    0,                                          //copy number
		    0);                                         //overlaps check

  double zpos_frame = -zcal/2+ModuleConstruction->GetBlockSizeZ()-fFrameWidth/2;
  new G4PVPlacement(0,                                          //no rotation
		    G4ThreeVector(0.,0.,zpos_frame),            //at position
		    frameLog,                                   //logical volume
		    "back frame",                               //its name
		    fCalorimeter,                               //mother  volume
		    false,                                      //no bool. op.
		    0,                                          //copy number
		    0);                                         //overlaps check
  */

  G4Box* voidBox = new G4Box("voidBox", xmod/2, ymod/2, fFrameWidth/2+0.5*mm);
  G4Box* frameBox = new G4Box("frameBox", (xmod+fFrameThick)/2,
			                  (ymod+fFrameThick)/2,
			                  fFrameWidth/2);
  G4SubtractionSolid* frameSolid =
    new G4SubtractionSolid("frameSolid", frameBox, voidBox);

  G4LogicalVolume* FrontFrameLog = new G4LogicalVolume(frameSolid,
				   carbon_composite, "front frame", 0, 0, 0);

  G4LogicalVolume* RearFrameLog = new G4LogicalVolume(frameSolid,
				  mu_metal, "rear frame", 0, 0, 0);

  copy_number = 0;

  ypos = -ycal/2 + ystep/2;
  for (int irow = 0; irow<fNROW; irow++) {

    ///    double xpos = -xcal/2 + xstep/2;
    double xpos = xcal/2 - xstep/2;
    for (int icol = 0; icol<fNCOL; icol++) {
  
      new G4PVPlacement(0,                                //no rotation
			G4ThreeVector(xpos,ypos,-zcal/2+fFrameWidth/2), //at pos
			FrontFrameLog,                    //its logical volume
			"front frame",                    //its name
			fCalorimeter,                     //its mother  volume
			false,                            //no boolean operation
			copy_number,                      //copy number
			0);                               //overlaps checking

      double zpos_frame =
	     -zcal/2+ModuleConstruction->GetBlockSizeZ()-fFrameWidth/2;
      new G4PVPlacement(0,                                //no rotation
			G4ThreeVector(xpos,ypos,zpos_frame), //at pos
			RearFrameLog,                     //its logical volume
			"rear frame",                     //its name
			fCalorimeter,                     //its mother  volume
			false,                            //no boolean operation
			copy_number,                      //copy number
			0);                               //overlaps checking

      copy_number++;
      ///      xpos += xstep;
      xpos -= xstep;
    }

    ypos += ystep;
  }

  // PE front plate for sensosrs. ----------------------------------------------

  ///  const double sensor_plateThick = 15.*mm;

  ///  G4Box* sensor_plateBox = new G4Box("sensor_plateBox", xcal/2., ycal/2.,
  ///				     sensor_plateThick/2.);

  ///  G4Material* PE = man->FindOrBuildMaterial("G4_POLYETHYLENE");

  ///  G4LogicalVolume* sensor_plateLog = new G4LogicalVolume(sensor_plateBox,
  ///					 PE, "sensor_plate_log", 0, 0, 0);

  ///  double sensor_plate_Z = -zcal/2. - 1.*cm - sensor_plateThick/2.;

  ///  new G4PVPlacement(0,                                //no rotation
  ///		    G4ThreeVector(0.,0.,sensor_plate_Z), //at pos
  ///		    sensor_plateLog,                  //its logical volume
  ///		    "SensorPlate phys",               //its name
  ///		    fCalorimeter,                     //its mother  volume
  ///		    false,                            //no boolean operation
  ///		    0,                                //copy number
  ///		    0);                               //overlaps checking

  //  G4Colour yellow(0.5, 0.5, 0.);
  //  G4VisAttributes*PEVisAttributes= new G4VisAttributes(yellow);
  //  sensor_plateLog->SetVisAttributes(PEVisAttributes);

  // Cooling plate sizes. ------------------------------------------------------

  const double cool_plate_dX =  12.*mm;
  const double cool_plate_dZ = 165.*mm;

  // Front Aluminum frame. -----------------------------------------------------

  //  const double front_frame_dZ = 20.*mm;
  //  const double front_frame_dY = 36.*mm;

  ///  G4Box* front_frameBox = new G4Box("front_frameBox",
  ///				    //(2.*front_frame_dY + xcal)/2.,
  ///				    (2.*cool_plate_dX + xcal)/2.,
  ///				    (2.*front_frame_dY + ycal)/2.,
  ///				    front_frame_dZ/2.);
  ///G4Box* front_frameVoid = new G4Box("front_frameVoid",xcal/2.,ycal/2.,zcal/2.);
  ///  G4SubtractionSolid* front_frameSolid = new G4SubtractionSolid(
  ///		      "front_frameSolid", front_frameBox, front_frameVoid);

  G4Material* Al = man->FindOrBuildMaterial("G4_Al");

  ///G4LogicalVolume* front_frameLog = new G4LogicalVolume(front_frameSolid, Al,
  ///"front_frameLog",0,0,0);

  //Thickness of platic frame (to hold front reflector)
  ///  const double deltaZ = 2.5*mm;
  ///  const double front_frame_Z = -zcal/2. - deltaZ + front_frame_dZ/2.;

  ///  new G4PVPlacement(0,                                //no rotation
  ///		    G4ThreeVector(0.,0.,front_frame_Z), //at pos
  ///		    front_frameLog,                   //its logical volume
  ///		    "FrontFrame phys",                //its name
  ///		    fCalorimeter,                     //its mother  volume
  ///		    false,                            //no boolean operation
  ///		    0,                                //copy number
  ///		    0);                               //overlaps checking

  // Cooling plates. -----------------------------------------------------------

  G4Box* cool_plateBox = new G4Box("cool_plateBox",
				   (xcal+cool_plate_dX)/2.,
				   (ycal+cool_plate_dX)/2.,
				   cool_plate_dZ/2.);
  G4Box* cool_plateVoid  = new G4Box("cool_plateVoid",xcal/2.,ycal/2.,zcal/2.);

  G4SubtractionSolid* cool_plateSolid =
    new G4SubtractionSolid("cool_plateSolid", cool_plateBox, cool_plateVoid);

  G4Material* Cu  = man->FindOrBuildMaterial("G4_Cu");

  G4LogicalVolume* cool_plateLog = new G4LogicalVolume(cool_plateSolid,
				       Cu, "cool_plate_log", 0,0,0);

  ///double cool_plate_Z = front_frame_Z + front_frame_dZ/2. + cool_plate_dZ/2.;

  ///  new G4PVPlacement(0,                                //no rotation
  ///		    G4ThreeVector(0.,0.,cool_plate_Z), //at pos
  ///		    cool_plateLog,                    //its logical volume
  ///		    "CoolingPlates phys",             //its name
  ///		    fCalorimeter,                     //its mother  volume
  ///		    false,                            //no boolean operation
  ///		    0,                                //copy number
  ///		    0);                               //overlaps checking

  G4Colour brown(0.7, 0.4, 0.1);
  G4VisAttributes*copperVisAttributes= new G4VisAttributes(brown);
  cool_plateLog->SetVisAttributes(copperVisAttributes);

  // Rear aluminum frame, assume identical to the Front frame.------------------

  /// double rear_frame_Z = cool_plate_Z + cool_plate_dZ/2. + front_frame_dZ/2.;

  ///  new G4PVPlacement(0,                                //no rotation
  ///		    G4ThreeVector(0.,0.,rear_frame_Z), //at pos
  ///		    front_frameLog,                   //its logical volume
  ///		    "RearFrame phys",                //its name
  ///		    fCalorimeter,                     //its mother  volume
  ///		    false,                            //no boolean operation
  ///		    0,                                //copy number
  ///		    0);                               //overlaps checking

  // Case. ---------------------------------------------------------------------

  //  G4Material* Al = man->FindOrBuildMaterial("G4_Al");
  ///  G4Material* PVC = man->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");

  G4Box* caseBox = new G4Box("caseBox", xcase/2, ycase/2, zcase/2);

  const double caseThick = 1.*mm;
  //  const double caseThick = 10.*mm;   //10 or 15 mm (from Carlos)
  //  const double caseThick = 5.*mm;

  G4Box* caseVoid = new G4Box("caseVoid", xcase/2-caseThick, ycase/2-caseThick,
			      zcase/2-caseThick);
  G4SubtractionSolid* caseSolid =
    new G4SubtractionSolid("caseSolid", caseBox, caseVoid);

  ///  G4LogicalVolume* caseLog = new G4LogicalVolume(caseSolid, PVC, "case_log",
  G4LogicalVolume* caseLog = new G4LogicalVolume(caseSolid, Al, "case_log",
						 0, 0, 0);

  new G4PVPlacement(0,                                //no rotation
		    G4ThreeVector(),                  //at pos
		    caseLog,                          //its logical volume
		    "CaloCase phys",                  //its name
		    fCalorimeter,                     //its mother  volume
		    false,                            //no boolean operation
		    0,                                //copy number
		    0);                               //overlaps checking

  ///  fCalorimeter->SetVisAttributes (G4VisAttributes::Invisible);
}
