#include "TCSHodoscopeConstruction.hh"
//#include "NPSModuleConstruction.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4SystemOfUnits.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
//#include "G4SubtractionSolid.hh"

#include <iostream>

using namespace std;

TCSHodoscopeConstruction::TCSHodoscopeConstruction() {
  //  Construct();  call separately
}

TCSHodoscopeConstruction::~TCSHodoscopeConstruction() {;}

G4LogicalVolume* TCSHodoscopeConstruction::GetHodoscope() {
  return fHodoscope;
}

void TCSHodoscopeConstruction::Construct() {

  G4NistManager* man = G4NistManager::Instance();
  //  man->SetVerbose(1);

  // Module construction.

  G4Material* PS = man->FindOrBuildMaterial("G4_POLYSTYRENE");

  G4Box* modBox = new G4Box("moduleBox",fModSizeX/2,fModSizeY/2,fModSizeZ/2);
  G4LogicalVolume* modLV = new G4LogicalVolume(modBox, PS, "hodo_module_LV",
					       0, 0, 0);

  // Hodoscope construction.

  double xstep = fModSizeX;
  double ystep = fModSizeY;

  double xhodo = xstep*fNCOL;
  double yhodo = ystep*fNROW;
  double zhodo = fModSizeZ;
  cout << "TCSHodoscopeConstruction::GetHodoscope: hodoscope matrix sizes:"
       << endl;
  cout << "  xhodo = " << xhodo/cm << " cm" << endl;
  cout << "  yhodo = " << yhodo/cm << " cm" << endl;
  cout << "  zhodo = " << zhodo/cm << " cm" << endl;

  //  getchar();

  G4Material* Air = man->FindOrBuildMaterial("G4_AIR");

  G4Box* hodoBox = new G4Box("hodoBox", xhodo/2, yhodo/2, zhodo/2);
  fHodoscope = new G4LogicalVolume(hodoBox, Air, "Hodoscope_LV", 0, 0, 0);

  // Positioning of modules in the hodoscope. Numbering of rows from bottom
  // top, columns from left to right.

  int copy_number = 0;
  const double zpos = 0.;

  double ypos = -yhodo/2 + ystep/2;
  for (int irow = 0; irow<fNROW; irow++) {

    double xpos = xhodo/2 - xstep/2;
    for (int icol = 0; icol<fNCOL; icol++) {

      new G4PVPlacement(0,                                //no rotation
			G4ThreeVector(xpos,ypos,zpos),    //at position
			modLV,                       //its logical volume
			"hModule_PV",                    //its name
			fHodoscope,                     //its mother  volume
			false,                            //no boolean operation
			copy_number,                      //copy number
			0);                               //overlaps checking

      copy_number++;
      xpos -= xstep;
    }

    ypos += ystep;
  }

  //  fHodoscope->SetVisAttributes (G4VisAttributes::Invisible);
}
