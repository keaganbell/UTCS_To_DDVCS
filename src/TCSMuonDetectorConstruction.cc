//Muon Detector Construction// 

#include "TCSMuonDetectorConstruction.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4SystemOfUnits.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"

#include <iostream>

using namespace std; 

TCSMuonDetectorConstruction::TCSMuonDetectorConstruction(){

}

TCSMuonDetectorConstruction::~TCSMuonDetectorConstruction() {;}

G4LogicalVolume* TCSMuonDetectorConstruction::GetMuon() {
  return fMuon;
}

void TCSMuonDetectorConstruction::Construct(){

 //Materials can also be defined using the internal Geant4 database, based on NIST
 //NIST : National Institute of Standards and Technologies
// G4NistManager* man = G4NistManager::Instance();

cout <<"------------------------------------------------"<<endl;
cout <<"                                                 "<<endl;
cout <<"check the muon detector construction class" << endl;
cout <<"                                                 "<<endl;
//Get nist material manager
 G4NistManager* man = G4NistManager::Instance();

 //module construction 

 G4Material* PS = man->FindOrBuildMaterial("G4_POLYSTYRENE");
 
 G4Box* modBox = new G4Box("moduleBox", fMuonModSizeL,fMuonModSizeW,fMuonModSizeH); //size of the each module


 cout << "muon module size =" << endl; 
 cout << fMuonModSizeL <<" by "<< fMuonModSizeW<< " by " <<fMuonModSizeH <<  endl; // this will print out in mm

 G4LogicalVolume* modLV = new G4LogicalVolume(modBox, PS, "muon_module_LV",

					       0, 0, 0);

 //90 degree rot

 // G4Material* PS1 = man->FindOrBuildMaterial("G4_POLYSTYRENE");
 
 G4Box* modBox1 = new G4Box("moduleBox1", fMuonModSizeW,fMuonModSizeL,fMuonModSizeH); //size of the each module

 G4LogicalVolume* modLV1 = new G4LogicalVolume(modBox1, PS, "muon_module_LV1",

					       0, 0, 0);
  
 //define the size of the big box

  double xmuon = fMuonModSizeL;
  double ymuon = fMuonModSizeL;
  double zmuon = fNZROW*fMuonModSizeH; //size of box

  cout << "TCSMuonscopeConstruction::GetMuonscope: muonscope matrix sizes:"
       << endl;
  cout << "  xmuon = " << xmuon/cm << " cm" << endl; // "xmuon/cm" means show xmuon in cm 
  cout << "  ymuon = " << ymuon/cm << " cm" << endl;
  cout << "  zmuon = " << zmuon/cm << " cm" << endl;


 G4Material* Air = man->FindOrBuildMaterial("G4_AIR");
 
 //G4Box* muonBox = new G4Box("muonBox", xmuon/2, ymuon/2, zmuon/2);
 G4Box* muonBox = new G4Box("muonBox", xmuon, ymuon, zmuon); //size of bigger muon box 
 
 fMuon = new G4LogicalVolume(muonBox, Air, "Muon_LV", 0, 0, 0);


  //Muon Detector construction.

  double xstep = fMuonModSizeL;
  double ystep = fMuonModSizeL;
  double zstep = fMuonModSizeH;



  
 // this is for parts inside the big volume 
 
 // Positioning of modules in the muon detector. Numbering of rows from bottom
 // top, columns from left to right.
 
 // double zpos = -zmuon/2 + zstep/2;
  double zpos = -60;

 
 for (int izrow =0; izrow<fNZROW; izrow++) {
    cout << "zpos = " << zpos/cm << endl;
 

   if(izrow %2 == 0  ){

     double ypos = ymuon/2 - ystep;
  
     //double xpos = xmuon/2 - xstep;
     double xpos = 0;
     
     for (int irow =0; irow<fNROW; irow++) {
       
       new G4PVPlacement(0,                              //no rotation of mother frame 
			 G4ThreeVector(xpos, ypos,zpos), //at position of inside pannels
			 modLV,                          //its logical volume
			 "mModule_PV",                   //its name
			 fMuon,                          //its mother  volume
			 false,                          //no boolean operation
			 0,                              //copy number : 0 now, check later
			 0);                             //overlap check 
       
       ypos += ystep;
      
     }
   } else{
     
     double xpos = xmuon/2 - xstep;
     
     // double ypos = xmuon/2 - xstep;
     double ypos =0;

        for (int icol = 0; icol<fNCOL; icol++) {
     
     
           new G4PVPlacement(0,                              //no rotation
      			 G4ThreeVector(xpos, ypos,zpos), //at position of inside pannels
      			 modLV1,                          //its logical volume
      			 "mModule_PV1",                   //its name
      			 fMuon,                          //its mother  volume
      			 false,                          //no boolean operation
      			 0,                              //copy number : 0 now, check later
      			 0);                             //overlaps checking
          xpos += xstep;
       }
   } 
   
   zpos += (zstep*2); 
   
 }
 
 
 //G4Material* Air = man->FindOrBuildMaterial("G4_AIR");
 
 // G4Box* muonBox = new G4Box("muonBox", xcase/2, ycase/2, zcase/2);
 //G4Box* muonBox = new G4Box("muonBox", 52.31, 57.31, 35.227);
 
 //fMuon = new G4LogicalVolume(muonBox, Air, "muon_LV", 0, 0, 0);
 
 
 cout <<"------------------------------------------------"<<endl;
 
 
 //fMuonDetector = new G4LogicalVolume
}

 
