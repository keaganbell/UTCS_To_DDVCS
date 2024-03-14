//Muon Detector Construction //
#ifndef TCSMuonDetectorConstruction_h
#define TCSMuonDetectorConstruction_h 1

#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"

class G4LogicalVolume;
class G4Material;

class TCSMuonDetectorConstruction {

public: 
  TCSMuonDetectorConstruction();
  ~TCSMuonDetectorConstruction(); 
  G4LogicalVolume* GetMuon();

  void Construct();


  int GetNCOL() {return fNCOL;};
  int GetNROW() {return fNROW;};
  int GetNZROW() {return fNZROW;};
  
  //just copied from Hodoscope Construction header file 
private:
  
  const int fNCOL = 2;
  const int fNROW = 2;
  const int fNZROW = 4; 
  
  const double fMuonModSizeH = 2.*cm;
  const double fMuonModSizeL = 25.*cm;   //
  const double fMuonModSizeW = 12.*cm;

  //     fMuonModSizeL
  //  ........................
  //  .    scint bar         . fMuonModSizeW
  //  ........................

  // const double fMuonModSizeX = 5.0*cm;   
  //const double fMuonModSizeY = 25.0*cm;
  
  G4LogicalVolume* fMuon;
};

#endif
