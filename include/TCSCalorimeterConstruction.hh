#ifndef TCSCalorimeterConstruction_h
#define TCSCalorimeterConstruction_h 1

#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"

class G4LogicalVolume;
class G4Material;

class TCSCalorimeterConstruction {

public:
  TCSCalorimeterConstruction();
  ~TCSCalorimeterConstruction();
  G4LogicalVolume* GetCalorimeter();

  void Construct();

  int GetNCOL() {return fNCOL;};
  int GetNROW() {return fNROW;};

private:

  //const int fNCOL = 23;
  //const int fNROW = 23;

  const int fNCOL = 50;
  const int fNROW = 100;

  //  const double fFrameThick = 1.*mm;    //test
  const double fFrameThick = 0.5*mm;
  const double fFrameWidth = 2.*cm;

  G4LogicalVolume* fCalorimeter;
};

#endif
