#ifndef TCSHodoscopeConstruction_h
#define TCSHodoscopeConstruction_h 1

#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"

class G4LogicalVolume;
class G4Material;

class TCSHodoscopeConstruction {

public:
  TCSHodoscopeConstruction();
  ~TCSHodoscopeConstruction();
  G4LogicalVolume* GetHodoscope();

  void Construct();

  int GetNCOL() {return fNCOL;};
  int GetNROW() {return fNROW;};

private:

  const int fNCOL = 23;
  const int fNROW = 23;

  const double fModSizeZ = 5.*cm;
  const double fModSizeX = 2.12*cm;   //match calorimeter module + mesh size
  const double fModSizeY = 2.12*cm;

  G4LogicalVolume* fHodoscope;
};

#endif
