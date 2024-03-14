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

#ifndef TCSDetectorConstruction_h
#define TCSDetectorConstruction_h 1

#include "globals.hh"
#include "SimpleField.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SystemOfUnits.hh"
#include <G4SubtractionSolid.hh>

// attempt to use the GDML parser
//#include <vector>
#include "G4GDMLParser.hh"

class G4LogicalVolumeStore;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;
class G4ChordFinder;
class G4UniformMagField;
class G4UserLimits;

/// Detector construction class to define materials and geometry.

class TCSDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  TCSDetectorConstruction();
  ///  TCSDetectorConstruction(G4VPhysicalVolume *aworld)
  ///  {
  ///    G4cout << "We are using GDML..."<<G4endl;
  ///    physWorld =  aworld;
  ///  };
  ////  TCSDetectorConstruction(G4VPhysicalVolume *aworld,
  ////			  G4LogicalVolume *scoring)
  ////  {
  ////    G4cout << "We are using GDML..."<<G4endl;
  ////    physWorld =  aworld;
  ////    fScoringVolume = scoring;
  ////  };
  
  virtual ~TCSDetectorConstruction();

  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();
  G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }
  
protected:
  G4LogicalVolume*  fScoringVolume;

private:

  void ConstructField();
  
  SimpleField *fField;
  ////G4UniformMagField *fField;
  G4Mag_UsualEqRhs *fEquation;
  G4double minStepMagneticField;
  G4MagIntegratorStepper* fStepper;
  G4ChordFinder*          fChordFinder;

  struct {
    ///    const double Distance = 150*cm + 9*cm;       //no carbon fiber mesh
    //const double Distance = 159.2*cm + 10.1*cm;   //NPS like calorimeter

    const double Distance = 359*cm; // need to chamnge after talking with Vardan 
    //For quarters.
    //const int NQuarter = 4;
    //const double PositionAngle = 13.835*degree;
    //const double TiltAngle     = 13.835*degree;
    //const double RotationAngle = 10.034*degree;

    //For single caloromiter facing beam.
        const int NQuarter = 1;
        const double PositionAngle = 0*degree;
        const double TiltAngle     = 0*degree;
        const double RotationAngle = 0*degree;

  } Calo;

  struct {
    const double Distance = 150.*cm; //may29, temp off
    const double PositionAngle = 13.835*degree; //may 29, temp off
    const double TiltAngle     = 13.835*degree; //may 29, temp off
    const double RotationAngle = 10.034*degree; //may 29, temp off
    
    
  } Hodo;

  struct {
    const double PositionAngle = 13.835*degree;
    const double TiltAngle     = 13.835*degree;
    const double RotationAngle = 10.034*degree;
    const double Distance[4] {120*cm, 130*cm, 140*cm, 52*cm};
  } Tracker;

  struct { //DB
    const double Distance = 300.*cm; //z position of muon detector
    const double PositionAngle = 13.835*degree;
    const double TiltAngle     = 13.835*degree;
    const double RotationAngle = 10.034*degree;
  } Muon;
  struct {
    const double Distance = 123.46*cm; // 5 cm between the face of the scattering chamber and the face of the magnet
    const double PositionAngle = 0.0*degree;
    const double mag_vol_x = 46.99*cm;
    const double mag_vol_y = 121.92*cm;
    const double mag_vol_z = 121.92*cm;
  } SBS_mag;
  struct{
    const double len1 = 90.00*cm;
    const double len2 = 37.00*cm;
    const double len3 = 23.00*cm;
    const double len4 = 100.00*cm;
    const double rin_cham_len = (41.0/2*2.54)*cm;
    const double thick_cham_len = (2*2.54)*cm;
    
  }beam_pipe_pos; 

  void PositionCalorimeter(G4LogicalVolume* Calorimeter_log, int quarter);
  void PositionMuon(G4LogicalVolume* Muon_log, int quarter); //DB
  void PositionHodoscope(G4LogicalVolume* Hodoscope_log, int quarter);
  void PositionTracker(G4LogicalVolume* Tracker_log, int quarter, int layer);

  G4LogicalVolume* GetGDMLVolume(const string file_name, const string vol_name);

  // const double fXWorld = 3.*m; temp off may23, 2023
  // const double fYWorld = 3.*m; temp off may23, 2023
  // const double fZWorld = 7.*m; temp off may23, 2023
  
  const double fXWorld = 7.*m; 
  const double fYWorld = 7.*m; 
  const double fZWorld = 20.*m; 

  G4GDMLParser fParser;

  G4VPhysicalVolume* physWorld;

  G4UserLimits* fStepLimit;            // pointer to user step limits

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
