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

#ifndef TCSHistoManager_h
#define TCSHistoManager_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"

#include "G4SystemOfUnits.hh"    
#include <vector>
#include <fstream>

using namespace std;

class TFile;
class TTree;
class TH1D;

struct TrackerHitContainer {
  vector<double> X;
  vector<double> Y;
  vector<double> Edep;
  vector<double> Length;
  vector<int> Det;
  vector<int> Layer;
  vector<int> PID;
  vector<int> PIDOrig;
  vector<int> trackID;
  vector<int> Nstep;
};

struct TargetHitContainer {
  vector<double> Edep;
  vector<int> PID;
};

const G4int MaxHisto = 5;

class TCSHistoManager {

public:
  
  TCSHistoManager();
  TCSHistoManager(string,string);
  ~TCSHistoManager();
   
  void book();
  void save();
  void autosave();    //save intermittent data.

  void FillHisto(G4int id, G4double bin, G4double weight = 1.0);
  void Normalize(G4int id, G4double fac);    

  void SetBeam(double e, G4ThreeVector p, G4ThreeVector orig, int pid);
  void SetTCSVertex(double e, G4ThreeVector p, G4ThreeVector orig, int pid);
  void AddHit(double edep, int pid);
  void AddHit(int det, uint col, uint row, double edep, int pid);
  void AddHit(int det, uint col, uint row, int npe, int pid);
  //  void AddHit(int det, uint chan, double edep, int pid,
  //	      HodoHitContainer &HodoHitCont);
  void AddHit(double x, double y, double edep, double length,
	      int det, int layer, int pid, int pidorig, int trackid,
	      TrackerHitContainer &TrackerHitCont);
  //For Hodo-s, to be changed by AddHit later on.
  void AddHodoHit(int det, uint col, uint row, double edep, int pid);

  bool CheckTargetHitCont() {
    return (fTargetHitCont.Edep.size() == fTargetHitCont.PID.size());
  }

  bool CheckCaloHitCont() {
    uint sz = fCaloHitCont.Det.size(); 
    return (fCaloHitCont.Col.size() != sz || fCaloHitCont.Row.size() != sz ||
	    fCaloHitCont.Edep.size() != sz || fCaloHitCont.Npe.size() != sz ||
	    fCaloHitCont.PID.size() != sz
	    ? false : true);
  }

  bool CheckHodoHitCont() {
    uint sz = fHodoHitCont.Det.size();
    return (fHodoHitCont.Col.size() != sz || fHodoHitCont.Row.size() != sz ||
	    fHodoHitCont.Edep.size() != sz ||
	    //fHodoHitCont.Npe.size() != sz ||
	    fHodoHitCont.PID.size() != sz
	    ? false : true);
  }

  bool CheckTrackerHitCont(TrackerHitContainer &TrackerHitCont) {
    uint sz = TrackerHitCont.Det.size();
    return (TrackerHitCont.Nstep.size() != sz ||
	    TrackerHitCont.Layer.size() != sz ||
	    TrackerHitCont.PID.size() != sz ||
	    TrackerHitCont.PIDOrig.size() != sz ||
	    TrackerHitCont.trackID.size() != sz ||
	    TrackerHitCont.X.size() != sz ||
	    TrackerHitCont.Y.size() != sz ||
	    TrackerHitCont.Edep.size() != sz ||
	    TrackerHitCont.Length.size() != sz
	    ? false : true);
  }

  void ResetBeam() {
    fBeam.E =0;
    fBeam.Px = 0.;
    fBeam.Py = 0.;
    fBeam.Pz = 0.;
    fBeam.X = 0.;
    fBeam.Y = 0.;
    fBeam.Z = 0.;
    fBeam.PID = 0;
  };

  void ResetTarget() {
    fTargetHitCont.Edep.clear();
    fTargetHitCont.PID.clear();
  };

  void ResetCalo() {
    fCaloHitCont.Det.clear();
    fCaloHitCont.Col.clear();
    fCaloHitCont.Row.clear();
    fCaloHitCont.Edep.clear();
    fCaloHitCont.Npe.clear();
    fCaloHitCont.PID.clear();
  };

  void ResetHodo() {
    fHodoHitCont.Det.clear();
    fHodoHitCont.Col.clear();
    fHodoHitCont.Row.clear();
    fHodoHitCont.Edep.clear();
    //    fHodoHitCont.Npe.clear();
    fHodoHitCont.PID.clear();
  };

  void ResetTracker(TrackerHitContainer &TrackerHitCont) {
    TrackerHitCont.X.clear();
    TrackerHitCont.Y.clear();
    TrackerHitCont.Edep.clear();
    TrackerHitCont.Length.clear();
    TrackerHitCont.Det.clear();
    TrackerHitCont.Layer.clear();
    TrackerHitCont.PID.clear();
    TrackerHitCont.PIDOrig.clear();
    TrackerHitCont.trackID.clear();
    TrackerHitCont.Nstep.clear();
  };

  void ResetKinVar() {
    fKinVar.Q2 = 0.;
    fKinVar.t = 0.;
    fKinVar.s = 0.;
    fKinVar.xi = 0.;
    fKinVar.tau = 0.;
    fKinVar.eta = 0.;
    fKinVar.phi_cm = 0.;
    fKinVar.the_cm = 0.;
    fKinVar.psf = 0.;
    //    fKinVar.flux_factor = 0.;
    fKinVar.FlagSing = 0;
    fKinVar.crs_BH = 0.;
    fKinVar.Eg = 0.;
    for (int i=0; i<4; i++) {
      fKinVar.P_minus_lab[i] = 0.;
      fKinVar.P_plus_lab[i] = 0.;
      fKinVar.P_recoil_lab[i] = 0.;
    };
    for (int i=0; i<3; i++)
      fKinVar.vertexXYZ[i] = 0.;
  }

  void Reset() {
    ResetBeam();
    ResetTarget();
    ResetCalo();
    ResetHodo();
    ResetTracker(fTrackerHitCont);
    ResetKinVar();
  };

  void FillTrees();

  void PrintStatistic();
        
private:

  ////  char*    fRootFileName;
  ////  char*    fKinFileName;
  string    fKinFileName;
  string    fRootFileName;

  ifstream fKinFile;
  TFile*   fRootFile;

  TTree*   fBeamTree;
  TTree*   fTargetTree;
  TTree*   fCaloTree;
  TTree*   fHodoTree;
  TTree*   fTrackerTree;
  TTree*   fKinTree;

  TH1D*    fHisto[MaxHisto];            

  struct {
    double E;
    double Px, Py, Pz;
    double X, Y, Z;
    int PID;
  } fBeam;

  TargetHitContainer fTargetHitCont;

  struct CaloHitContainer {
    vector<int> Det;
    vector<uint> Col;
    vector<uint> Row;
    vector<double> Edep;
    vector<int> Npe;
    vector<int> PID;
  } fCaloHitCont;

  struct HodoHitContainer {
    vector<int> Det;
    vector<uint> Col;
    vector<uint> Row;
    vector<double> Edep;
    //    vector<int> Npe;
    vector<int> PID;
  } fHodoHitCont;

  TrackerHitContainer fTrackerHitCont;

  struct KinVar {
    double Q2;
    double t;
    double s;
    double xi;
    double tau;
    double eta;
    double phi_cm;
    double the_cm;
    double psf;
    //    double flux_factor;
    int FlagSing;   //BH singularity flag
    double crs_BH;
    double Eg;
    double P_minus_lab[4];   //momenta at the TCS vertex
    double P_plus_lab[4];
    double P_recoil_lab[4];
    double vertexXYZ[3];
  } fKinVar;

  friend class TCSEventAction;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
