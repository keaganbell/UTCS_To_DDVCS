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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <CLHEP/Units/SystemOfUnits.h>

#include "TCSHistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"

#include "TCSPrimaryGeneratorAction.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSHistoManager::TCSHistoManager() : fKinFile(0), fRootFile(0),
				     fBeamTree(0), fTargetTree(0) , 
				     fCaloTree(0)//, fHodoTree(0), //temp jul14, 2023
				     // fTrackerTree(0) //temp jul14 2023
{
  fKinFileName ="ddvcs_gen.kin_data";
  fRootFileName="ddvcs_setup.root";
  // histogram(s)
  for (G4int k=0; k<MaxHisto; k++) fHisto[k] = 0;
    
  // Trees
  //  fCaloTree = 0;
  //  fHodoXTree = 0;
  //  fHodoYTree = 0;
  //  fTrackerXTree = 0;
  //  fTrackerYTree = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

////TCSHistoManager::TCSHistoManager(char *kname, char *rname) :
TCSHistoManager::TCSHistoManager(string kname, string rname) :
  fKinFile(0), fRootFile(0), fBeamTree(0), fTargetTree(0), fCaloTree(0) //, temp jul14 2023
  //fHodoTree(0), fTrackerTree(0) //temp jul14 2023
{
  fKinFileName = kname;  
  fRootFileName= rname;  

  // histogram(s)
  for (G4int k=0; k<MaxHisto; k++) fHisto[k] = 0;
    
  // Trees
  //  fCaloTree = 0;
  //  fHodoXTree = 0;
  //  fHodoYTree = 0;
  //  fTrackerXTree = 0;
  //  fTrackerYTree = 0;
  //  fKinTree = 0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSHistoManager::~TCSHistoManager()
{
  if ( fRootFile ) delete fRootFile;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSHistoManager::book()
{ 
  /*
  TCSPrimaryGeneratorAction* gen_action;
  enum mode {beam, tcs};
  switch (gen_action->GetMode()) {
  case beam : cout << "TCSHistoManager::book: beam mode is beam " << endl;
    break;
  case tcs  : cout << "TCSHistoManager::book: beam mode is tcs" << endl;
    break;
  default : cout <<  "TCSHistoManager::book: unknown beam mode" << endl;
  }
  */
  
  fKinFile.open(fKinFileName.c_str());  //to read per event kinematic quantities
  if(!fKinFile) {
   G4cout << " ***** HistoManager::book :" 
          << " problem openiing Kin file *****"
          << G4endl;
   return;
  }

  // Creating tree containers to handle histograms and Trees.
  // These trees are associated to an output file.
  //
  //// fRootFile = new TFile(fRootFileName,"RECREATE");
  fRootFile = new TFile(fRootFileName.c_str(),"RECREATE");
  if(!fRootFile) {
    G4cout << " ***** HistoManager::book :" 
	   << " problem creating the ROOT TFile *****"
	   << G4endl;
    return;
  }
   
 fHisto[1] = new TH1D("h1", "Total Edep(MeV)", 100, 0., 5000.);
 if (!fHisto[1]) G4cout << "\n can't create histo 1" << G4endl;

 fBeamTree = new TTree("beam", "TCS beam");
 fBeamTree->Branch("e", &(fBeam.E));
 fBeamTree->Branch("px", &(fBeam.Px));
 fBeamTree->Branch("py", &(fBeam.Py));
 fBeamTree->Branch("pz", &(fBeam.Pz));
 fBeamTree->Branch("x", &(fBeam.X));
 fBeamTree->Branch("y", &(fBeam.Y));
 fBeamTree->Branch("z", &(fBeam.Z));
 fBeamTree->Branch("pid", &(fBeam.PID));

 fTargetTree = new TTree("target", "TCS Target per event hit collections");
 fTargetTree->Branch("edepcont", &(fTargetHitCont.Edep));
 fTargetTree->Branch("pidcont", &(fTargetHitCont.PID));

 fCaloTree = new TTree("calo", "TCS Calorimeters' per event hit collections");
 fCaloTree->Branch("detcont", &(fCaloHitCont.Det));
 fCaloTree->Branch("colcont", &(fCaloHitCont.Col));
 fCaloTree->Branch("rowcont", &(fCaloHitCont.Row));
 fCaloTree->Branch("edepcont", &(fCaloHitCont.Edep));
 fCaloTree->Branch("npecont", &(fCaloHitCont.Npe));
 fCaloTree->Branch("pidcont", &(fCaloHitCont.PID));
 /*
 fHodoTree = new TTree("hodo", "TCS Hodoscopes' per event hit collections");
 fHodoTree->Branch("detcont", &(fHodoHitCont.Det));
 fHodoTree->Branch("colcont", &(fHodoHitCont.Col));
 fHodoTree->Branch("rowcont", &(fHodoHitCont.Row));
 fHodoTree->Branch("edepcont", &(fHodoHitCont.Edep));
 // fHodoTree->Branch("npecont", &(fHodoHitCont.Npe));
 fHodoTree->Branch("pidcont", &(fHodoHitCont.PID));

 fTrackerTree = new TTree("tracker", "TCS trackers' per event hit collections");
 fTrackerTree->Branch("xcont", &(fTrackerHitCont.X));
 fTrackerTree->Branch("ycont", &(fTrackerHitCont.Y));
 fTrackerTree->Branch("edepcont", &(fTrackerHitCont.Edep));
 fTrackerTree->Branch("lengthcont", &(fTrackerHitCont.Length));
 fTrackerTree->Branch("detcont", &(fTrackerHitCont.Det));
 fTrackerTree->Branch("layercont", &(fTrackerHitCont.Layer));
 fTrackerTree->Branch("pidcont", &(fTrackerHitCont.PID));
 fTrackerTree->Branch("pidorigcont", &(fTrackerHitCont.PIDOrig));
 fTrackerTree->Branch("trackidcont", &(fTrackerHitCont.trackID));
 fTrackerTree->Branch("nstepcont", &(fTrackerHitCont.Nstep));
 */ // temp jul 14 2023
 fKinTree = new TTree("kin","TCS kinematics");
 fKinTree->Branch("Q2",&fKinVar.Q2);
 fKinTree->Branch("t",&fKinVar.t);
 fKinTree->Branch("s",&fKinVar.s);
 fKinTree->Branch("xi",&fKinVar.xi);
 fKinTree->Branch("tau",&fKinVar.tau);
 fKinTree->Branch("eta",&fKinVar.eta);
 fKinTree->Branch("phi_cm",&fKinVar.phi_cm);
 fKinTree->Branch("the_cm",&fKinVar.the_cm);
 fKinTree->Branch("psf",&fKinVar.psf);
 // fKinTree->Branch("flux_factor",&fKinVar.flux_factor);
 fKinTree->Branch("FlagSing",&fKinVar.FlagSing);
 fKinTree->Branch("crs_BH",&fKinVar.crs_BH);
 fKinTree->Branch("Eg",&fKinVar.Eg);
 fKinTree->Branch("pminus",fKinVar.P_minus_lab,"pminus[4]/D");
 fKinTree->Branch("pplus",fKinVar.P_plus_lab,"pplus[4]/D");
 fKinTree->Branch("precoil",fKinVar.P_recoil_lab,"precoil[4]/D");
 fKinTree->Branch("vertex",fKinVar.vertexXYZ,"vertex[3]/D");

 G4cout << "\n----> Root file is opened in " << fRootFileName << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSHistoManager::save()
{ 
  if (fRootFile) {
    fRootFile->Write();       // Writing the histograms to the file
    fRootFile->Close();        // and closing the tree (and the file)
    G4cout << "\n----> Histogram Tree is saved \n" << G4endl;
  }
  fKinFile.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSHistoManager::autosave() {
  TDirectory* savedir = gDirectory;
  fRootFile->cd();
  fBeamTree->AutoSave();
  fTargetTree->AutoSave();
  fCaloTree->AutoSave(); // save tree to file
  //fHodoTree->AutoSave(); //temp jul14 2023
  //fTrackerTree->AutoSave(); //temp jul14 2023
  fKinTree->AutoSave();
  fRootFile->SaveSelf();  // save file directory containing this tree
  savedir->cd();
  cout << "TCSHistoManager::autosave: saved data in root file" << endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSHistoManager::FillHisto(G4int ih, G4double xbin, G4double weight)
{
  if (ih >= MaxHisto) {
    G4cout << "---> warning from HistoManager::FillHisto() : histo " << ih
           << " does not exist. (xbin=" << xbin << " weight=" << weight << ")"
           << G4endl;
    return;
  }
 if  (fHisto[ih]) { fHisto[ih]->Fill(xbin, weight); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSHistoManager::Normalize(G4int ih, G4double fac)
{
  if (ih >= MaxHisto) {
    G4cout << "---> warning from HistoManager::Normalize() : histo " << ih
           << " does not exist. (fac=" << fac << ")" << G4endl;
    return;
  }
  if (fHisto[ih]) fHisto[ih]->Scale(fac);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSHistoManager::FillTrees()
{
  if (fBeamTree) {
    //G4cout <<"Filling Tree right now! fHitList size = " << fHitList.size()
    //<< "\n";
    fBeamTree->Fill();
  }

  if (fTargetTree) {
    //G4cout <<"Filling Tree right now! fHitList size = " << fHitList.size()
    //<< "\n";
    fTargetTree->Fill();
  }

  if (fCaloTree) {
    //G4cout <<"Filling Tree right now! fHitList size = " << fHitList.size()
    //<< "\n";
    fCaloTree->Fill();
  }
  /*
  if (fHodoTree) {
    //G4cout <<"Filling Tree right now! fHitList size = " << fHitList.size()
    //<< "\n";
    fHodoTree->Fill();
  }

  if (fTrackerTree) {
    //G4cout <<"Filling Tree right now! fHitList size = " << fHitList.size()
    //<< "\n";
    fTrackerTree->Fill();
  }
  */ // Temp jul14 2023
  fKinFile >> fKinVar.Q2 >> fKinVar.t >> fKinVar.s >> fKinVar.xi >> fKinVar.tau
	   >> fKinVar.eta >> fKinVar.phi_cm >> fKinVar.the_cm
    //	   >> fKinVar.psf >> fKinVar.flux_factor >> fKinVar.crs_BH
	   >> fKinVar.psf >> fKinVar.FlagSing >> fKinVar.crs_BH
	   >> fKinVar.Eg;

  if (fKinTree) {
    fKinTree->Fill();
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSHistoManager::SetBeam(double e, G4ThreeVector p, G4ThreeVector orig,
			      int pid) {
  fBeam.E = e;
  fBeam.Px = p.getX();
  fBeam.Py = p.getY();
  fBeam.Pz = p.getZ();
  fBeam.X = orig.getX();
  fBeam.Y = orig.getY();
  fBeam.Z = orig.getZ();
  fBeam.PID = pid;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSHistoManager::SetTCSVertex(double e, G4ThreeVector p,
				   G4ThreeVector orig, int pid) {

  //  cout << "TCSHistoManager::SetTCSVertex:" << endl;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* electron=particleTable->FindParticle("e-");
  G4ParticleDefinition* positron=particleTable->FindParticle("e+");
  G4ParticleDefinition* proton=particleTable->FindParticle("proton");

  //Pions for background calcualtions.
  G4ParticleDefinition* pi_plus=particleTable->FindParticle("pi+");
  G4ParticleDefinition* pi_minus=particleTable->FindParticle("pi-");

  //  cout << "   Got particle definitions..." << endl;

  int electron_pid = electron->GetPDGEncoding();
  int positron_pid = positron->GetPDGEncoding();
  int proton_pid = proton->GetPDGEncoding();

  int pi_plus_pid = pi_plus->GetPDGEncoding();
  int pi_minus_pid = pi_minus->GetPDGEncoding();

  //  cout << "  electron_pid = " << electron_pid << endl;
  //  cout << "  positron_pid = " << positron_pid << endl;
  //  cout << "  proton_pid   = " << proton_pid << endl;
  //  getchar();

  double P[4] = {e, p[0], p[1], p[2]};
  //  cout << "   P: " << P[0] << " " << P[1] << " " << P[2] << " " << P[3]
  //  << endl;
  //  getchar();

  if (pid == electron_pid || pid == pi_minus_pid)
    copy(P, P+4, fKinVar.P_minus_lab);
  else if (pid == positron_pid || pid == pi_plus_pid)
    copy(P, P+4, fKinVar.P_plus_lab);
  else if (pid == proton_pid)
    copy(P, P+4, fKinVar.P_recoil_lab);

  for (int i=0; i<3; i++)
    fKinVar.vertexXYZ[i] = orig[i];

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSHistoManager::AddHit(double edep, int pid) {

  bool found = false;

  if (!found) {
    fTargetHitCont.Edep.push_back(edep);
    fTargetHitCont.PID.push_back(pid);
    //    cout << "TCSHistoManager::AddHit: pushed back Target edep = " << edep
    //	 << " pid = " << pid << endl;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSHistoManager::AddHit(int det, uint col,uint row, double edep, int pid) {

  bool found = false;

  vector<uint>::iterator ic = fCaloHitCont.Col.begin();
  vector<uint>::iterator ir = fCaloHitCont.Row.begin();
  vector<double>::iterator ie = fCaloHitCont.Edep.begin();
  vector<int>::iterator ip = fCaloHitCont.PID.begin();

  for (vector<int>::iterator id=fCaloHitCont.Det.begin();
       id != fCaloHitCont.Det.end(); id++) {
    //cout << "id = " << *id << " ,det= "<< det <<endl;
    if (*id == det && *ic == col && *ir == row && *ip == pid) {
      //      cout << "AddHit: *ie = " << *ie << "  edep = " << edep << endl;
        *ie += edep;
      //      cout << "AddHit: *ie = " << *ie << endl;
      //      getchar();
      found = true;
      break;
    }

    ic++; ir++; ie++; ip++;
  }

  if (!found) {
    fCaloHitCont.Det.push_back(det);
    fCaloHitCont.Col.push_back(col);
    fCaloHitCont.Row.push_back(row);
    fCaloHitCont.Edep.push_back(edep);
    fCaloHitCont.Npe.push_back(0);
    fCaloHitCont.PID.push_back(pid);
    //cout << "Pushed back cal hit " << det << " " << col << " " << row << " "
    //<< edep << endl;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSHistoManager::AddHit(int det, uint col,uint row, int npe, int pid) {

  bool found = false;

  vector<uint>::iterator ic = fCaloHitCont.Col.begin();
  vector<uint>::iterator ir = fCaloHitCont.Row.begin();
  vector<int>::iterator in = fCaloHitCont.Npe.begin();
  vector<int>::iterator ip = fCaloHitCont.PID.begin();
 
  for (vector<int>::iterator id=fCaloHitCont.Det.begin();
       id != fCaloHitCont.Det.end(); id++) {
   
    if (*id == det && *ic == col && *ir == row && *ip == pid) {
      //      cout << "AddHit: *ie = " << *ie << "  npe = " << npe << endl;
     
      *in += npe;
      //      cout << "AddHit: *ie = " << *ie << endl;
      //      getchar();
      found = true;
      break;
    }

    ic++; ir++; in++; ip++;
  }

  if (!found) {
    fCaloHitCont.Det.push_back(det);
    fCaloHitCont.Col.push_back(col);
    fCaloHitCont.Row.push_back(row);
    fCaloHitCont.Edep.push_back(0.);
    fCaloHitCont.Npe.push_back(npe);
    fCaloHitCont.PID.push_back(pid);
    //cout << "Pushed back cal hit " << det << " " << col << " " << row << " "
    //<< npe << endl;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
void TCSHistoManager::AddHit(int det, uint chan, double edep, int pid,
                             HodoHitContainer& HodoHitCont) {

  // Add hit to a hodoscope hit container.

  bool found = false;

  vector<uint>::iterator ic = HodoHitCont.Chan.begin();
  vector<double>::iterator ie = HodoHitCont.Edep.begin();
  vector<int>::iterator ip = HodoHitCont.PID.begin();

  for (vector<int>::iterator id=HodoHitCont.Det.begin();
       id != HodoHitCont.Det.end(); id++) {

    if (*id == det && *ic == chan && *ip == pid) {
      //      cout << "          TCSHistoManager::AddHit: *ie = " << *ie
      //           << "  edep = " << edep << endl;
      *ie += edep;
      //    cout << "          TCSHistoManager::AddHit: *ie = " << *ie << endl;
      //      getchar();
      found = true;
      break;
    }

    ic++; ie++; ip++;
  }

  if (!found) {
    HodoHitCont.Det.push_back(det);
    HodoHitCont.Chan.push_back(chan);
    HodoHitCont.Edep.push_back(edep);
    HodoHitCont.PID.push_back(pid);
    //G4cout << "          TCSHistoManager::AddHit: hit pushed back" << G4endl;
  }

}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSHistoManager::AddHit(double x, double y, double edep, double length,
			     int det, int layer, int pid,
			     int pidorig, int trackid,
                             TrackerHitContainer& TrackerHitCont) {

  // Add hit to a tracker hit container.

  bool new_track = true;

  for (uint i=0; i<TrackerHitCont.trackID.size(); i++) {

    if (pidorig == TrackerHitCont.PIDOrig.at(i) &&
	trackid == TrackerHitCont.trackID.at(i) &&
	det     == TrackerHitCont.Det.at(i) &&
	layer   == TrackerHitCont.Layer.at(i)) {
  
      new_track = false;

      double edep_old = TrackerHitCont.Edep.at(i);
      double length_old = TrackerHitCont.Length.at(i);
      double x_old = TrackerHitCont.X.at(i);
      double y_old = TrackerHitCont.Y.at(i);
      double edep_new = edep_old + edep;
      double length_new = length_old + length;

      if (edep_new > 0.) {
	double x_new = (x_old*edep_old + x*edep)/edep_new;
	double y_new = (y_old*edep_old + y*edep)/edep_new;
	TrackerHitCont.X.at(i) = x_new;
	TrackerHitCont.Y.at(i) = y_new;
      }

      TrackerHitCont.Edep.at(i) = edep_new;
      TrackerHitCont.Length.at(i) = length_new;
      TrackerHitCont.Nstep.at(i)++;

      if (pid != TrackerHitCont.PID.at(i))
	G4cout << "*** TCSHistoManager::AddHit: pid != TrackerHitCont.PID ***"
	       << G4endl;
      break;
    }

  }

  if (new_track) {
      TrackerHitCont.X.push_back(x);
      TrackerHitCont.Y.push_back(y);
      TrackerHitCont.Edep.push_back(edep);
      TrackerHitCont.Length.push_back(length);
      TrackerHitCont.Det.push_back(det);
      TrackerHitCont.Layer.push_back(layer);
      TrackerHitCont.PID.push_back(pid);
      TrackerHitCont.PIDOrig.push_back(pidorig);
      TrackerHitCont.trackID.push_back(trackid);
      TrackerHitCont.Nstep.push_back(1);
  }

  //G4cout << "          TCSHistoManager::AddHit: hit pushed back" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSHistoManager::AddHodoHit(int det, uint col,uint row, double edep,
				 int pid) {

  bool found = false;

  vector<uint>::iterator ic = fHodoHitCont.Col.begin();
  vector<uint>::iterator ir = fHodoHitCont.Row.begin();
  vector<double>::iterator ie = fHodoHitCont.Edep.begin();
  vector<int>::iterator ip = fHodoHitCont.PID.begin();

  for (vector<int>::iterator id=fHodoHitCont.Det.begin();
       id != fHodoHitCont.Det.end(); id++) {

    if (*id == det && *ic == col && *ir == row && *ip == pid) {
      //      cout << "AddHit: *ie = " << *ie << "  edep = " << edep << endl;
      *ie += edep;
      //      cout << "AddHit: *ie = " << *ie << endl;
      //      getchar();
      found = true;
      break;
    }

    ic++; ir++; ie++; ip++;
  }

  if (!found) {
    fHodoHitCont.Det.push_back(det);
    fHodoHitCont.Col.push_back(col);
    fHodoHitCont.Row.push_back(row);
    fHodoHitCont.Edep.push_back(edep);
    //    fHodoHitCont.Npe.push_back(0);
    fHodoHitCont.PID.push_back(pid);
    //cout << "Pushed back cal hit " << det << " " << col << " " << row << " "
    //<< edep << endl;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSHistoManager::PrintStatistic()
{
  if(fHisto[1]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    /*
    G4cout 
    << " EAbs : mean = " << G4BestUnit(fHisto[1]->GetMean(), "Energy") 
    << " rms = " << G4BestUnit(fHisto[1]->GetRMS(),  "Energy") << G4endl;
    */
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
