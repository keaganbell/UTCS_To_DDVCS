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

#include "TCSEventAction.hh"
#include "TCSRun.hh"
#include "TCSRunAction.hh"
#include "TCSHistoManager.hh"
#include "TCSCalorimeterHit.hh"
#include "TCSCalorimeterPMTHit.hh"
#include "TCSHodoscopeHit.hh"
#include "TCSTrackerHit.hh"
#include "TCSTargetHit.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4SDManager.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4SystemOfUnits.hh"    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSEventAction::TCSEventAction(TCSHistoManager *histo)
  : G4UserEventAction(), fEdep(0.), fNTar_epos(0), fNTar_eneg(0),
    fHistoManager(histo), fPrintModulo(0),
    fTargetCollID(-1), fCalorimeterCollID(-1), fCalorimeterPMTCollID(-1),
    //fHodoscopeCollID(-1), //temp jun14 
    //fTrackerCollID(-1), //temp jun14
    fEvtNo(-1)
{
  fPrintModulo = 1000;
  //  fPrintModulo = 1;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSEventAction::~TCSEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSEventAction::BeginOfEventAction(const G4Event* evt)
{
  fEvtNo = evt->GetEventID();
  
  G4SDManager * SDman = G4SDManager::GetSDMpointer();

  if(fCalorimeterCollID<0)
  {
    //    G4cout << "  Getting calorimeter collection id..." << G4endl;
    //G4String colNam="CalorimeterHitsCollection";
    //fCalorimeterCollID = SDman->GetCollectionID(colNam[0]);
    fCalorimeterCollID = SDman->GetCollectionID("CalorimeterHitsCollection");
    //G4cout << "  Calorimeter collection id = " << fCalorimeterCollID << G4endl;
  }

  if(fCalorimeterPMTCollID<0)
  {
    //    G4cout << "  Getting calorimeter PMT collection id..." << G4endl;
    fCalorimeterPMTCollID =
      SDman->GetCollectionID("CalorimeterPMTHitsCollection");
    //G4cout << "  Calorimeter PMT collection id = " << fCalorimeterPMTCollID << G4endl;
  }
  /*
  if(fHodoscopeCollID<0)
  {
    fHodoscopeCollID = SDman->GetCollectionID("HodoscopeHitsCollection");
  }

  if(fTrackerCollID<0)
    {
      fTrackerCollID = SDman->GetCollectionID("TrackerHitsCollection");
    }
  */ //temp jun 14
  if(fTargetCollID<0)
  {
    fTargetCollID = SDman->GetCollectionID("TargetHitsCollection");
    //    G4cout << "  Target collection id = " << fTargetCollID << G4endl;
  }

  // initialization of per event quantities

  fEdep = 0.;   //Per event total energy deposition in the target.
  fNTar_epos = 0;
  fNTar_eneg = 0;
    
  fHistoManager->Reset();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSEventAction::EndOfEventAction(const G4Event* event)
{   
  // Accumulate statistics in TCSRun.
  TCSRun* run = static_cast<TCSRun*>(
		       G4RunManager::GetRunManager()->GetNonConstCurrentRun() );
  run->AddEdep(fEdep);
  run->AddNTar_eneg(fNTar_eneg);
  run->AddNTar_epos(fNTar_epos);

  // Fill histogram with event's energy deposition.
  fHistoManager->FillHisto(1, fEdep/MeV);

  // get number of stored trajectories

  //G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
  //  G4int n_trajectories = 0;
  //  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

  //  G4cout << "TCSEventAction::EndOfEventAction: n_trajectories = " 
  //	 << n_trajectories << G4endl;
  //  getchar();

  //  G4int eventID = event->GetEventID();

  // Hit collection for this event.

  G4HCofThisEvent * HCE = event->GetHCofThisEvent();

  // Target hits.

  TCSTargetHitsCollection* TC = 0;
  if(HCE) {
    TC = (TCSTargetHitsCollection*)(HCE->GetHC(fTargetCollID));
    //    G4cout << "  Found target hit collection." << G4endl;
  }

  if(TC) {
    int n_hit = TC->entries();
    //    G4cout << "  target n_hit = " << n_hit << G4endl;
    for(int i=0;i<n_hit;i++) {
      G4int pid =(*TC)[i]->GetPID();
      G4double energy=(*TC)[i]->GetEnergy();
      fHistoManager->AddHit(energy/MeV, pid);
    }
  }

  //Check hit container's consistency first.
  if (!fHistoManager->CheckTargetHitCont()) {
    cout <<"*** TCSEventAction::EndOfEventAction: "
	 << "target hit container inconsistent! ***" << endl;
    //    getchar();
  }

  // Calorimeter hits.

  TCSCalorimeterHitsCollection* CC = 0;
  if(HCE) {
    CC = (TCSCalorimeterHitsCollection*)(HCE->GetHC(fCalorimeterCollID));
    //    G4cout << "  Found calorimeter hit collection." << G4endl;
  }

  if(CC) {
    int n_hit = CC->entries();
    //    G4cout << "  n_hit = " << n_hit << G4endl;

    for(int i=0;i<n_hit;i++) {
      G4int boundary_flag=(*CC)[i]->GetBoundaryFlag();
      //Fill Tree if track is within the calorimeter.
      if (boundary_flag == 0) {
	G4ThreeVector pos=(*CC)[i]->GetPos();
	//	G4int detpos = pos.getY() > 0. ? 1 : -1;
	G4int detpos = GetQuarter(pos.getX(), pos.getY());
	G4int col =(*CC)[i]->GetCol();
	G4int row =(*CC)[i]->GetRow();
	G4int pid =(*CC)[i]->GetPID();
	G4double energy=(*CC)[i]->GetEnergy();
	fHistoManager->AddHit(detpos, col, row, energy/MeV, pid);
      }
    }

    //Check hit container's consistency first.
    if (!fHistoManager->CheckCaloHitCont())
      cout <<"*** TCSEventAction::EndOfEventAction: "
	   << "calorimeter hit container inconsistent! ***" << endl;

    //    getchar();
  }

  //

  TCSCalorimeterPMTHitsCollection* CCPMT = 0;
  if(HCE) {
    CCPMT=(TCSCalorimeterPMTHitsCollection*)(HCE->GetHC(fCalorimeterPMTCollID));
    //    G4cout << "  Found calorimeter PMT hit collection." << G4endl;
  }

  if(CCPMT) {
    int n_hit = CCPMT->entries();
    //    G4cout << "  n_hit = " << n_hit << G4endl;
    for(int i=0;i<n_hit;i++) {
      G4ThreeVector pos=(*CCPMT)[i]->GetPos();
      G4int detpos = GetQuarter(pos.getX(), pos.getY());
      G4int col =(*CCPMT)[i]->GetCol();
      G4int row =(*CCPMT)[i]->GetRow();
      G4int pid =(*CCPMT)[i]->GetPID();
      int npe=(*CCPMT)[i]->GetNpe();
      fHistoManager->AddHit(detpos, col, row, npe, pid);
    }

    //Check hit container's consistency first.
    if (!fHistoManager->CheckCaloHitCont())
      cout <<"*** TCSEventAction::EndOfEventAction: "
	   << "calorimeter hit container inconsistent! ***" << endl;

    //    getchar();
  }
  /*
  // Hodoscope hits.

  TCSHodoscopeHitsCollection* HC = 0;
  if(HCE) {
    HC = (TCSHodoscopeHitsCollection*)(HCE->GetHC(fHodoscopeCollID));
    //    G4cout << "  Found hodoscope hit collection." << G4endl;
  }

  if(HC) {
    int n_hit = HC->entries();
    //    G4cout << "  n_hit = " << n_hit << G4endl;

    for(int i=0;i<n_hit;i++) {
      G4int boundary_flag=(*HC)[i]->GetBoundaryFlag();
      //Fill Tree if track is within the hodoscope.
      if (boundary_flag == 0) {
	G4ThreeVector pos=(*HC)[i]->GetPos();
	G4int detpos = GetQuarter(pos.getX(), pos.getY());
	G4int col =(*HC)[i]->GetCol();
	G4int row =(*HC)[i]->GetRow();
	G4int pid =(*HC)[i]->GetPID();
	G4double energy=(*HC)[i]->GetEnergy();
	fHistoManager->AddHodoHit(detpos, col, row, energy/MeV, pid);
      }
    }

    //Check hit container's consistency first.
    if (!fHistoManager->CheckHodoHitCont())
      cout <<"*** TCSEventAction::EndOfEventAction: "
	   << "hodoscope hit container inconsistent! ***" << endl;

    //    getchar();
  }

  // Tracker hits.

  TCSTrackerHitsCollection* TrC = 0;
  if(HCE) {
    TrC = (TCSTrackerHitsCollection*)(HCE->GetHC(fTrackerCollID));
    //    G4cout << "  Found tracker hit collection." << G4endl;

    if(TrC) {
      //      G4cout << "    Add Tracker hits:" << G4endl;
      AddTrackerHit(TrC, fHistoManager->fTrackerHitCont);
    }

  }
  */
  int nvertex =  event->GetNumberOfPrimaryVertex();

  //  cout << "TCSEventAction::EndOfEventAction: nvertex = " << nvertex << endl;

  for (int iv=0; iv<nvertex; iv++) {

    G4PrimaryVertex* vertex = event->GetPrimaryVertex(iv);
    //    vertex->Print();
    //    getchar();

    G4ThreeVector origin = vertex->GetPosition();

    int nparticles = vertex->GetNumberOfParticle();
    //    cout << "   number of paticles = " << nparticles << endl;

    //Alow externally generated multi-particle events.
    //    if (nparticles != 1 && nparticles != 3)
    //      cout << "TCSEventAction::EndOfEventAction: wrong nparticles = "
    //	   << nparticles << " ! *****" << endl;

    for (int ip=0; ip<nparticles; ip++) {

      G4PrimaryParticle* particle = vertex->GetPrimary(ip);
      double E = particle->GetTotalEnergy();
      G4ThreeVector P = particle->GetMomentum();
      int PID = particle->GetPDGcode();

      fHistoManager->SetTCSVertex(E, P, origin, PID);

      if (iv==0 && ip==0)
	fHistoManager->SetBeam(E, P, origin, PID);

    }

  }

  // if (CC || HC || TrC) { temp jun14
  if (CC) {
    fHistoManager->FillTrees();
    //    getchar();
  }
  

  //  if (fEvtNo%1000 == 0) fHistoManager->autosave();
  if (fEvtNo%fPrintModulo == 0)
    G4cout << "event: " << fEvtNo
	   << "  e- : " << run->GetNTar_eneg()
      	   << "  e+ : " << run->GetNTar_epos()
	   << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
void TCSEventAction::AddHodoHit(TCSHodoHitsCollection* HC,
                                HodoHitContainer& HodoHitCont)
{
      int n_hit = HC->entries();
      //      G4cout << "      HC n_hit = " << n_hit << G4endl;

      for(int i=0;i<n_hit;i++) {
        G4int boundary_flag=(*HC)[i]->GetBoundaryFlag();
	//     G4cout << "        boundary_flag = " << boundary_flag << G4endl;
        //Fill Tree if track is within the hodoscope.
        if (boundary_flag == 0) {
          G4ThreeVector pos=(*HC)[i]->GetPos();
	  //          G4int detpos = pos.getY() > 0. ? 1 : -1;
	  G4int detpos = GetQuarter(pos.getX(), pos.getY());
          G4int chan =(*HC)[i]->GetChannel();
	  G4int pid =(*HC)[i]->GetPID();
          G4double energy=(*HC)[i]->GetEnergy();
          fHistoManager->AddHit(detpos, chan, energy/MeV, pid, HodoHitCont);
        }
      }

      //Check hit container's consistency first.
      if (!fHistoManager->CheckHodoHitCont(HodoHitCont))
        cout <<"*** TCSEventAction::EndOfEventAction: "
             << "hodoscope hit container inconsistent! ***" << endl;
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSEventAction::AddTrackerHit(TCSTrackerHitsCollection* HC,
				   TrackerHitContainer& TrackerHitCont)
{
      int n_hit = HC->entries();
      //G4cout << "TCSEventAction::AddTrackerHit: HC n_hit = " << n_hit<<G4endl;

      for(int i=0;i<n_hit;i++) {
	G4double x     = (*HC)[i]->GetX();
	G4double y     = (*HC)[i]->GetY();
	G4double edep  = (*HC)[i]->GetEdep();
	G4double length  = (*HC)[i]->GetLength();
	G4int    det   = (*HC)[i]->GetQuarter();
	G4int    layer = (*HC)[i]->GetLayer();
	G4int    pid = (*HC)[i]->GetPID();
	G4int    pidorig = (*HC)[i]->GetPIDOrig();
	G4int    trackid = (*HC)[i]->GetTrackID();
	fHistoManager->AddHit(x, y, edep/keV, length/mm, det, layer,
			      pid, pidorig, trackid, TrackerHitCont);
//	G4cout << "x = " << x << "  y = " << y << G4endl;
      }

      //Check hit container's consistency first.
      if (!fHistoManager->CheckTrackerHitCont(TrackerHitCont))
        cout <<"*** TCSEventAction::EndOfEventAction: "
             << "tracker hit container inconsistent! ***" << endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////
/*
void TCSEventAction::AddTrackerHit(TCSTrackerHitsCollection* HC,
				   TrackerHitContainer& TrackerHitCont,
				   TrackerHitContainer& TrackerFluxCont)
{
      int n_hit = HC->entries();
      //      G4cout << "      HC n_hit = " << n_hit << G4endl;

      for(int i=0;i<n_hit;i++) {
        G4int boundary_flag=(*HC)[i]->GetBoundaryFlag();
	//	G4cout << "        boundary_flag = " << boundary_flag << G4endl;
	G4ThreeVector pos=(*HC)[i]->GetPos();
	G4int detpos = pos.getY() > 0. ? 1 : -1;
	G4int chan =(*HC)[i]->GetChannel();
	G4int pid =(*HC)[i]->GetPID();
	G4double energy=(*HC)[i]->GetEnergy();

        if (boundary_flag == 0)
          fHistoManager->AddHit(detpos, chan, energy/MeV, pid, TrackerHitCont);
	else
          fHistoManager->AddFlux(detpos,chan, energy/MeV, pid, TrackerFluxCont);

      }

      //Check hit container's consistency first.
      if (!fHistoManager->CheckTrackerHitCont(TrackerHitCont))
        cout <<"*** TCSEventAction::EndOfEventAction: "
             << "hodoscope hit container inconsistent! ***" << endl;
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
int TCSEventAction::GetQuarter(double x, double y) {
  // Quarter in X-Y plane: left top, right top, right bottom, left bottom.
  // x>0, y>0 : 0
  // x<0, y>0 : 1
  // x<0> y<0 : 2
  // x>0, y<0 : 3

  int q = -1;
  if (x>0.)
    if (y>0.)
      q=0;
    else
      q=3;
  else
    if (y>0.)
      q=1;
    else
      q=2;

  return q;
}
*/ // commented as not a very good method to find the copy number or q 

int TCSEventAction::GetQuarter(double x, double y) {
  // Quarter in X-Y plane: left top, right top, right bottom, left bottom.
  // x>0, y>0 : 0
  // x<0, y>0 : 1
  // x<0> y<0 : 2
  // x>0, y<0 : 3

  int q = 0;

  return q;
}
