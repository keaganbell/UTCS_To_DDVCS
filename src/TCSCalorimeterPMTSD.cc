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

#include "TCSCalorimeterPMTSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "TCSTrackInformation.hh"
#include "TCSCalorimeterConstruction.hh"

#include "G4OpBoundaryProcess.hh"

#include "G4ProcessManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSCalorimeterPMTSD::TCSCalorimeterPMTSD(const G4String& name,
					 const G4String& hitsCollectionName) 
 : G4VSensitiveDetector(name),
   fHitsCollection(NULL)
{
  collectionName.insert(hitsCollectionName);

  G4cout << "TCSCalorimeterPMTSD::TCSCalorimeterPMTSD: constructed" << G4endl;
  G4cout << "   name = " << name << G4endl;
  G4cout << "   hitsCollectionName = " << hitsCollectionName << G4endl;
  //  getchar();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCSCalorimeterPMTSD::~TCSCalorimeterPMTSD() 
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSCalorimeterPMTSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  //  G4cout << "TCSCalorimeterPMTSD::Initialize: creating hit collection. "
  //  << G4endl;
  //  getchar();

  fHitsCollection 
  = new TCSCalorimeterPMTHitsCollection(SensitiveDetectorName,
					collectionName[0]); 

  // Add this collection in hce

  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 

  //  G4cout << "TCSCalorimeterPMTSD::Initialize: initialized" << G4endl;
  //  getchar();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TCSCalorimeterPMTSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  //  G4cout << "TCSCalorimeterPMTSD::ProcessHits: entered" << G4endl;
  
  G4ParticleDefinition* particleType = step->GetTrack()->GetDefinition();
  if(particleType==G4OpticalPhoton::OpticalPhotonDefinition()) {
    //Optical photon only
  
    //Check to see if the partcile was actually at a boundary
    //Otherwise the boundary status may not be valid
    //Prior to Geant4.6.0-p1 this would not have been enough to check

    if(step->GetPostStepPoint()->GetStepStatus()==fGeomBoundary){

      G4VPhysicalVolume* postvol = 
	step->GetPostStepPoint()->GetPhysicalVolume();
      G4VPhysicalVolume* prevol = 
	step->GetPreStepPoint()->GetPhysicalVolume();

      G4OpBoundaryProcessStatus boundaryStatus=Undefined;
      static G4OpBoundaryProcess* boundary=NULL;
  
      //find the boundary process only once
      if(!boundary){
	G4ProcessManager* pm 
	  = step->GetTrack()->GetDefinition()->GetProcessManager();
	G4int nprocesses = pm->GetProcessListLength();
	G4ProcessVector* pv = pm->GetProcessList();
	G4int i;
	for( i=0;i<nprocesses;i++){
	  if((*pv)[i]->GetProcessName()=="OpBoundary"){
	    boundary = (G4OpBoundaryProcess*)(*pv)[i];
	    break;
	  }
	}
      }

      boundaryStatus=boundary->GetStatus();

      if (prevol->GetName() == "PMTWindow_phys" &&
	  postvol->GetName() == "Cathode_phys" &&
	  boundaryStatus == Detection) {

	G4StepPoint* postStepPoint = step->GetPostStepPoint();
	G4TouchableHandle theTouchable = postStepPoint->GetTouchableHandle();

	int nmod = 
	  step->GetPostStepPoint()->GetTouchableHandle()->GetCopyNumber(1);

	TCSCalorimeterConstruction cc;
	int ncol = cc.GetNCOL();
	int nrow = cc.GetNROW();
	int row = nmod/ncol;
	int col = nmod - row*ncol;

	if (row < 0 || row >= nrow)
	  G4cout << "*** TCSCalorimeterPMTSD::ProcessHits: wrong row = "
		 << row << G4endl;
	if (col < 0 || col >= ncol)
	  G4cout << "*** TCSCalorimeterPMTSD::ProcessHits: wrong col = "
		 << col << G4endl;

	TCSTrackInformation* info =
	  (TCSTrackInformation*)(step->GetTrack()->GetUserInformation());
	G4int pid = info->GetOriginalParticle()->GetPDGEncoding();

	G4ThreeVector pos = step->GetTrack()->GetPosition();

	TCSCalorimeterPMTHit* hit =
	  new TCSCalorimeterPMTHit(col, row, pid, 1, pos);
	fHitsCollection->insert( hit );

      }   //detection

    }     //boundary

  }       //optical photon

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TCSCalorimeterPMTSD::EndOfEvent(G4HCofThisEvent*)
{
  //  G4cout << "TCSCalorimeterPMTSD::EndOfEvent:" << G4endl;
  //  G4cout << "  Number of entries in hit collection = "
  //	 << fHitsCollection->entries() << G4endl;

  ////  for ( G4int i=0; i<fHitsCollection->entries(); i++ )
  ////    (*fHitsCollection)[i]->delete();

  //  getchar();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
