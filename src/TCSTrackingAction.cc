#include "TCSTrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4TrackVector.hh"
#include "TCSTrackInformation.hh"

TCSTrackingAction::TCSTrackingAction()
{;}

TCSTrackingAction::~TCSTrackingAction()
{;}

void TCSTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  //  G4cout << "TCSTrackingAction::PreUserTrackingAction:" << G4endl;
  if(aTrack->GetParentID()==0 && aTrack->GetUserInformation()==0)
  {
    //    G4cout << "   ParentID = 0, No Info " << G4endl;
    TCSTrackInformation* anInfo = new TCSTrackInformation(aTrack);
    G4Track* theTrack = (G4Track*)aTrack;
    theTrack->SetUserInformation(anInfo);
  }
  //  getchar();
}

void TCSTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
  if(secondaries)
  {
    TCSTrackInformation* info = (TCSTrackInformation*)(aTrack->GetUserInformation());
    size_t nSeco = secondaries->size();
    if(nSeco>0)
    {
      for(size_t i=0;i<nSeco;i++)
      { 
        TCSTrackInformation* infoNew = new TCSTrackInformation(info);
        (*secondaries)[i]->SetUserInformation(infoNew);
      }
    }
  }
}
