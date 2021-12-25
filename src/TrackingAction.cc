#include "TrackingAction.hh"
#include "G4Track.hh"

void TrackingAction::PostUserTrackingAction(const G4Track* track) {

  G4int    trackID  = track->GetTrackID(); 
//G4double tracklen = track->GetTrackLength();
 
 if (trackID == 1) { // Primary particle

   G4int status = 0; // absorbed, transmitted, reflected
  
   if (!track->GetNextVolume()) {
     
     if (track->GetMomentumDirection().x() > 0.) status = 1;
     else                                        status = 2;
   }

   if (0) G4cout << " status: " << status << G4endl;
 }
 else { // Secondary particle
   
   //if (track->GetDefinition()->GetPDGCharge() != 0.) 
 }
}
