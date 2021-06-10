#ifdef CMSSW_GIT_HASH
#include "HGCalTileSim/Tile/interface/LYSimAnalysis.hh"
#include "HGCalTileSim/Tile/interface/LYSimDetectorConstruction.hh"
#include "HGCalTileSim/Tile/interface/LYSimPMTSD.hh"
#include "HGCalTileSim/Tile/interface/LYSimSteppingAction.hh"
#include "HGCalTileSim/Tile/interface/LYSimSteppingMessenger.hh"
#else
#include "LYSimAnalysis.hh"
#include "LYSimDetectorConstruction.hh"
#include "LYSimPMTSD.hh"
#include "LYSimSteppingAction.hh"
#include "LYSimSteppingMessenger.hh"
#endif

#include "G4EventManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4UnitsTable.hh"

using namespace CLHEP;

LYSimSteppingAction::LYSimSteppingAction( LYSimAnalysis* ana )
  : G4UserSteppingAction(),
  maxtracklength( 500000.*mm ),
  messenger( new LYSimSteppingMessenger( this ) ),
  analysis( ana )
{
}

LYSimSteppingAction::~LYSimSteppingAction()
{
  delete messenger;
}

void
LYSimSteppingAction::UserSteppingAction( const G4Step* step )
{

  G4StepPoint *thePrePoint = step->GetPreStepPoint();
  G4StepPoint *thePostPoint = step->GetPostStepPoint();
  G4VPhysicalVolume *thePrePV = thePrePoint->GetPhysicalVolume();
  G4VPhysicalVolume *thePostPV = thePostPoint->GetPhysicalVolume();
  G4String thePrePVName = "";
  if (thePrePV) thePrePVName = thePrePV->GetName();
  G4String thePostPVName = "";
  if (thePostPV) thePostPVName = thePostPV->GetName();
  if (thePrePVName.contains("Tile"))
  LYSimAnalysis::GetInstance()->addenergy(step->GetTotalEnergyDeposit()/ GeV, step->GetNonIonizingEnergyDeposit()/ GeV);

  const G4ThreeVector &thePrePosition = thePrePoint->GetPosition();
  const G4ThreeVector &thePostPosition = thePostPoint->GetPosition();
  G4double global_x = thePrePosition.x() / mm;
  G4double global_y = thePrePosition.y() / mm;
  G4double global_z = thePrePosition.z() / mm;

  G4double global_x1 = thePostPosition.x() / mm;
  G4double global_y1 = thePostPosition.y() / mm;
  G4double global_z1 = thePostPosition.z() / mm;

  G4Track *theTrack = step->GetTrack();
  G4ParticleDefinition *particleType = theTrack->GetDefinition();
  // find the photon: 1) generated in Tile; 2) avoid double counting
  if (particleType == G4OpticalPhoton::OpticalPhotonDefinition()){
     if(thePrePVName.contains("Tile") && theTrack->GetCurrentStepNumber()==1 ){
         LYSimAnalysis::GetInstance()->addgenphoton();
         LYSimAnalysis::GetInstance()->addgenz(global_z, global_z1);
         //std::cout<<"get one"<<std::endl;
    }
  }
  //find the other process besides scintillation 
  if (particleType == G4OpticalPhoton::OpticalPhotonDefinition()){
     if(LYSimAnalysis::GetInstance()->IsGenProton() && theTrack->GetCreatorProcess()->GetProcessName()!="Scintillation")  std::cout<<theTrack->GetCreatorProcess()->GetProcessName()<<std::endl;
  }


  //try to catch the track of proton
  if( LYSimAnalysis::GetInstance()->IsGenProton() ){
    G4int TrPDGid = theTrack->GetDefinition()->GetPDGEncoding();
    if(thePrePVName.contains("TrackerFront") && TrPDGid==2212){
      LYSimAnalysis::GetInstance()->AddTracker(0,global_x,global_y,global_z);
    }
    if(thePrePVName.contains("TrackerBack") && TrPDGid==2212){
      LYSimAnalysis::GetInstance()->AddTracker(1,global_x,global_y,global_z);
    }
  }


  G4OpBoundaryProcessStatus boundaryStatus = Undefined;
  static G4OpBoundaryProcess* boundary     = NULL;

  // find the boundary process only once
  if( !boundary ){
    G4ProcessManager* pm = step->GetTrack()->GetDefinition()->GetProcessManager();
    G4int nprocesses     = pm->GetProcessListLength();
    G4ProcessVector* pv  = pm->GetProcessList();
    G4int i;

    for( i = 0; i < nprocesses; i++ ){
      if( ( *pv )[i]->GetProcessName() == "OpBoundary" ){
        boundary = (G4OpBoundaryProcess*)( *pv )[i];
        break;
      }
    }
  }

  // Stop if boundary is not found
   if( boundary ){
     boundaryStatus = boundary->GetStatus();

  switch( boundaryStatus ){
  case Detection:
    // Note, this assumes that the volume causing detection is the photocathode
    // because it is the only one with non-zero efficiency. Trigger sensitive
    // detector manually since photon is absorbed but status was Detection
  {
    //some absorbed photons also appear here, use this to avoid them
    if(thePostPVName == ("SiPM")){
      G4EventManager::GetEventManager()->KeepTheCurrentEvent();
      G4SDManager* SDman = G4SDManager::GetSDMpointer();
      LYSimPMTSD* pmtSD
        = (LYSimPMTSD*)SDman->FindSensitiveDetector( "/LYSimPMT" );
      if( pmtSD ){
        pmtSD->ProcessHits_constStep( step, NULL );
      }
    }
    break;
  }
  default:
    break;
  }

  G4double tracklength = step->GetTrack()->GetTrackLength();
  if( tracklength > maxtracklength ){
    G4cout << "Track length exceeded limit of " << maxtracklength/mm
           << "mm" << G4endl;
    step->GetTrack()->SetTrackStatus( fStopAndKill );
  }


}
}


