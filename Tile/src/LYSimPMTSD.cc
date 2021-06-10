#ifdef CMSSW_GIT_HASH
#include "HGCalTileSim/Tile/interface/LYSimAnalysis.hh"
#include "HGCalTileSim/Tile/interface/LYSimPMTHit.hh"
#include "HGCalTileSim/Tile/interface/LYSimPMTSD.hh"
#else
#include "LYSimAnalysis.hh"
#include "LYSimPMTHit.hh"
#include "LYSimPMTSD.hh"
#endif

#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"

LYSimPMTSD::LYSimPMTSD( G4String name )
  : G4VSensitiveDetector( name ),
  fPMTHitsCollection( 0 )
{
  collectionName.insert( "PMTHitsCollection" );
}

LYSimPMTSD::~LYSimPMTSD(){}

void
LYSimPMTSD::Initialize( G4HCofThisEvent* )
{
  fPMTHitsCollection = new LYSimPMTHitsCollection( GetName(), collectionName[0] );
}

G4bool
LYSimPMTSD::ProcessHits( G4Step* aStep, G4TouchableHistory* )
{
  return false;// ProcessHits should not be called
}

G4bool
LYSimPMTSD::ProcessHits_constStep( const G4Step*       aStep,
                                   G4TouchableHistory* ROhist )
{
  G4Track* theTrack = aStep->GetTrack();

  // need to know if this is an optical photon
  if( theTrack->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition() ){ return false;}

  LYSimPMTHit* hit = new LYSimPMTHit();

G4StepPoint *thePostPointsd = aStep->GetPostStepPoint();
G4VPhysicalVolume *thePostPVsd = thePostPointsd->GetPhysicalVolume();
G4String thePostPVNamesd = thePostPVsd->GetName();
if(!(thePostPVNamesd=="SiPM")){ std::cout<<"Weird SD! "<<thePostPVNamesd<<std::endl;}

  hit->AddEnergy( theTrack->GetTotalEnergy() );
  hit->IncPhotonCount();// increment hit for the selected pmt
  hit->SetTime( aStep->GetPostStepPoint()->GetGlobalTime() );
  hit->SetLength( theTrack->GetTrackLength() );
  hit->SetBounceCount( theTrack->GetCurrentStepNumber() );

  fPMTHitsCollection->insert( hit );

  return true;

}

void
LYSimPMTSD::EndOfEvent( G4HCofThisEvent* HCE )
{
  // Store collection with event and keep ID
  static G4int HCID = -1;
  if( HCID < 0 ){ HCID = GetCollectionID( 0 );}
  HCE->AddHitsCollection( HCID, fPMTHitsCollection );

  G4int nHits = fPMTHitsCollection->entries();
  if( verboseLevel >= 1 ){
    G4cout << "[LYSim] PMT collection: " << nHits << " hits" << G4endl;
    if( verboseLevel >= 2 ){
      fPMTHitsCollection->PrintAllHits();
    }
  }
  fPMTHitsCollection->PrintAllHits();
}

void
LYSimPMTSD::clear(){}

void
LYSimPMTSD::DrawAll(){}

void
LYSimPMTSD::PrintAll(){}
