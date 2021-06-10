#ifdef CMSSW_GIT_HASH
#include "HGCalTileSim/Tile/interface/LYSimAnalysis.hh"
#include "HGCalTileSim/Tile/interface/LYSimDetectorConstruction.hh"
#include "HGCalTileSim/Tile/interface/LYSimFormat.hh"
#include "HGCalTileSim/Tile/interface/LYSimPrimaryGeneratorAction.hh"
#include "HGCalTileSim/Tile/interface/LYSimProtonGeneratorAction.hh"
#include "HGCalTileSim/Tile/interface/LYSimScintillation.hh"
#include "HGCalTileSim/Tile/interface/LYSimTrajectoryPoint.hh"
#else
#include "LYSimAnalysis.hh"
#include "LYSimDetectorConstruction.hh"
#include "LYSimFormat.hh"
#include "LYSimPrimaryGeneratorAction.hh"
#include "LYSimProtonGeneratorAction.hh"
#include "LYSimScintillation.hh"
#include "LYSimTrajectoryPoint.hh"
#endif

#include "g4root.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4SDManager.hh"
#include "G4TransportationManager.hh"
#include "G4UnitsTable.hh"

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

#include <iostream>
#include <signal.h>

using namespace std;
using namespace CLHEP;

// Static Helper fuction
static void
CloseAll( int signal )
{
  LYSimAnalysis::GetInstance()->EndOfExperiment();
  std::exit( signal );
}

static bool IsSiPMTrajectory( G4Navigator*, const G4ThreeVector& );

// LYSimAnalysis Programs
LYSimAnalysis* LYSimAnalysis::singleton = 0;

LYSimAnalysis::LYSimAnalysis() :
  generatorAction( nullptr ),
  protonAction( nullptr )
{
}

LYSimAnalysis::~LYSimAnalysis()
{
}

void
LYSimAnalysis::PrepareExperiment()
{
  file      = TFile::Open( filename.c_str(), "RECREATE" );
  tree      = new TTree( "LYSim", "LYSim" );
  runtree   = new TTree( "LYSimRun", "LYSimRun" );
  format    = new LYSimFormat();
  runformat = new LYSimRunFormat();
  format->AddToTree( tree );
  runformat->AddToTree( runtree );

  // Preparing the signal termination function
  signal( SIGABRT, CloseAll );
  // signal( SIGFPE,  CloseAll );
  // signal( SIGILL,  CloseAll );
  signal( SIGINT,  CloseAll );
  // signal( SIGSEGV, CloseAll );
  signal( SIGTERM, CloseAll );
}

void
LYSimAnalysis::PrepareNewRun( const G4Run* )
{
  runformat->tile_x = DetectorConstruction->GetTileX();
  runformat->tile_y = DetectorConstruction->GetTileY();
  runformat->tile_z = DetectorConstruction->GetTileZ();

  runformat->sipm_width = DetectorConstruction->GetSiPMX();
  runformat->sipm_rim   = DetectorConstruction->GetSiPMRim();
  runformat->sipm_stand = DetectorConstruction->GetSiPMStand();

  runformat->dimple_rad  = DetectorConstruction->GetDimpleRadius();
  runformat->dimple_ind  = DetectorConstruction->GetDimpleIndent();
  runformat->dimple_type = DetectorConstruction->GetDimpleType();

  runformat->abs_mult     = DetectorConstruction->GetTileAbsMult();
  runformat->wrap_ref     = DetectorConstruction->GetWrapReflect();
  runformat->tile_alpha   = DetectorConstruction->GetTileAlpha();
  runformat->dimple_alpha = DetectorConstruction->GetDimpleAlpha();

  runformat->pcb_rad = DetectorConstruction->GetPCBRadius();
  runformat->pcb_ref = DetectorConstruction->GetPCBReflect();

  if( generatorAction ){
    runformat->beam_x = generatorAction->GetBeamX();
    runformat->beam_y = generatorAction->GetBeamY();
    runformat->beam_w = generatorAction->GetWidth();
  } else if( protonAction ){
    runformat->beam_x = protonAction->GetBeamX();
    runformat->beam_y = protonAction->GetBeamY();
    runformat->beam_w = protonAction->GetWidth();
  } else {
    runformat->beam_x = 0;
    runformat->beam_y = 0;
    runformat->beam_w = 0;
  }

#ifdef CMSSW_GIT_HASH
  runformat->UpdateHash();
#endif
}

void
LYSimAnalysis::PrepareNewEvent( const G4Event* event )
{
  // All primary vertex's in the event share the same x,y values
  format->beam_x   = event->GetPrimaryVertex()->GetX0();
  format->beam_y   = event->GetPrimaryVertex()->GetY0();
  format->run_hash = runformat->run_hash;
}

void
LYSimAnalysis::EndOfEvent( const G4Event* event )
{
  if( generatorAction ){
    format->genphotons = generatorAction->NSources();
  } else {
    format->genphotons = 1.0;// TODO: Somehow recover the number of generatoed photons.
  }

  format->nphotons     = GetNPhotons( event );
  format->savedphotons = std::min( format->genphotons
                                 , (unsigned)LYSIMFORMAT_MAX_PHOTONS );

  G4TrajectoryContainer* trajectory_list = event->GetTrajectoryContainer();
  G4Navigator* navigator                 =
    G4TransportationManager::GetTransportationManager()
    ->GetNavigator( "World" );

  assert( format->genphotons == trajectory_list->size() );

  unsigned nhits     = format->nphotons;
  unsigned saveindex = 0;

  // Getting the number of saved photons.
  const unsigned target_hit_photons
    = std::max( unsigned(1)
              , format->savedphotons * format->nphotons / format->genphotons );

  unsigned num_hit_photons = 0;
  unsigned num_los_photons = 0;


  for( size_t i = 0; i < trajectory_list->size()
       && num_hit_photons + num_los_photons < format->savedphotons; ++i ){
    G4VTrajectory* trajectory = ( *trajectory_list )[i];

    unsigned wrapbounce = 0;
    unsigned pcbbounce  = 0;
    double tracklength  = 0;
    bool isdetected     = false;

    for( int j = 0; j < trajectory->GetPointEntries(); ++j ){
      const G4ThreeVector pos_end   = trajectory->GetPoint( j )->GetPosition();
      const G4ThreeVector pos_start = j == 0 ? pos_end :
                                      trajectory->GetPoint( j-1 )->GetPosition();

      G4VPhysicalVolume* volume
        = navigator->LocateGlobalPointAndSetup( pos_end );
      if( volume->GetName() == "Wrap" ){
        ++wrapbounce;
      } else if( volume->GetName() == "PCB" ){
        ++pcbbounce;
      }

      tracklength += ( pos_end - pos_start ).mag();
    }

    const G4ThreeVector endpoint
      = trajectory->GetPoint( trajectory->GetPointEntries()-1 )->GetPosition();

    if( IsSiPMTrajectory( navigator, endpoint ) && nhits ){
      --nhits;
      isdetected = true;
    }

    if( isdetected ){
      if( num_hit_photons < target_hit_photons ){
        num_hit_photons++;
      } else {
        continue;
      }
    } else {
      if( num_los_photons < format->savedphotons - target_hit_photons ){
        num_los_photons++;
      } else {
        continue;
      }
    }

    format->NumWrapReflection[saveindex] = wrapbounce;
    format->NumPCBReflection[saveindex]  = pcbbounce;
    format->IsDetected[saveindex]        = isdetected;
    format->OpticalLength[saveindex]
      = tracklength / LYSimFormat::opt_length_unit;
    format->EndX[saveindex]
      = endpoint.x() / LYSimFormat::end_pos_unit;
    format->EndY[saveindex]
      = endpoint.y() / LYSimFormat::end_pos_unit;
    ++saveindex;
    assert( saveindex == num_hit_photons + num_los_photons );
  }

#ifdef CMSSW_GIT_HASH
  format->UpdateHash();
#endif

  // Filling the tree
  tree->Fill();
  // tree->Write( NULL, TObject::kOverwrite );

#ifdef CMSSW_GIT_HASH// Disabling event saving for non-interactive stuff
  G4EventManager::GetEventManager()
  ->GetNonconstCurrentEvent()
  ->KeepTheEvent( false );
#endif
}


void
LYSimAnalysis::EndOfRun( const G4Run* )
{
  runtree->Fill();
  runtree->Write( NULL, TObject::kOverwrite );
}

void
LYSimAnalysis::EndOfExperiment()
{
  // Flushing run tree if no entries is present
  if( runtree->GetEntries() == 0 ){
    runtree->Fill();
  }

  runtree->Write( NULL, TObject::kOverwrite );
  tree->Write( NULL, TObject::kOverwrite );
  file->Close();
}

int
LYSimAnalysis::GetNPhotons( const G4Event* event )
{
  static G4SDManager* SDman = G4SDManager::GetSDMpointer();
  static G4int hitCollID    = SDman->GetCollectionID( "PMTHitsCollection" );

  G4HCofThisEvent* hitsCollections = event->GetHCofThisEvent();
  if( !hitsCollections ){
    G4cerr << "hitsCollection not found" << G4endl;
    return -1;
  }

  LYSimPMTHitsCollection* hits
    = static_cast<LYSimPMTHitsCollection*>(
        hitsCollections->GetHC( hitCollID ) );

  unsigned EventPhotonCount = 0;

  for( size_t i = 0; i < hits->entries(); ++i ){
    assert( ( *hits )[i]->GetPhotonCount() == 1 );
    ++EventPhotonCount;
  }

  return EventPhotonCount;
}

// Implementation of helper functions
static bool
IsSiPMTrajectory( G4Navigator* nav, const G4ThreeVector& endpoint )
{
  const G4ThreeVector zvec( 0, 0, 1 );
  const G4VPhysicalVolume* endvol
    = nav->LocateGlobalPointAndSetup( endpoint, &zvec, false, false );
  // Expecting SiPM hits to be pointing in the z direction.
  if( endvol ){
    return endvol->GetName() == "SiPM";
  } else {// Rare occassion that endvol is not found
    return false;
  }
}
