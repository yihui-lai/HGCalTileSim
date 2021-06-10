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

  runformat->tile_layer = DetectorConstruction->GetTileN();
  runformat->is_ESR = DetectorConstruction->GetWrapESR();
  runformat->sipm_width = DetectorConstruction->GetSiPMX();
  runformat->sipm_rim   = DetectorConstruction->GetSiPMRim();
  runformat->sipm_stand = DetectorConstruction->GetSiPMStand();

  runformat->dimple_rad = DetectorConstruction->GetDimpleRadius();
  runformat->dimple_ind = DetectorConstruction->GetDimpleIndent();

  runformat->abs_mult = DetectorConstruction->GetTileAbsMult();
  runformat->wrap_ref = DetectorConstruction->GetWrapReflect();

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
  format->E_dep_tot = 0;
  format->E_dep_nonion = 0;
  format->genphotons = 1.0;
  format->TrackerBX = -99;
  format->TrackerBY = -99;
  format->TrackerBZ = -99;
  format->TrackerFX = -99;
  format->TrackerFY = -99;
  format->TrackerFZ = -99;
  for(int ip=0;ip<10000;ip++){format->gen_z[ip]=-99;format->gen_z2[ip]=-99;}
  format->index = 0;
}

void
LYSimAnalysis::EndOfEvent( const G4Event* event )
{
  if( generatorAction ){
     format->genphotons = generatorAction->NSources();
   } else {
      //directly cout photons in step action, maybe a better way to do
     //format->genphotons = 1.0;// TODO: Somehow recover the number of generatoed photons. 
   }

  format->nphotons     = GetNPhotons( event );
  format->savedphotons = std::min( format->genphotons
                                 , (unsigned)LYSIMFORMAT_MAX_PHOTONS );

  G4TrajectoryContainer* trajectory_list = event->GetTrajectoryContainer();
  G4Navigator* navigator                 =
    G4TransportationManager::GetTransportationManager()->GetNavigator( "World" );

  if( generatorAction ){
    assert( format->genphotons == trajectory_list->size() );
  }

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

// add for counting average reflectivity
    unsigned pcb_ref  = 0;

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




    unsigned sipm_touch  = 0;
    unsigned sipm_ref  = 0;    
    bool touch_sipm=false;
 
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
    format->NumPCBHitandRef[saveindex]  = pcb_ref;
    format->NumSiPMTouch[saveindex]  = sipm_touch;
    format->NumSiPMHitandRef[saveindex]  = sipm_ref;

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
//  tree->Write( NULL, TObject::kOverwrite );

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

  //for(int ip=0;ip<adc_sample;ip++) format->arr_time[ip]=0;
  for( int i = 0; i < hits->entries(); ++i ){
    assert( ( *hits )[i]->GetPhotonCount() == 1 );

    //int arrival_time = ( *hits )[i]->GetTime()/ns;
    //format->arr_time[arrival_time+5]++;

    ++EventPhotonCount;
  }

  return EventPhotonCount;
}

// Implementation of helper functions
static bool IsSiPMTrajectory( G4Navigator* nav, const G4ThreeVector& endpoint )
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

void LYSimAnalysis::addenergy(double tot, double nonion){
  format->E_dep_tot +=tot;
  format->E_dep_nonion +=nonion;
}
void LYSimAnalysis::addgenphoton(){
  format->genphotons++;
}
bool LYSimAnalysis::IsGenProton( ){
  if(generatorAction) return false; 
  else return true;
}
void LYSimAnalysis::AddTracker(int id, double xp, double yp, double zp){
  if(id==0){
    format->TrackerFX = xp;
    format->TrackerFY = yp;
    format->TrackerFZ = zp;
  }else if(id==1){
    format->TrackerBX = xp;
    format->TrackerBY = yp;
    format->TrackerBZ = zp;
  }
}
void LYSimAnalysis::addgenz(double z, double z2){
  format->gen_z[format->index] = z;
  format->gen_z2[format->index] = z2;
  if(format->index<10000-1) format->index++;
  //else cout<<"too many photons"<<endl;
}


