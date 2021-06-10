#include "HGCalTileSim/Tile/interface/LYSimAnalysis.hh"
#include "HGCalTileSim/Tile/interface/LYSimDetectorConstruction.hh"
#include "HGCalTileSim/Tile/interface/LYSimPhysicsList.hh"
#include "HGCalTileSim/Tile/interface/LYSimPrimaryGeneratorAction.hh"
#include "HGCalTileSim/Tile/interface/LYSimProtonGeneratorAction.hh"
#include "HGCalTileSim/Tile/interface/LYSimSteppingAction.hh"
#include "HGCalTileSim/Tile/interface/LYSimTrackingAction.hh"

#include "UserUtils/Common/interface/ArgumentExtender.hpp"

#include "G4RunManager.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"

int
main( int argc, char** argv )
{
  usr::po::options_description desc( "Running a run with a certain geometry" );
  desc.add_options()
    ( "beamx,x",       usr::po::defvalue<double>( 0 ),     "x center of beam [mm]" )
    ( "beamy,y",       usr::po::defvalue<double>( 0 ),     "y center of beam [mm]" )
    ( "tilez,z",       usr::po::defvalue<double>( 3 ),     "thick of tile [mm]" )
    ( "tilewidth,L",   usr::po::defvalue<double>( 30 ),    "Length of tile [mm]" )
    ( "beamwidth,w",   usr::po::defvalue<double>( 1.5 ),   "width of beam [mm]" )
    ( "dimplerad,r",   usr::po::defvalue<double>( 3.0 ),   "Dimple radius [mm]" )
    ( "dimpleind,d",   usr::po::defvalue<double>( 1.5 ),   "Dimple indent [mm]" )
    ( "dimpletype,T",  usr::po::defvalue<int>( 0 ),        "Dimple Type (0:Spherical, 1:Ellipsoid, 2:FlatDome, 3:Cylinder" )
    ( "sipmwidth,W",   usr::po::defvalue<double>( 2.0 ),   "SiPM active width [mm]" )
    ( "sipmrim,R",     usr::po::defvalue<double>( 0.1 ),   "SiPM Rim width [mm]" )
    ( "sipmstand,S",   usr::po::defvalue<double>( 0.3 ),   "SiPM stand [mm]" )
    ( "absmult,a",     usr::po::defvalue<double>( 1 ),     "Multple of inbuilt absorption length" )
    ( "wrapreflect,m", usr::po::defvalue<double>( 0.985 ), "Wrap reflectivity" )
    ( "tilealpha,A",   usr::po::defvalue<double>( 0.01 ),  "Tile surface alpha value" )
    ( "dimplealpha,D", usr::po::defvalue<double>( 0.1 ),   "Dimple surface alpha value" )
    ( "pcbreflect,p",  usr::po::defvalue<double>( 0.8 ),   "Average reflectivity of Exposed PCB back plane" )
    ( "pcbradius,b",   usr::po::defvalue<double>( 2.5 ),   "Radius of Exposed PCB back plane" )
    ( "NEvents,N",     usr::po::defvalue<unsigned>( 1 ),   "Number of events to run" )
    ( "useProton,U",    usr::po::defvalue<int>( 1 ),        "Flag to switch the source to a true proton source" )
    ( "output,o",      usr::po::defvalue<std::string>( "test.root" ), "output file" )
  ;

  usr::ArgumentExtender args;
  args.AddOptions( desc );
  args.ParseOptions( argc, argv );

  const double x_center  = args.Arg<double>( "beamx"       );
  const double y_center  = args.Arg<double>( "beamy"       );
  const double tilewidth = args.Arg<double>( "tilewidth"   );
  const double tilez = args.Arg<double>( "tilez"   );
  const double width     = args.Arg<double>( "beamwidth"   );
  const double dimplerad = args.Arg<double>( "dimplerad"   );
  const double dimpleind = args.Arg<double>( "dimpleind"   );
  const int dimpletype   = args.Arg<int>( "dimpletype"  );
  const double sipmwidth = args.Arg<double>( "sipmwidth"   );
  const double sipmrim   = args.Arg<double>( "sipmrim"     );
  const double sipmstand = args.Arg<double>( "sipmstand"   );
  const double absmult   = args.Arg<double>( "absmult"     );
  const double wrapref   = args.Arg<double>( "wrapreflect" );
  const double tilealpha = args.Arg<double>( "tilealpha"   );
  const double dimpalpha = args.Arg<double>( "dimplealpha" );
  const double pcbref    = args.Arg<double>( "pcbreflect"  );
  const double pcbrad    = args.Arg<double>( "pcbradius"   );
  const unsigned N       = args.Arg<unsigned>( "NEvents" );
  const int useProton   = args.Arg<int>( "useProton" );
  std::string filename   = args.Arg<std::string>( "output" );
  filename.insert( filename.length()-5, "_" + usr::RandomString( 6 ) );

  // Must initialize Run Manager first
  G4RunManager* runManager   = new G4RunManager;
  LYSimPhysicsList* physlist = new LYSimPhysicsList();
  // Overriding the detector parameters
  LYSimDetectorConstruction* detector = new LYSimDetectorConstruction();
  detector->SetDimpleRadius( dimplerad );
  detector->SetDimpleIndent( dimpleind );
  detector->SetDimpleType( dimpletype );
  detector->SetTileX( tilewidth );
  detector->SetTileY( tilewidth );
  detector->SetTileZ( tilez );
  detector->SetTileAbsMult( absmult );
  detector->SetWrapReflect( wrapref );
  detector->SetTileAlpha( tilealpha );
  detector->SetDimpleAlpha( dimpalpha );
  detector->SetSiPMX( sipmwidth );
  detector->SetSiPMY( sipmwidth );
  detector->SetSiPMRim( sipmrim );
  detector->SetSiPMStand( sipmstand );
  detector->SetPCBReflect( pcbref );
  detector->SetPCBRadius( pcbrad );

  runManager->SetUserInitialization( detector );
  runManager->SetUserInitialization( physlist );

  // Overriding the generator parameters
  if( useProton==0 ){
    LYSimPrimaryGeneratorAction* genaction
      = new LYSimPrimaryGeneratorAction( detector );
    genaction->SetBeamX( x_center );
    genaction->SetBeamY( y_center );
    genaction->SetWidth( width );
    LYSimAnalysis::GetInstance()->SetGeneratorAction( genaction );
    runManager->SetUserAction( genaction );
  } else {
    LYSimProtonGeneratorAction* genaction
      = new LYSimProtonGeneratorAction();
    genaction->SetBeamX( x_center );
    genaction->SetBeamY( y_center );
    genaction->SetWidth( width );
    LYSimAnalysis::GetInstance()->SetProtonGeneratorAction( genaction );
    runManager->SetUserAction( genaction );
  }

  // Construct LYSimAnalysis class
  LYSimAnalysis::GetInstance()->SetDetector( detector );
  LYSimAnalysis::GetInstance()->SetOutputFile( filename );


  runManager->SetUserAction( new LYSimAnalysis::RunAction() );
  runManager->SetUserAction( new LYSimAnalysis::EventAction() );
  runManager->SetUserAction( new LYSimTrackingAction() );
  runManager->SetUserAction(
    new LYSimSteppingAction( LYSimAnalysis::GetInstance() ) );

  // Preparing the LYSimAnalysis
  LYSimAnalysis::GetInstance()->PrepareExperiment();

  // Initialize G4 kernel
  runManager->Initialize();
  G4VisManager* visManager = new G4VisExecutive( "Quiet" );
  visManager->Initialize();
  G4UImanager* UIManager = G4UImanager::GetUIpointer();

  char cmd[1024];
  UIManager->ApplyCommand( "/control/verbose 1" );
  UIManager->ApplyCommand( "/run/verbose 0" );
  UIManager->ApplyCommand( "/process/setVerbose 0" );
  UIManager->ApplyCommand( "/tracking/verbose 0" );

  srand( time( NULL ) );

  sprintf( cmd, "/random/setSeeds %d %d", rand(), rand() );
  UIManager->ApplyCommand( cmd );

  std::cout << N << std::endl;
  runManager->BeamOn( N );

  LYSimAnalysis::GetInstance()->EndOfExperiment();

  // Job termination Free the store: user actions, physics_list and
  // detector_description are owned and deleted by the run manager, so they
  // should not be deleted in the main() program !
  delete visManager;
  delete runManager;

  return 0;
}
