#include "HGCalTileSim/Tile/interface/LYSimAnalysis.hh"
#include "HGCalTileSim/Tile/interface/LYSimDetectorConstruction.hh"
#include "HGCalTileSim/Tile/interface/LYSimPhysicsList.hh"
#include "HGCalTileSim/Tile/interface/LYSimPrimaryGeneratorAction.hh"
#include "HGCalTileSim/Tile/interface/LYSimSteppingAction.hh"
#include "HGCalTileSim/Tile/interface/LYSimTrackingAction.hh"
#include "HGCalTileSim/Tile/interface/LYSimProtonGeneratorAction.hh"
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
    ( "beamz,z",       usr::po::defvalue<double>( 0 ),     "z center of beam [mm]" )
    ( "tileZ,l",   usr::po::defvalue<double>( 30 ),    "Length of tile [mm]" )
    ( "beamwidth,w",   usr::po::defvalue<double>( 1.5 ),   "width of beam [mm]" )
    ( "fiberZ,f",   usr::po::defvalue<float>( 5.2 ),   "fiber length [m]" )
    ( "fiberZshift,s",   usr::po::defvalue<double>( 1.7 ),   "fiber shift [m]" )
    ( "absmult,a",     usr::po::defvalue<double>( 1000 ),     "absorption length at 425nm, unit mm" )
    ( "yield,y",     usr::po::defvalue<double>( 10 ),     "light yield / keV" )
    ( "wrapreflect,m", usr::po::defvalue<float>( 0.985 ), "Wrap reflectivity" )
    ( "NEvents,N",     usr::po::defvalue<unsigned>( 1 ),   "Number of events to run" )
    ( "useProton,P",     usr::po::defvalue<int>( 1 ),  "Flag to switch the source to a true proton source" )
    ( "handwrap,H",      usr::po::defvalue<int>( 0 ),    "Flag to switch to handwrap" )
    ( "output,o",      usr::po::defvalue<std::string>( "test.root" ), "output file" )
  ;

  usr::ArgumentExtender args;
  args.AddOptions( desc );
  args.ParseOptions( argc, argv );

  const double x_center  = args.Arg<double>( "beamx"       );
  const double z_center  = args.Arg<double>( "beamz"       );
  const double tileZ = args.Arg<double>( "tileZ"   );
  const double width     = args.Arg<double>( "beamwidth"   );
  const double fiberZ     = args.Arg<float>( "fiberZ"   );
  const double fiberZshift     = args.Arg<double>( "fiberZshift"   );
  const double absmult   = args.Arg<double>( "absmult"     );
  const double yield   = args.Arg<double>( "yield"     );
  const double wrapref   = args.Arg<float>( "wrapreflect" );
  //const double tilealpha = args.Arg<double>( "tilealpha"   );
  //const double dimpalpha = args.Arg<double>( "dimplealpha" );
  //const double pcbref    = args.Arg<double>( "pcbreflect"  );
  //const double pcbrad    = args.Arg<double>( "pcbradius"   );
  const unsigned N       = args.Arg<unsigned>( "NEvents" );
  const bool useProton   = args.Arg<int>( "useProton" );
  const bool handwrap   = args.Arg<int>( "handwrap" );
  std::string filename   = args.Arg<std::string>( "output" );
  filename.insert( filename.length()-5, "_" + usr::RandomString( 6 ) );

  // Must initialize Run Manager first
  G4RunManager* runManager   = new G4RunManager;
  LYSimPhysicsList* physlist = new LYSimPhysicsList();
  // Overriding the detector parameters
  LYSimDetectorConstruction* detector = new LYSimDetectorConstruction();
  //detector->SetDimpleRadius( dimplerad );
  //detector->SetDimpleIndent( dimpleind );
  //detector->SetDimpleType( dimpletype );
  detector->SetTileZ( tileZ );
  detector->SetFiberZ( fiberZ );
  detector->SetFiberZoff( fiberZshift );
  detector->SetTileAbsMult( absmult );
  detector->SetTileScintillation(yield);
  detector->SetWrapReflect( wrapref );
  //detector->SetTileAlpha( tilealpha );
  //detector->SetDimpleAlpha( dimpalpha );
  //detector->SetSiPMX( sipmwidth );
  //detector->SetSiPMY( sipmwidth );
  //detector->SetSiPMRim( sipmrim );
  //detector->SetSiPMStand( sipmstand );
  //detector->SetPCBReflect( pcbref );
  //detector->SetPCBRadius( pcbrad );
  detector->Set_handwrap( handwrap );

  runManager->SetUserInitialization( detector );
  runManager->SetUserInitialization( physlist );

  // Overriding the generator parameters
  if( !useProton ){
    LYSimPrimaryGeneratorAction* genaction
      = new LYSimPrimaryGeneratorAction( detector );
    genaction->SetBeamX( x_center );
    genaction->SetBeamY( z_center );
    genaction->SetWidth( width );
    LYSimAnalysis::GetInstance()->SetGeneratorAction( genaction );
    runManager->SetUserAction( genaction );
  } else {
    LYSimProtonGeneratorAction* genaction
      = new LYSimProtonGeneratorAction();
    genaction->SetBeamX( x_center );
    genaction->SetBeamZ( z_center );
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
