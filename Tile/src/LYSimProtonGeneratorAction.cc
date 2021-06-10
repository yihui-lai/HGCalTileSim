#ifdef CMSSW_GIT_HASH
#include "HGCalTileSim/Tile/interface/LYSimAnalysis.hh"
#include "HGCalTileSim/Tile/interface/LYSimProtonGeneratorAction.hh"
#include "HGCalTileSim/Tile/interface/ProjectPath.hh"
#else
#include "LYSimAnalysis.hh"
#include "LYSimProtonGeneratorAction.hh"
#include "ProjectPath.hh"
#endif

#include "G4GeneralParticleSource.hh"
#include "G4ParticleTypes.hh"
#include "G4SystemOfUnits.hh" 
#include "Randomize.hh"

LYSimProtonGeneratorAction::LYSimProtonGeneratorAction() :
   particleSource( new G4GeneralParticleSource() ),
   messenger( new LYSimProtonGeneratorMessenger( this ) )
 {
   std::cout << "[LYSIM] entering LYSIMProtonGeneratorAction" << std::endl;

particleSource->SetParticleDefinition(
     G4Proton::ProtonDefinition() );
   particleSource->GetCurrentSource()->SetNumberOfParticles( 1 );


   G4SPSAngDistribution* ang = particleSource->GetCurrentSource()->GetAngDist();
   ang->SetAngDistType( "iso" );
   ang->SetMinTheta( 0 );
   ang->SetMaxTheta( 0 );
   ang->SetMinPhi( 0 );
   ang->SetMaxPhi( 0 );
  G4SPSEneDistribution* ene = particleSource->GetCurrentSource()->GetEneDist();
   ene->SetEnergyDisType( "Mono" );
   ene->SetMonoEnergy( 120 * GeV );
 }

LYSimProtonGeneratorAction::~LYSimProtonGeneratorAction()
 {
   delete particleSource;
 }

 void
 LYSimProtonGeneratorAction::GeneratePrimaries( G4Event* anEvent )
 {
 G4SPSPosDistribution* pos = particleSource->GetCurrentSource()->GetPosDist();
   pos->SetPosDisType( "Plane" );
   pos->SetPosDisShape( "Square" );
   //pos->SetCentreCoords( G4ThreeVector(_beamx, _beamy, -95.0 *mm) );
   //pos->SetHalfX( 0.5 * _width );
   //pos->SetHalfY( 0.5 * _width );
   pos->SetCentreCoords( G4ThreeVector( ( 2*G4UniformRand()-1 )*_width + _beamx, ( 2*G4UniformRand()-1 )*_width + _beamy, -95.0 *mm) );

   particleSource->GeneratePrimaryVertex( anEvent );
 }

#include "G4UIcmdWithADouble.hh"
 #include "G4UIcmdWithADoubleAndUnit.hh"


 LYSimProtonGeneratorMessenger::LYSimProtonGeneratorMessenger( LYSimProtonGeneratorAction* act )
   : genaction( act )
 {
   detDir = new G4UIdirectory( "/LYSim/" );
   detDir->SetGuidance( " Proton Particle source " );

   SetBeamXCmd = new G4UIcmdWithADoubleAndUnit( "/LYSim/SetBeamX", this );
   SetBeamXCmd->SetGuidance( "Set central X position of the beam" );
   SetBeamXCmd->SetParameterName( "BeamX", false );
   SetBeamXCmd->SetUnitCategory( "Length" );
   SetBeamXCmd->SetDefaultUnit( "mm" );
   SetBeamXCmd->AvailableForStates( G4State_PreInit, G4State_Idle );

   SetBeamYCmd = new G4UIcmdWithADoubleAndUnit( "/LYSim/SetBeamY", this );
   SetBeamYCmd->SetGuidance( "Set central Y position of the beam" );
   SetBeamYCmd->SetParameterName( "BeamY", false );
   SetBeamYCmd->SetUnitCategory( "Length" );
   SetBeamYCmd->SetDefaultUnit( "mm" );
   SetBeamYCmd->AvailableForStates( G4State_PreInit, G4State_Idle );

   SetWidthCmd = new G4UIcmdWithADoubleAndUnit( "/LYSim/SetBeamWidth", this );
   SetWidthCmd->SetGuidance( "Set Width of Beam" );
   SetWidthCmd->SetParameterName( "BeamWidth", false );
   SetWidthCmd->SetUnitCategory( "Length" );
   SetWidthCmd->SetDefaultUnit( "mm" );
   SetWidthCmd->SetRange( "BeamWidth>=0" );
   SetWidthCmd->AvailableForStates( G4State_PreInit, G4State_Idle );
 }

 LYSimProtonGeneratorMessenger::~LYSimProtonGeneratorMessenger()
 {
   delete detDir;
   delete SetBeamXCmd;
   delete SetBeamYCmd;
   delete SetWidthCmd;
 }

 void
 LYSimProtonGeneratorMessenger::SetNewValue( G4UIcommand* command, G4String val )
 {
   if( command == SetBeamXCmd ){
     genaction->SetBeamX( G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( val ) );
   } else if( command == SetBeamYCmd ){
     genaction->SetBeamY( G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( val ) );
   } else if( command == SetWidthCmd ){
     genaction->SetWidth( G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( val ) );
   }
 }




