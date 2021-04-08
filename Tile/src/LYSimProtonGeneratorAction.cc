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
#include "TMath.h"

LYSimProtonGeneratorAction::LYSimProtonGeneratorAction() :
  particleSource( new G4GeneralParticleSource() ),
  messenger( new LYSimProtonGeneratorMessenger( this ) )
{
  std::cout << "[LYSIM] entering LYSIMProtonGeneratorAction" << std::endl;

  /**
   * This particle source will be directly controlled via the /gps commands In
   * user scripts and the ApplyCommand calls in the main function. We are going
   * to set up the common stuff here. For setting the via the interactive
   * session, look at:
   * http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/ForApplicationDeveloper/html/GettingStarted/generalParticleSource.html
   */
  //particleSource->SetParticleDefinition( G4Proton::ProtonDefinition() );
  particleSource->SetParticleDefinition( G4MuonPlus::MuonPlusDefinition() );
  particleSource->GetCurrentSource()->SetNumberOfParticles( 1 );


  G4SPSAngDistribution* ang = particleSource->GetCurrentSource()->GetAngDist();
  ang->SetAngDistType( "iso" );
  ang->SetMinTheta( 0 );
  ang->SetMaxTheta( 0 );
  //ang->SetMinPhi( TMath::Pi()/2 );
  //ang->SetMaxPhi( TMath::Pi()/2 );

  // Energy distribution.
  G4SPSEneDistribution* ene = particleSource->GetCurrentSource()->GetEneDist();
  ene->SetEnergyDisType( "Mono" );
  ene->SetMonoEnergy( 10 * GeV );
}

LYSimProtonGeneratorAction::~LYSimProtonGeneratorAction()
{
  delete particleSource;
}

void
LYSimProtonGeneratorAction::GeneratePrimaries( G4Event* anEvent )
{
  // We still need to generate the primary vertex.
  // Position distribution
  G4SPSPosDistribution* pos = particleSource->GetCurrentSource()->GetPosDist();
  pos->SetPosDisType( "Plane" );
  pos->SetPosDisShape( "Square" );
  pos->SetHalfX( 0.5 * _width );
  pos->SetHalfY( 0.5 * _width );
  pos->SetPosRot2(G4ThreeVector( 0., 0, -1 ));
  pos->SetCentreCoords( G4ThreeVector(_beamx, 90*mm, _beamz) );

  particleSource->GeneratePrimaryVertex( anEvent );
}

// ------------------------------------------------------------------------------

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

  SetBeamZCmd = new G4UIcmdWithADoubleAndUnit( "/LYSim/SetBeamZ", this );
  SetBeamZCmd->SetGuidance( "Set central Z position of the beam" );
  SetBeamZCmd->SetParameterName( "BeamZ", false );
  SetBeamZCmd->SetUnitCategory( "Length" );
  SetBeamZCmd->SetDefaultUnit( "mm" );
  SetBeamZCmd->AvailableForStates( G4State_PreInit, G4State_Idle );

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
  delete SetBeamZCmd;
  delete SetWidthCmd;
}

void
LYSimProtonGeneratorMessenger::SetNewValue( G4UIcommand* command, G4String val )
{
  if( command == SetBeamXCmd ){
    genaction->SetBeamX( G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( val ) );
  } else if( command == SetBeamZCmd ){
    genaction->SetBeamZ( G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( val ) );
  } else if( command == SetWidthCmd ){
    genaction->SetWidth( G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( val ) );
  }
}
