#ifdef CMSSW_GIT_HASH
#include "HGCalTileSim/Tile/interface/LYSimAnalysis.hh"
#include "HGCalTileSim/Tile/interface/LYSimDebugGeneratorAction.hh"
#include "HGCalTileSim/Tile/interface/ProjectPath.hh"
#else
#include "LYSimAnalysis.hh"
#include "LYSimDebugGeneratorAction.hh"
#include "ProjectPath.hh"
#endif

#include "G4GeneralParticleSource.hh"
#include "G4ParticleTypes.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

LYSimDebugGeneratorAction::LYSimDebugGeneratorAction() :
  particleSource( new G4GeneralParticleSource() ),
  messenger( new LYSimDebugGeneratorMessenger( this ) )
{
  std::cout << "[LYSIM] entering LYSIMDebugGeneratorAction" << std::endl;

  /**
   * This particle source will be directly controlled via the /gps commands In
   * user scripts and the ApplyCommand calls in the main function. We are going
   * to set up the common stuff here. For setting the via the interactive
   * session, look at:
   * http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/
   * ForApplicationDeveloper/html/GettingStarted/generalParticleSource.html
   */
  particleSource->SetParticleDefinition(
    G4OpticalPhoton::OpticalPhotonDefinition() );
  particleSource->GetCurrentSource()->SetNumberOfParticles( 1 );

  G4SPSPosDistribution* pos = particleSource->GetCurrentSource()->GetPosDist();
  pos->SetPosDisType( "Point" );
  pos->SetCentreCoords( G4ThreeVector( 0, 0, 0 ) );

  G4SPSAngDistribution* ang = particleSource->GetCurrentSource()->GetAngDist();
  ang->SetAngDistType( "iso" );
  ang->SetMinTheta( 0 );
  ang->SetMaxTheta( 0 );
  ang->SetMinPhi( 0 );
  ang->SetMaxPhi( 0 );

  // Energy distribution.
  G4SPSEneDistribution* ene = particleSource->GetCurrentSource()->GetEneDist();
  ene->SetEnergyDisType( "Arb" );
  ene->ArbEnergyHistoFile( project_base + "/data/PhotonSpectrum_ej200.dat" );
  ene->ArbInterpolate( "Lin" );
}

LYSimDebugGeneratorAction::~LYSimDebugGeneratorAction()
{
  delete particleSource;
}

void
LYSimDebugGeneratorAction::GeneratePrimaries( G4Event* anEvent )
{
  assert( particleSource->GetParticleDefinition()->GetParticleName()
    == "opticalphoton" );

  // We still need to generate the primary vertex.
  particleSource->GeneratePrimaryVertex( anEvent );
}

void
LYSimDebugGeneratorAction::SetParticlePosition( const G4ThreeVector& v )
{
  particleSource->GetCurrentSource()->GetPosDist()->SetCentreCoords( v );
}

void
LYSimDebugGeneratorAction::SetParticleDirection( const G4ThreeVector& v )
{
  G4SPSAngDistribution* ang = particleSource->GetCurrentSource()->GetAngDist();
  ang->SetMinTheta( ( -v ).theta() );
  ang->SetMaxTheta( ( -v ).theta() );
  ang->SetMinPhi( ( -v ).phi() );
  ang->SetMaxPhi( ( -v ).phi() );

  // Randomizing the polarization.
  const double angle = G4UniformRand() * 360.0*deg;
  const G4ThreeVector normal(
    G4UniformRand(), G4UniformRand(), G4UniformRand() );
  const G4ThreeVector kphoton   = v.unit();
  const G4ThreeVector product   = normal.cross( kphoton );
  const double modul2           = product*product;
  const G4ThreeVector e_perpend = ( 1./std::sqrt( modul2 ) )*product;
  const G4ThreeVector e_paralle = e_perpend.cross( kphoton );
  const G4ThreeVector polar     = std::cos( angle )*e_paralle
                                  + std::sin( angle )*e_perpend;
  particleSource->SetParticlePolarization( polar );
}

// ------------------------------------------------------------------------------

#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

LYSimDebugGeneratorMessenger::LYSimDebugGeneratorMessenger( LYSimDebugGeneratorAction* act )
  : genaction( act )
{
  detDir = new G4UIdirectory( "/LYSim/Debug/" );
  detDir->SetGuidance( " Debugging Particle source " );

  SetParticlePosition = new G4UIcmdWith3VectorAndUnit(
    "/LYSim/Debug/SetParticlePosition", this );
  SetParticlePosition->SetGuidance( "Set initial position of debugging photon" );
  SetParticlePosition->SetUnitCategory( "Length" );
  SetParticlePosition->SetDefaultUnit( "mm" );
  SetParticlePosition->AvailableForStates( G4State_PreInit, G4State_Idle );

  SetParticleDirection = new G4UIcmdWith3Vector(
    "/LYSim/Debug/SetParticleDirection", this );
  SetParticleDirection->SetGuidance(
    "Set initial direction of debugging photon" );
  SetParticleDirection->AvailableForStates( G4State_PreInit, G4State_Idle );

}

LYSimDebugGeneratorMessenger::~LYSimDebugGeneratorMessenger()
{
  delete detDir;
  delete SetParticlePosition;
  delete SetParticleDirection;
}

void
LYSimDebugGeneratorMessenger::SetNewValue( G4UIcommand* command, G4String val )
{
  if( command == SetParticlePosition ){
    genaction->SetParticlePosition(
      G4UIcmdWith3VectorAndUnit::GetNew3VectorValue( val ) );
  } else if( command == SetParticleDirection ){
    genaction->SetParticleDirection(
      G4UIcmdWith3Vector::GetNew3VectorValue( val ) );
  }
}
