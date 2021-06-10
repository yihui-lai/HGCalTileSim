#ifdef CMSSW_GIT_HASH
#include "HGCalTileSim/Tile/interface/LYSimAnalysis.hh"
#include "HGCalTileSim/Tile/interface/LYSimDetectorConstruction.hh"
#include "HGCalTileSim/Tile/interface/LYSimPrimaryGeneratorAction.hh"
#include "HGCalTileSim/Tile/interface/LYSimPrimaryGeneratorMessenger.hh"
#include "HGCalTileSim/Tile/interface/ProjectPath.hh"
#else
#include "LYSimAnalysis.hh"
#include "LYSimDetectorConstruction.hh"
#include "LYSimPrimaryGeneratorAction.hh"
#include "LYSimPrimaryGeneratorMessenger.hh"
#include "ProjectPath.hh"
#endif

#include "Randomize.hh"

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTypes.hh"
#include "G4SystemOfUnits.hh"
#include "G4Poisson.hh"

// Static helper functions.
static double CalcNumPhotons( const double thickness, const int mm );

LYSimPrimaryGeneratorAction::LYSimPrimaryGeneratorAction( LYSimDetectorConstruction* det ) :
  particleSource( new G4GeneralParticleSource() ),
  fDetector( det ),
  messenger( new LYSimPrimaryGeneratorMessenger( this ) ),
  _beamx( 0 ),
  _beamy( 0 ),
  _width( 0 ),
  _photon_multiplier( 1 ),
  _open_angle( CLHEP::pi )
{
  // For our use case, we will need to randomize on a event by event basis
  particleSource->SetMultipleVertex( true );
}

LYSimPrimaryGeneratorAction::~LYSimPrimaryGeneratorAction()
{
  delete particleSource;
}

void
LYSimPrimaryGeneratorAction::GeneratePrimaries( G4Event* anEvent )
{
  RandomizePosition();

  assert( particleSource->GetParticleDefinition()->GetParticleName()
    == "opticalphoton" );

  // We still need to generate the primary vertex.
  particleSource->GeneratePrimaryVertex( anEvent );
}

void
LYSimPrimaryGeneratorAction::RandomizePosition()
{
  /**
   * This particle source will be directly controlled via the /gps commands In
   * user scripts and the ApplyCommand calls in the main function. We are going
   * to set up the common stuff here. For setting the via the interactive
   * session, look at:
   * http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/
   * ForApplicationDeveloper/html/GettingStarted/generalParticleSource.html
   */
  particleSource->ClearAll();
  // Randomizing position of the the source for each event.
  const double x = ( 2*G4UniformRand()-1 )*_width + _beamx;
  const double y = ( 2*G4UniformRand()-1 )*_width + _beamy;
  const double z = fDetector->GetTileZ();
  const double tgap = fDetector->Gettgap();
  const int tn = fDetector->Gettn();
  const double t = fDetector->LocalTileZ( x, y );
  const int mmaterial = fDetector->GetMaterial();
  const int np   = G4Poisson( _photon_multiplier * CalcNumPhotons( t, mmaterial ) );

  for( int i = 0; i < np; ++i ){
    // Creating a new particle source with just 1 photon
    particleSource->AddaSource( 1 );
    particleSource->GetCurrentSource()->SetNumberOfParticles( 1 );
    particleSource->GetCurrentSource()->SetParticleDefinition(
      G4OpticalPhoton::OpticalPhotonDefinition() );

    // Setting default spacial distribution
    G4SPSPosDistribution* pos = particleSource->GetCurrentSource()->GetPosDist();
    pos->SetPosDisType( "Volume" );
    pos->SetPosDisShape( "Cylinder" );
    pos->SetRadius( 0.0001*CLHEP::mm );
    pos->SetHalfZ( t/2 );
    pos->SetCentreCoords( G4ThreeVector( x, y, -(z*(tn-0.5)/tn + tgap*(tn-1)/2/tn)+ t/2  ) );

    // Setting default angular distrution of particle  (isotropic)
    G4SPSAngDistribution* ang = particleSource->GetCurrentSource()->GetAngDist();
    ang->SetAngDistType( "iso" );
    ang->SetMinTheta( CLHEP::pi - _open_angle );
    ang->SetMaxTheta( CLHEP::pi );

    // Energy distribution.
    G4SPSEneDistribution* ene = particleSource->GetCurrentSource()->GetEneDist();
    ene->SetEnergyDisType( "Arb" );
    if(mmaterial==1){
      ene->ArbEnergyHistoFile( project_base + "/data/PhotonSpectrum_ej200.dat");
    }else{
      ene->ArbEnergyHistoFile( project_base + "/data/PhotonSpectrum_scsn81.dat");
//      ene->ArbEnergyHistoFile( project_base + "/data/PhotonSpectrum_scsn81_new.dat");
    }
    ene->ArbInterpolate( "Lin" );

    // Randomizing the polarization.
    const G4ThreeVector normal(
      G4UniformRand(), G4UniformRand(), G4UniformRand() );
    const G4ThreeVector kphoton
      = particleSource->GetParticleMomentumDirection();
    const G4ThreeVector product = normal.cross( kphoton.unit() );
    particleSource->SetParticlePolarization( product.unit() );
  }

}

double
LYSimPrimaryGeneratorAction::NPhotons() const
{
  const double t = fDetector->LocalTileZ( _beamx, _beamy );
  const int mmaterial = fDetector->GetMaterial();
  return _photon_multiplier * CalcNumPhotons( t, mmaterial );
}

unsigned
LYSimPrimaryGeneratorAction::NSources() const
{
  return particleSource->GetNumberofSource();
}

double
CalcNumPhotons( const double thickness, const int mm )
{
  if(mm==1){
//    return 2.2219*10000 * thickness/10; //ej200 4GeV Muon
//    return 2.4031*10000 * thickness/10; //ej200 14GeV Muon
    return 2.55936*10000 * thickness/10;
  }else{
//    return 2.1462*(0.5/0.64)*10000 * thickness/10; //scsn81 4GeV
//    return 2.3217*(0.5/0.64)*10000 * thickness/10; //scsn81
    return 2.60866*(0.5/0.64)*10000 * thickness/10;
  }
}



