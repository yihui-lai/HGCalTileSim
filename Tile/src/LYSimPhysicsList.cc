#ifdef CMSSW_GIT_HASH
#include "HGCalTileSim/Tile/interface/LYSimPhysicsList.hh"
#include "HGCalTileSim/Tile/interface/LYSimScintillation.hh"
#else
#include "LYSimPhysicsList.hh"
#include "LYSimScintillation.hh"
#endif

#include "globals.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"

#include "G4ProcessManager.hh"

#include "G4Cerenkov.hh"
#include "G4OpAbsorption.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpMieHG.hh"
#include "G4OpRayleigh.hh"
#include "G4OpWLS.hh"

#include "G4EmSaturation.hh"
#include "G4LossTableManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

LYSimPhysicsList::LYSimPhysicsList() :  G4VUserPhysicsList()
{
  defaultCutValue              = 0.01*mm;
  theCerenkovProcess           = NULL;
  theScintillationProcess      = NULL;
  theAbsorptionProcess         = NULL;
  theRayleighScatteringProcess = NULL;
  theMieHGScatteringProcess    = NULL;
  theBoundaryProcess           = NULL;
  theWLSProcess                = NULL;

  SetVerboseLevel( 0 );
}

LYSimPhysicsList::~LYSimPhysicsList(){}

void
LYSimPhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called for all particles
  // which you want to use. This ensures that objects of these particle types
  // will be created in the program.

  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBaryons();
}

void
LYSimPhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
  G4Gamma::GammaDefinition();

  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
}

void
LYSimPhysicsList::ConstructLeptons()
{
  // leptons
  //  e+/-
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  // mu+/-
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();
  // nu_e
  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  // nu_mu
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}

void
LYSimPhysicsList::ConstructMesons()
{
  //  mesons
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
}

void
LYSimPhysicsList::ConstructBaryons()
{
  //  barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();

  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();
}

void
LYSimPhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructGeneral();
  ConstructEM();
  ConstructOp();
}

#include "G4Decay.hh"

void
LYSimPhysicsList::ConstructGeneral()
{
  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();
  G4ParticleTable::G4PTblDicIterator* theParticleIterator
    = theParticleTable->GetIterator();
  theParticleIterator->reset();

  while( ( *theParticleIterator )() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager     = particle->GetProcessManager();
    if( theDecayProcess->IsApplicable( *particle ) ){
      pmanager->AddProcess( theDecayProcess );
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager->SetProcessOrdering( theDecayProcess, idxPostStep );
      pmanager->SetProcessOrdering( theDecayProcess, idxAtRest );
    }
  }
}

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"

#include "G4eBremsstrahlung.hh"
#include "G4eIonisation.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuBremsstrahlung.hh"
#include "G4MuIonisation.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"

void
LYSimPhysicsList::ConstructEM()
{
  G4ParticleTable::G4PTblDicIterator* theParticleIterator
    = theParticleTable->GetIterator();
  theParticleIterator->reset();

  while( ( *theParticleIterator )() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager     = particle->GetProcessManager();
    G4String particleName          = particle->GetParticleName();

    if( particleName == "gamma" ){
      // gamma
      // Construct processes for gamma
      pmanager->AddDiscreteProcess( new G4GammaConversion() );
      pmanager->AddDiscreteProcess( new G4ComptonScattering() );
      pmanager->AddDiscreteProcess( new G4PhotoElectricEffect() );

    } else if( particleName == "e-" ){
      // electron
      // Construct processes for electron
      pmanager->AddProcess( new G4eMultipleScattering(), -1, 1, 1 );
      pmanager->AddProcess( new G4eIonisation(),         -1, 2, 2 );
      pmanager->AddProcess( new G4eBremsstrahlung(),     -1, 3, 3 );

    } else if( particleName == "e+" ){
      // positron
      // Construct processes for positron
      pmanager->AddProcess( new G4eMultipleScattering(), -1, 1, 1 );
      pmanager->AddProcess( new G4eIonisation(),         -1, 2, 2 );
      pmanager->AddProcess( new G4eBremsstrahlung(),     -1, 3, 3 );
      pmanager->AddProcess( new G4eplusAnnihilation(),   0, -1, 4 );

    } else if( particleName == "mu+" ||
               particleName == "mu-"    ){
      // muon
      // Construct processes for muon
      pmanager->AddProcess( new G4MuMultipleScattering(), -1, 1, 1 );
      pmanager->AddProcess( new G4MuIonisation(),         -1, 2, 2 );
      pmanager->AddProcess( new G4MuBremsstrahlung(),     -1, 3, 3 );
      pmanager->AddProcess( new G4MuPairProduction(),     -1, 4, 4 );

    } else {
      if( ( particle->GetPDGCharge() != 0.0 ) &&
          ( particle->GetParticleName() != "chargedgeantino" ) ){
        // all others charged particles except geantino
        pmanager->AddProcess( new G4hMultipleScattering(), -1, 1, 1 );
        pmanager->AddProcess( new G4hIonisation(),         -1, 2, 2 );
      }
    }
  }
}

void
LYSimPhysicsList::ConstructOp()
{
  theCerenkovProcess           = new G4Cerenkov( "Cerenkov" );
  theScintillationProcess      = new LYSimScintillation( "Scintillation" );
  theAbsorptionProcess         = new G4OpAbsorption();
  theRayleighScatteringProcess = new G4OpRayleigh();
  theMieHGScatteringProcess    = new G4OpMieHG();
  theBoundaryProcess           = new G4OpBoundaryProcess();
  theWLSProcess                = new G4OpWLS();

  // Set verbose level to 0
  SetVerbose( 1 );

  theCerenkovProcess->SetMaxNumPhotonsPerStep( 20 );
  theCerenkovProcess->SetMaxBetaChangePerStep( 10.0 );
  theCerenkovProcess->SetTrackSecondariesFirst( true );

  theScintillationProcess->SetScintillationYieldFactor( 1. );
  theScintillationProcess->SetTrackSecondariesFirst( true );

  theWLSProcess->UseTimeProfile( "delta" );

  G4ParticleTable::G4PTblDicIterator* theParticleIterator
    = theParticleTable->GetIterator();
  theParticleIterator->reset();

  while( ( *theParticleIterator )() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager     = particle->GetProcessManager();
    G4String particleName          = particle->GetParticleName();
//    if( theCerenkovProcess->IsApplicable( *particle ) ){
//      pmanager->AddProcess( theCerenkovProcess );
//      pmanager->SetProcessOrdering( theCerenkovProcess, idxPostStep );
//    }
    if( theScintillationProcess->IsApplicable( *particle ) ){
      pmanager->AddProcess( theScintillationProcess );
      pmanager->SetProcessOrderingToLast( theScintillationProcess, idxAtRest );
      pmanager->SetProcessOrderingToLast( theScintillationProcess, idxPostStep );
    }
    if( particleName == "opticalphoton" ){
      pmanager->AddDiscreteProcess(         theAbsorptionProcess );
      pmanager->AddDiscreteProcess( theRayleighScatteringProcess );
      pmanager->AddDiscreteProcess(    theMieHGScatteringProcess );
      pmanager->AddDiscreteProcess(           theBoundaryProcess );
      pmanager->AddDiscreteProcess(                theWLSProcess );
    }
  }
}


void
LYSimPhysicsList::ConstructIdealOp()
{
  theCerenkovProcess           = new G4Cerenkov( "Cerenkov" );
  theScintillationProcess      = new LYSimScintillation( "Scintillation" );
  theAbsorptionProcess         = new G4OpAbsorption();
  theRayleighScatteringProcess = new G4OpRayleigh();
  theMieHGScatteringProcess    = new G4OpMieHG();
  theBoundaryProcess           = new G4OpBoundaryProcess();
  theWLSProcess                = new G4OpWLS();

  // Verbose level set to 0
  SetVerbose( 0 );

  theCerenkovProcess->SetMaxNumPhotonsPerStep( 20 );
  theCerenkovProcess->SetMaxBetaChangePerStep( 10.0 );
  theCerenkovProcess->SetTrackSecondariesFirst( true );

  theScintillationProcess->SetScintillationYieldFactor( 1. );
  theScintillationProcess->SetTrackSecondariesFirst( true );

  G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
  theScintillationProcess->AddSaturation( emSaturation );

  theWLSProcess->UseTimeProfile( "delta" );

  G4ParticleTable::G4PTblDicIterator* theParticleIterator
    = theParticleTable->GetIterator();
  theParticleIterator->reset();

  while( ( *theParticleIterator )() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager     = particle->GetProcessManager();
    G4String particleName          = particle->GetParticleName();
    if( particleName == "opticalphoton" ){
      pmanager->AddDiscreteProcess( theBoundaryProcess );
      pmanager->AddDiscreteProcess(      theWLSProcess );
    }
  }
}

void
LYSimPhysicsList::SetVerbose( G4int verbose )
{
  theCerenkovProcess->SetVerboseLevel( verbose );
  theScintillationProcess->SetVerboseLevel( verbose );
  theAbsorptionProcess->SetVerboseLevel( verbose );
  theRayleighScatteringProcess->SetVerboseLevel( verbose );
  theMieHGScatteringProcess->SetVerboseLevel( verbose );
  theBoundaryProcess->SetVerboseLevel( verbose );
  theWLSProcess->SetVerboseLevel( 1 );
}

void
LYSimPhysicsList::SetNbOfPhotonsCerenkov( G4int MaxNumber )
{
  theCerenkovProcess->SetMaxNumPhotonsPerStep( MaxNumber );
}

void
LYSimPhysicsList::SetCuts()
{
  SetCutsWithDefault();
}
