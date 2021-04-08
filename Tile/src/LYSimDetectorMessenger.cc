#ifdef CMSSW_GIT_HASH
#include "HGCalTileSim/Tile/interface/LYSimDetectorMessenger.hh"
#else
#include "LYSimDetectorMessenger.hh"
#endif

#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIdirectory.hh"

LYSimDetectorMessenger::LYSimDetectorMessenger( LYSimDetectorConstruction* Det )
  : Detector( Det )
{
  detDir = new G4UIdirectory( "/LYSim/" );
  detDir->SetGuidance( " Geometry Setup " );

  UpdateCmd = new G4UIcmdWithoutParameter( "/LYSim/Update", this );
  UpdateCmd->SetGuidance( "Update geometry" );
  UpdateCmd->SetGuidance( "This command MUST be applied before \"beamOn\" " );
  UpdateCmd->SetGuidance( "if you changed geometrical value(s)." );
  UpdateCmd->AvailableForStates( G4State_Idle );

  // Tile Geometry commands
  SetTileXCmd = new G4UIcmdWithADoubleAndUnit( "/LYSim/SetTileX", this );
  SetTileXCmd->SetGuidance( "Set the scintillator transverse dimensions" );
  SetTileXCmd->SetParameterName( "TileX", false );
  SetTileXCmd->SetUnitCategory( "Length" );
  SetTileXCmd->SetDefaultUnit( "mm" );
  SetTileXCmd->SetRange( "TileX>=0." );
  SetTileXCmd->AvailableForStates( G4State_PreInit, G4State_Idle );

  SetTileYCmd = new G4UIcmdWithADoubleAndUnit( "/LYSim/SetTileY", this );
  SetTileYCmd->SetGuidance( "Set the scintillator transverse dimensions" );
  SetTileYCmd->SetParameterName( "TileY", false );
  SetTileYCmd->SetUnitCategory( "Length" );
  SetTileYCmd->SetDefaultUnit( "mm" );
  SetTileYCmd->SetRange( "TileY>=0." );
  SetTileYCmd->AvailableForStates( G4State_PreInit, G4State_Idle );

  SetTileZCmd = new G4UIcmdWithADoubleAndUnit( "/LYSim/SetTileZ", this );
  SetTileZCmd->SetGuidance( "Set the scintillator thickness" );
  SetTileZCmd->SetParameterName( "TileZ", false );
  SetTileZCmd->SetUnitCategory( "Length" );
  SetTileZCmd->SetDefaultUnit( "mm" );
  SetTileZCmd->SetRange( "TileZ>=0." );
  SetTileZCmd->AvailableForStates( G4State_PreInit, G4State_Idle );

  // SiPM geometry commands
  SetSiPMXCmd = new G4UIcmdWithADoubleAndUnit( "/LYSim/SetSiPMX", this );
  SetSiPMXCmd->SetGuidance( "Set the SiPM transverse dimensions" );
  SetSiPMXCmd->SetParameterName( "SiPMX", false );
  SetSiPMXCmd->SetUnitCategory( "Length" );
  SetSiPMXCmd->SetDefaultUnit( "mm" );
  SetSiPMXCmd->SetRange( "SiPMX>=0." );
  SetSiPMXCmd->AvailableForStates( G4State_PreInit, G4State_Idle );

  SetSiPMYCmd = new G4UIcmdWithADoubleAndUnit( "/LYSim/SetSiPMY", this );
  SetSiPMYCmd->SetGuidance( "Set the SiPM transverse dimensions" );
  SetSiPMYCmd->SetParameterName( "SiPMY", false );
  SetSiPMYCmd->SetUnitCategory( "Length" );
  SetSiPMYCmd->SetDefaultUnit( "mm" );
  SetSiPMYCmd->SetRange( "SiPMY>=0." );
  SetSiPMYCmd->AvailableForStates( G4State_PreInit, G4State_Idle );

  SetSiPMRimCmd = new G4UIcmdWithADoubleAndUnit( "/LYSim/SetSiPMRim", this );
  SetSiPMRimCmd->SetGuidance( "Set the SiPM transverse dimensions" );
  SetSiPMRimCmd->SetParameterName( "SiPMRim", false );
  SetSiPMRimCmd->SetUnitCategory( "Length" );
  SetSiPMRimCmd->SetDefaultUnit( "mm" );
  SetSiPMRimCmd->SetRange( "SiPMRim>=0." );
  SetSiPMRimCmd->AvailableForStates( G4State_PreInit, G4State_Idle );

  SetSiPMStandCmd = new G4UIcmdWithADoubleAndUnit( "/LYSim/SetSiPMStand", this );
  SetSiPMStandCmd->SetGuidance( "Set the SiPM transverse dimensions" );
  SetSiPMStandCmd->SetParameterName( "SiPMStand", false );
  SetSiPMStandCmd->SetUnitCategory( "Length" );
  SetSiPMStandCmd->SetDefaultUnit( "mm" );
  SetSiPMStandCmd->SetRange( "SiPMStand>=0." );
  SetSiPMStandCmd->AvailableForStates( G4State_PreInit, G4State_Idle );

  // PCB Hole commands
  SetPCBRadiusCmd = new G4UIcmdWithADoubleAndUnit( "/LYSim/SetPCBRadius", this );
  SetPCBRadiusCmd->SetGuidance( "Set PCB Hole dimensions" );
  SetPCBRadiusCmd->SetParameterName( "PCBRadius", false );
  SetPCBRadiusCmd->SetUnitCategory( "Length" );
  SetPCBRadiusCmd->SetDefaultUnit( "mm" );
  SetPCBRadiusCmd->SetRange( "PCBRadius>=0." );
  SetPCBRadiusCmd->AvailableForStates( G4State_PreInit, G4State_Idle );

  SetPCBReflectCmd = new G4UIcmdWithADouble( "/LYSim/SetPCBReflect", this );
  SetPCBReflectCmd->SetGuidance( "Set reflectivity of exposed PCB" );
  SetPCBReflectCmd->SetParameterName( "PCBReflect", false );
  SetPCBReflectCmd->SetRange( "PCBReflect>=0." );
  SetPCBReflectCmd->AvailableForStates( G4State_PreInit, G4State_Idle );

  // Optical setting commands
  SetTileAbsMultCmd = new G4UIcmdWithADouble( "/LYSim/SetTileAbsMult", this );
  SetTileAbsMultCmd->SetGuidance( "Set the light absorption length multiplier" );
  SetTileAbsMultCmd->SetParameterName( "TileAbsMult", false );
  SetTileAbsMultCmd->SetRange( "TileAbsMult>=0." );
  SetTileAbsMultCmd->AvailableForStates( G4State_PreInit, G4State_Idle );

  SetWrapReflectCmd = new G4UIcmdWithADouble( "/LYSim/SetWrapReflect", this );
  SetWrapReflectCmd->SetGuidance( "Set Wrap reflectivity" );
  SetWrapReflectCmd->SetParameterName( "WrapReflect", false );
  SetWrapReflectCmd->SetRange( "WrapReflect>=0." );
  SetWrapReflectCmd->AvailableForStates( G4State_PreInit, G4State_Idle );

  SetTileAlphaCmd = new G4UIcmdWithADouble( "/LYSim/SetTileAlpha", this );
  SetTileAlphaCmd->SetGuidance( "Set micro facet distribution of tile" );
  SetTileAlphaCmd->SetParameterName( "TileAlpha", false );
  SetTileAlphaCmd->SetRange( "TileAlpha>=0." );
  SetTileAlphaCmd->AvailableForStates( G4State_PreInit, G4State_Idle );

  SetDimpleAlphaCmd = new G4UIcmdWithADouble( "/LYSim/SetDimpleAlpha", this );
  SetDimpleAlphaCmd->SetGuidance( "Set micro facet distribution of dimple" );
  SetDimpleAlphaCmd->SetParameterName( "DimpleAlpha", false );
  SetDimpleAlphaCmd->SetRange( "DimpleAlpha>=0." );
  SetDimpleAlphaCmd->AvailableForStates( G4State_PreInit, G4State_Idle );

  // Dimple Geometry commands
  SetDimpleIndentCmd = new G4UIcmdWithADoubleAndUnit(
    "/LYSim/SetDimpleIndent", this );
  SetDimpleIndentCmd->SetGuidance(
    "Setting the reccess depth of the dimple [mm]" );
  SetDimpleIndentCmd->SetParameterName( "DimpleIndent", false );
  SetDimpleIndentCmd->SetUnitCategory( "Length" );
  SetDimpleIndentCmd->SetDefaultUnit( "mm" );
  SetDimpleIndentCmd->SetRange( "DimpleIndent>=0." );
  SetDimpleIndentCmd->AvailableForStates( G4State_PreInit, G4State_Idle );

  SetDimpleRadiusCmd = new G4UIcmdWithADoubleAndUnit(
    "/LYSim/SetDimpleRadius", this );
  SetDimpleRadiusCmd->SetGuidance(
    "Setting the reccess depth of the dimple [mm]" );
  SetDimpleRadiusCmd->SetParameterName( "DimpleIndent", false );
  SetDimpleRadiusCmd->SetUnitCategory( "Length" );
  SetDimpleRadiusCmd->SetDefaultUnit( "mm" );
  SetDimpleRadiusCmd->SetRange( "DimpleIndent>=0." );
  SetDimpleRadiusCmd->AvailableForStates( G4State_PreInit, G4State_Idle );

  SetDimpleTypeCmd = new G4UIcmdWithAString( "/LYSim/SetDimpleType", this );
  SetDimpleTypeCmd->SetGuidance( "Setting the type of dimple to use" );
  SetDimpleTypeCmd->SetCandidates( "Spherical FlatDome Cylindrical Ellipsoid" );
  SetDimpleTypeCmd->AvailableForStates( G4State_PreInit, G4State_Idle );
}

LYSimDetectorMessenger::~LYSimDetectorMessenger()
{
  delete detDir;
  delete UpdateCmd;
  delete SetTileZCmd;
  delete SetTileXCmd;
  delete SetTileYCmd;

  delete SetTileAbsMultCmd;
  delete SetWrapReflectCmd;
  delete SetTileAlphaCmd;
  delete SetDimpleAlphaCmd;

  delete SetDimpleIndentCmd;
  delete SetDimpleRadiusCmd;
  delete SetDimpleTypeCmd;
}

void
LYSimDetectorMessenger::SetNewValue( G4UIcommand* command, G4String val )
{
  if( command == UpdateCmd ){
    Detector->UpdateGeometry();
  } else if( command == SetTileZCmd ){
    Detector->
    SetTileZ( G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( val ) );
  } else if( command == SetTileXCmd ){
    Detector->
    SetTileX( G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( val ) );
  } else if( command == SetTileYCmd ){
    Detector->SetTileY(
      G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( val ) );
  } else if( command == SetSiPMXCmd ){
    Detector->SetSiPMX( G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( val ) );
  } else if( command == SetSiPMYCmd ){
    Detector->SetSiPMY( G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( val ) );
  } else if( command == SetSiPMRimCmd ){
    Detector->SetSiPMRim( G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( val ) );
  } else if( command == SetSiPMStandCmd ){
    Detector->SetSiPMStand(
      G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( val ) );
  } else if( command == SetPCBRadiusCmd ){
    Detector->SetPCBRadius(
      G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( val ) );
  } else if( command == SetPCBReflectCmd ){
    Detector->SetSiPMReflect( G4UIcmdWithADouble::GetNewDoubleValue( val ) );
  } else if( command == SetTileAbsMultCmd ){
    Detector->SetTileAbsMult( G4UIcmdWithADouble::GetNewDoubleValue( val ) );
  } else if( command == SetWrapReflectCmd ){
    Detector->SetWrapReflect( G4UIcmdWithADouble::GetNewDoubleValue( val ) );
  } else if( command == SetTileAlphaCmd ){
    Detector->SetTileAlpha( G4UIcmdWithADouble::GetNewDoubleValue( val ) );
  } else if( command == SetDimpleAlphaCmd ){
    Detector->SetDimpleAlpha( G4UIcmdWithADouble::GetNewDoubleValue( val ) );
  } else if( command == SetDimpleIndentCmd ){
    Detector->SetDimpleIndent(
      G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( val ) );
  } else if( command == SetDimpleRadiusCmd ){
    Detector->SetDimpleRadius(
      G4UIcmdWithADoubleAndUnit::GetNewDoubleValue( val ) );
  } else if( command == SetDimpleTypeCmd ){
    if( val == "Spherical" ){
      Detector->SetDimpleType( LYSimDetectorConstruction::SPHERICAL );
    } else if( val == "FlatDome" ){
      Detector->SetDimpleType( LYSimDetectorConstruction::FLAT_DOME );
    } else if( val == "Cylindrical" ){
      Detector->SetDimpleType( LYSimDetectorConstruction::CYLINDRICAL );
    } else if( val == "Ellipsoid" ){
      Detector->SetDimpleType( LYSimDetectorConstruction::ELLIPSOID );
    }
  }
}
