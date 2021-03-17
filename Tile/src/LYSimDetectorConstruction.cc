#ifdef CMSSW_GIT_HASH
#include "HGCalTileSim/Tile/interface/LYSimDetectorConstruction.hh"
#include "HGCalTileSim/Tile/interface/LYSimDetectorMessenger.hh"
#include "HGCalTileSim/Tile/interface/Materials.hh"
#include "HGCalTileSim/Tile/interface/SurfaceProperty.hh"
#else
#include "LYSimDetectorConstruction.hh"
#include "LYSimDetectorMessenger.hh"
#include "Materials.hh"
#include "SurfaceProperty.hh"
#endif

#include <math.h>

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Ellipsoid.hh"
#include "G4GeometryManager.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Orb.hh"
#include "G4Orb.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4PVPlacement.hh"
#include "G4RegionStore.hh"
#include "G4RotationMatrix.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SolidStore.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4Torus.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"
#include "G4VPhysicalVolume.hh"

using std::cos;
using std::sin;
using std::tan;
using std::atan;
using std::exp;
using namespace CLHEP;

LYSimPMTSD* LYSimDetectorConstruction::fPMTSD = NULL;

LYSimDetectorConstruction::LYSimDetectorConstruction()
  : G4VUserDetectorConstruction()
{
  fdetectorMessenger = new LYSimDetectorMessenger( this );

  _tilex        = 50*mm;
  _tiley        = 10*mm;
  _tilez        = 50.0*mm;  //200*mm
  _tile_x1      = 0.0*mm;
  _tile_x2      = 0.0*mm;
  wrapgap       = 0.1*mm;
  wrapthickness = 0.1*mm;

  _absmult      = 1000; //SetTileAbsMult to 1m
  _ScintiN      = 10; //Set scintillation to 10 /keV
  _wrap_reflect = 0.985;
  _tile_alpha   = 0.01;
  _dimple_alpha = 0.1;

  _sipm_deadwidth  = 0.2*mm;
  _sipm_x          = 1.4*mm;
  _sipm_y          = 1.4*mm;
  _sipm_z          = 0.4*mm;
  _sipm_rimwidth   = 0.3*mm;
  _sipm_glasswidth = 0.1*mm;
  _sipm_standz     = 0.3*mm;
  _sipm_x          = 10*mm;
  _sipm_y          = 10*mm;


  // Default Dimple settings
  _dimple_type   = SPHERICAL;
  _dimple_indent = 1.0*mm;
  _dimple_radius = 6.0*mm;// 3.4409*mm
  
  // Default Hole settings
  _pcb_radius       = 2.5;
  _pcb_reflectivity = 0.8;


  // Defining material list.
  fBialkali = Make_Bialkali();
  fEpoxy    = Make_Epoxy();
  fAir      = Make_Custom_Air();
  fEJ200    = Make_EJ200();
  fResin    = Make_Resin();
  SetTileAbsMult( _absmult );
  SetTileScintillation( _ScintiN ); // N/keV

  // Defining surface list.
  fESROpSurface           = MakeS_RoughMirror();
  //fIdealPolishedOpSurface = MakeS_IdealPolished();
  fTileBulkSurface        = MakeS_RoughInterface( _tile_alpha );
  //fTileDimpleSurface      = MakeS_RoughInterface( _dimple_alpha );
  //fIdealWhiteOpSurface    = MakeS_IdealWhiteSurface();
  fSiPMSurface            = MakeS_SiPM();
  //fPCBSurface             = MakeS_PCBSurface();

  SetWrapReflect( _wrap_reflect );

//wls
_handwrap   = true;

_hole_radius = 1.0*mm;
_hole_x1 = -12.5*mm;
_hole_x2 = 12.5*mm;

_WLSfiberR = 0.7*mm;
_WLSfiber_clad_thick = 0.05*mm;

_WLSfiberZ = 5.2*m;
_WLS_zoff = 1.7*m;

_WLSfiberZ = 5.2*cm;
_WLS_zoff = 1.7*cm;
mfiber  = Make_Y11();
mfiber_clad = Make_Pethylene();
fcoating = Make_Coating();
fTiO2Surface = MakeS_TiO2Surface();


}

void
LYSimDetectorConstruction::UpdateGeometry()
{
  // clean-up previous geometry
  G4GeometryManager::GetInstance()->OpenGeometry();

  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  G4LogicalSkinSurface::CleanSurfaceTable();
  G4LogicalBorderSurface::CleanSurfaceTable();

  G4RunManager::GetRunManager()->DefineWorldVolume( Construct() );
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

LYSimDetectorConstruction::~LYSimDetectorConstruction()
{
  if( fdetectorMessenger ){ delete fdetectorMessenger; }
}

G4VPhysicalVolume*
LYSimDetectorConstruction::Construct()
{
  static const bool checkOverlaps = true;

  ///////////////////////////////////////////////////////////////////////////////
  // World volume
  ///////////////////////////////////////////////////////////////////////////////
  G4Box* solidWorld = new G4Box( "World", WorldHalfX()
                               , WorldHalfY(), WorldHalfZ() );

  G4LogicalVolume* logicWorld = new G4LogicalVolume( solidWorld, fAir, "World" );

  G4VPhysicalVolume* physWorld = new G4PVPlacement( 0
                                                  , G4ThreeVector( 0, 0, 0 )
                                                  , logicWorld
                                                  , "World"
                                                  , 0
                                                  , false
                                                  , 0
                                                  , checkOverlaps );
  ///////////////////////////////////////////////////////////////////////////////
  // wrapping or clad
  ///////////////////////////////////////////////////////////////////////////////
  G4VSolid* solidHoleBound = new G4Tubs( "TileDimpleBox"
                                         , 0, _hole_radius
                                         , _tilez*1.5
                                         , 0, 2*pi  );
  G4LogicalVolume* logicWrap;

  if(_handwrap){
  G4VSolid* solidWrap0 = ConstructHollowWrapSolid();
  G4VSolid* solidWrap1 = new G4SubtractionSolid( "solidWrap"
                            , solidWrap0, solidHoleBound
                            , 0, G4ThreeVector( _hole_x1, 0, 0 ) );
  G4VSolid* solidWrap = new G4SubtractionSolid( "solidWrap"
                            , solidWrap1, solidHoleBound
                            , 0, G4ThreeVector( _hole_x2, 0, 0 ) );
  logicWrap = new G4LogicalVolume( solidWrap, fEpoxy,  "Wrap" );

  G4VPhysicalVolume* physWrap = new G4PVPlacement( 0
                                                 , G4ThreeVector( 0, 0, 0 )
                                                 , logicWrap
                                                 , "Wrap"
                                                 , logicWorld
                                                 , false
                                                 , 0
                                                 , checkOverlaps );
  G4LogicalSkinSurface* WrapSurface =
    new G4LogicalSkinSurface( "WrapSurface"
                              , logicWrap, fESROpSurface );

  }else{
  G4VSolid* solidExtrusion = ConstructHollowWrapCladSolid();

  logicWrap =
    new G4LogicalVolume(solidExtrusion, fcoating, "Extrusion");

  G4VPhysicalVolume* physWrap = new G4PVPlacement( 0
                                                 , G4ThreeVector( 0, 0, 0 )
                                                 , logicWrap
                                                 , "Wrap"
                                                 , logicWorld
                                                 , false
                                                 , 0
                                                 , checkOverlaps );
                                                 
  G4LogicalSkinSurface* WrapSurface =
    new G4LogicalSkinSurface( "WrapSurface"
                              , logicWrap, fTiO2Surface );
  }
  ///////////////////////////////////////////////////////////////////////////////
  // extruded scintillator
  ///////////////////////////////////////////////////////////////////////////////
  G4VSolid* solidTile = ConstructTrapazoidSolid( "TileTrap"
                                               , _tilex
                                               , _tiley
                                               , _tilez
                                               , _tile_x1
                                               , _tile_x2 );
  G4VSolid* tileBulk1
    = new G4SubtractionSolid( "TileSolid_Bulk"
                            , solidTile, solidHoleBound
                            , 0, G4ThreeVector( _hole_x1, 0, 0 ) );
  G4VSolid* tileBulk
    = new G4SubtractionSolid( "TileSolid_Bulk"
                            , tileBulk1, solidHoleBound
                            , 0, G4ThreeVector( _hole_x2, 0, 0 ) );

  G4LogicalVolume* logicTileBulk = new G4LogicalVolume( tileBulk
                                                      , fEJ200
                                                      , "TileBulkLogic" );
  G4VPhysicalVolume* physTileBulk = new G4PVPlacement( 0
                                                     , G4ThreeVector( 0, 0, 0 )
                                                     , logicTileBulk
                                                     , "TileBulkPhysic"
                                                     , logicWorld
                                                     , false
                                                     , 0
                                                     , checkOverlaps );

  ///////////////////////////////////////////////////////////////////////////////
  // fiber. single layer clad
  ///////////////////////////////////////////////////////////////////////////////
  assert(_hole_radius>=_WLSfiberR+_WLSfiber_clad_thick); 
  G4VSolid* solidWLSfiber = new G4Tubs("WLSFiber", 0., _WLSfiberR, _WLSfiberZ, 0., 2*pi);
  G4VSolid* solidWLSfiber_clad = new G4Tubs("WLSFiber_clad", _WLSfiberR, _WLSfiberR+_WLSfiber_clad_thick, _WLSfiberZ, 0., 2*pi);
  G4LogicalVolume* logicWLSfiber = new G4LogicalVolume( solidWLSfiber , mfiber,  "logicWLSfiber" );
  G4LogicalVolume* logicWLSfiber_clad = new G4LogicalVolume( solidWLSfiber_clad , mfiber_clad,  "logicWLSfiber_clad" );
  G4VPhysicalVolume* physWLSfiber = new G4PVPlacement( 0, G4ThreeVector(_hole_x1, 0, -_WLS_zoff)
                                                      , logicWLSfiber
                                                      , "PhyhWLSfiber"
                                                      , logicWorld
                                                      , false
                                                      , 0
                                                      , checkOverlaps );
  G4VPhysicalVolume* physWLSfiber_clad = new G4PVPlacement( 0, G4ThreeVector(_hole_x1, 0, -_WLS_zoff)
                                                      , logicWLSfiber_clad
                                                      , "PhyhWLSfiber_cald"
                                                      , logicWorld
                                                      , false
                                                      , 0
                                                      , checkOverlaps );
  //TODO: surface properties
  /*
  G4LogicalSkinSurface* FiberSurface =
    new G4LogicalSkinSurface( "FiberSurface"
                              , logicWLSfiber_clad, fTiO2Surface );  

        new G4LogicalBorderSurface("surfaceClad1Out", physWLSfiber_clad,
                                   physWorld, opSurface);
        new G4LogicalBorderSurface("surfaceClad1In", physWorld, physWLSfiber_clad,
                                   opSurface);
                                   */
  ///////////////////////////////////////////////////////////////////////////////
  // realistic SiPM
  ///////////////////////////////////////////////////////////////////////////////
  
  /*
  G4Box* solidSiPMDead = new G4Box( "SiPMDead"
                                  , 0.5*_sipm_deadwidth, 0.5*_sipm_deadwidth
                                  , _sipm_z );

  G4Box* solidSiPMInnerBox = new G4Box( "SiPMInnerBox"
                                      , 0.5*_sipm_x, 0.5*_sipm_y,  0.8*_sipm_z );

  G4Box* solidSiPMOuter = new G4Box( "SiPMOuter"
                                   , 0.5*_sipm_x + _sipm_rimwidth
                                   , 0.5*_sipm_y + _sipm_rimwidth
                                   , 0.5*_sipm_z );
  G4Box* solidSiPMStand
    = new G4Box( "SiPMStand"
               , 0.5*_sipm_x+_sipm_rimwidth + _sipm_glasswidth
               , 0.5*_sipm_y+_sipm_rimwidth + _sipm_glasswidth
               , 0.5*_sipm_standz );

  G4Box* solidSiPMResinOuter
    = new G4Box( "SiPMResinOuter"
               , 0.5*_sipm_x + _sipm_rimwidth + _sipm_glasswidth
               , 0.5*_sipm_y + _sipm_rimwidth + _sipm_glasswidth
               , 0.5*_sipm_z + _sipm_glasswidth );

  G4VSolid* solidSiPMSubtract
    = new G4SubtractionSolid( "SiPMSubtract"
                            ,  solidSiPMInnerBox, solidSiPMDead
                            ,   0, G4ThreeVector( 0, 0, 0 ) );
  G4VSolid* solidSiPMCase
    = new G4SubtractionSolid( "SiPMCase"
                            , solidSiPMOuter, solidSiPMSubtract
                            , 0
                            , G4ThreeVector( 0, 0, -0.65 * _sipm_z ) );

  G4VSolid* solidSiPMInner
    = new G4IntersectionSolid( "SiPMInner"
                             , solidSiPMOuter, solidSiPMSubtract
                             , 0
                             , G4ThreeVector( 0, 0, -0.65*_sipm_z ) );

  G4VSolid* solidSiPMResin
    = new G4SubtractionSolid( "SiPMResin"
                            , solidSiPMResinOuter, solidSiPMOuter
                            ,  0
                            , G4ThreeVector( 0, 0, _sipm_glasswidth ) );


  G4LogicalVolume* logicSiPM = new G4LogicalVolume( solidSiPMInner
                                                  , fBialkali,  "SiPM" );

  G4LogicalVolume* logicSiPMCase = new G4LogicalVolume( solidSiPMCase
                                                      , fEpoxy, "SiPMBack" );

  G4LogicalVolume* logicSiPMResin = new G4LogicalVolume( solidSiPMResin
                                                       , fResin, "SiPMResin" );

  G4LogicalVolume* logicSiPMStand = new G4LogicalVolume( solidSiPMStand
                                                       , fEpoxy, "SiPMStand" );
  
  double Resinz = _WLSfiberZ - _WLS_zoff + 0.5*_sipm_z + _sipm_glasswidth;
  const G4ThreeVector ResinOffset( _hole_x1, 0, Resinz );

  const G4ThreeVector SiPMOffset( _hole_x1, 0, Resinz + _sipm_glasswidth);
  const G4ThreeVector StandOffset( _hole_x1, 0, Resinz + 0.5*_sipm_z + 0.5*_sipm_standz + _sipm_glasswidth);

  G4VPhysicalVolume* physSiPMStand = new G4PVPlacement( 0, StandOffset
                                                      , logicSiPMStand
                                                      , "SiPMStand"
                                                      , logicWorld
                                                      , false
                                                      , 0
                                                      , checkOverlaps );

  G4VPhysicalVolume* physSiPMCase = new G4PVPlacement( 0, SiPMOffset
                                                     , logicSiPMCase
                                                     , "Case"
                                                     , logicWorld
                                                     , false
                                                     , 0
                                                     , checkOverlaps );

  G4VPhysicalVolume* physSiPMResin = new G4PVPlacement( 0, ResinOffset
                                                      , logicSiPMResin
                                                      , "SiPMResin"
                                                      , logicWorld
                                                      , false
                                                      , 0
                                                      , checkOverlaps  );

  G4VPhysicalVolume* physSiPM = new G4PVPlacement( 0, SiPMOffset
                                                 , logicSiPM
                                                 , "SiPM"
                                                 , logicWorld
                                                 , false
                                                 , 0
                                                 , checkOverlaps );

  G4LogicalSkinSurface* CaseSurface
    = new G4LogicalSkinSurface( "SiPMCaseSurface"
                              , logicSiPMCase, fIdealWhiteOpSurface );
  G4LogicalSkinSurface* StandSurface
    = new G4LogicalSkinSurface( "SiPMStandSurface"
                              , logicSiPMStand, fIdealWhiteOpSurface );
  
  G4LogicalSkinSurface* PCBSurface
    = new G4LogicalSkinSurface( "PCBSurface", logicPCB, fPCBSurface );
  */

  
  ///////////////////////////////////////////////////////////////////////////////
  // Simple version of SiPM
  ///////////////////////////////////////////////////////////////////////////////
  const G4ThreeVector SiPMOffset_chan3( _hole_x1, 0, _WLSfiberZ - _WLS_zoff + 0.8*_sipm_z);  
  const G4ThreeVector SiPMOffset_chan4( _hole_x1, 0, -_WLSfiberZ - _WLS_zoff - 0.8*_sipm_z);
  //G4Box* solidSiPMInnerBox = new G4Box( "solidSiPMInnerBox", 0.5*_sipm_x, 0.5*_sipm_y,  0.8*_sipm_z );
  G4Tubs* solidSiPMInnerBox = new G4Tubs( "solidSiPMInnerBox", 0., _WLSfiberR+_WLSfiber_clad_thick, 0.8*_sipm_z , 0., 2*pi);
  G4LogicalVolume* logicSiPM = new G4LogicalVolume( solidSiPMInnerBox
                                                  , fBialkali,  "SiPM" );
                                                  
  G4VPhysicalVolume* physSiPM_chan3 = new G4PVPlacement( 0, SiPMOffset_chan3
                                                 , logicSiPM
                                                 , "physSiPM_chan3"
                                                 , logicWorld
                                                 , false
                                                 , 0
                                                 , checkOverlaps );
  G4VPhysicalVolume* physSiPM_chan4 = new G4PVPlacement( 0, SiPMOffset_chan4
                                                 , logicSiPM
                                                 , "physSiPM_chan4"
                                                 , logicWorld
                                                 , false
                                                 , 0
                                                 , checkOverlaps );

  G4LogicalSkinSurface* SiPMSurface
    = new G4LogicalSkinSurface( "SiPMSurface", logicSiPM, fSiPMSurface );
  ///////////////////////////////////////////////////////////////////////////////
  // Defining surfaces
  ///////////////////////////////////////////////////////////////////////////////
  
  // Tile surfaces
  G4LogicalBorderSurface* TileBulkSurface =
    new G4LogicalBorderSurface( "TileBulkSurface"
                              , physTileBulk
                              , physWorld
                              , fTileBulkSurface );
  std::cout << fTileBulkSurface << std::endl;

  // Other optical surfaces
  
  // Setting the sensitive detector
  if( !fPMTSD ){
    fPMTSD = new LYSimPMTSD( "/LYSimPMT" );
    G4SDManager* sdman = G4SDManager::GetSDMpointer();
    sdman->AddNewDetector( fPMTSD );
  }
  logicSiPM->SetSensitiveDetector( fPMTSD );

  // Visual attributes
  logicWorld->SetVisAttributes( G4VisAttributes::Invisible );

  // Avoid Colours Green and Yellow, since these are defaulted to the optical
  // Photons (I don't know how to change this)
  G4VisAttributes* SiPMVisAtt = new G4VisAttributes( G4Colour( 0, 0, 0 ) );
  SiPMVisAtt->SetForceSolid( true );
  SiPMVisAtt->SetVisibility( true );
  logicSiPM->SetVisAttributes( SiPMVisAtt );

  G4VisAttributes* CaseVisAtt = new G4VisAttributes( G4Colour( 0.8, 0.8, 0.8 ) );
  CaseVisAtt->SetForceSolid( true );
  CaseVisAtt->SetVisibility( true );
  //logicSiPMCase->SetVisAttributes( CaseVisAtt );
  //logicSiPMStand->SetVisAttributes( CaseVisAtt );

  G4VisAttributes* ResinVisAtt = new G4VisAttributes( G4Colour( 0., 1., 1. ) );
  ResinVisAtt->SetForceWireframe( true );
  ResinVisAtt->SetVisibility( true );
  //logicSiPMResin->SetVisAttributes( ResinVisAtt );

  G4VisAttributes* TileVisAtt = new G4VisAttributes( G4Colour( 0, 0, 1. ) );
  TileVisAtt->SetForceWireframe( true );
  TileVisAtt->SetVisibility( true );
  logicTileBulk->SetVisAttributes( TileVisAtt );

  G4VisAttributes* fiberVisAtt = new G4VisAttributes( G4Colour( 0.5, 0.5, 0 ) );
  fiberVisAtt->SetForceWireframe( true );
  fiberVisAtt->SetVisibility( true );
  fiberVisAtt->SetForceAuxEdgeVisible( true );
  //fiberVisAtt->SetForceLineSegmentsPerCircle( 30 );
  logicWLSfiber->SetVisAttributes( fiberVisAtt );
  logicWLSfiber_clad->SetVisAttributes( fiberVisAtt );

  G4VisAttributes* WrapVisAtt = new G4VisAttributes( G4Colour( 0.5, 0.5, 1.0 ) );
  WrapVisAtt->SetForceWireframe( true );
  WrapVisAtt->SetVisibility( true );
  logicWrap->SetVisAttributes( WrapVisAtt );


  G4VisAttributes* PCBVisAtt = new G4VisAttributes( G4Colour( 0.0, 0.4, 0.1 ) );
  PCBVisAtt->SetForceSolid( true );
  PCBVisAtt->SetVisibility( true );
  //logicPCB->SetVisAttributes( PCBVisAtt );

  return physWorld;
}

G4VSolid*
LYSimDetectorConstruction::ConstructTrapazoidSolid(
  const G4String& name,
  double          x,
  double          y,
  double          z,
  double          indent_x1,
  double          indent_x2 ) const
{
  const G4ThreeVector corners[8] = {
    G4ThreeVector( -x/2,           -y/2, -z/2 ),
    G4ThreeVector( x/2,            -y/2, -z/2 ),
    G4ThreeVector( -x/2+indent_x2, y/2,  -z/2 ),
    G4ThreeVector( x/2-indent_x1,  y/2,  -z/2 ),
    G4ThreeVector( -x/2,           -y/2, z/2 ),
    G4ThreeVector( x/2,            -y/2, z/2 ),
    G4ThreeVector( -x/2+indent_x2, y/2,  z/2 ),
    G4ThreeVector( x/2-indent_x1,  y/2,  z/2 )
  };

  return new G4Trap( name, corners );
}

G4VSolid*
LYSimDetectorConstruction::ConstructHollowWrapSolid() const
{
  G4VSolid* wrapOuter
    = ConstructTrapazoidSolid( "WrapOuter"
                             , _tilex + 2*wrapgap + 2*wrapthickness
                             , _tiley + 2*wrapgap + 2*wrapthickness
                             , _tilez + 2*wrapgap + 2*wrapthickness
                             , 0, 0 );
  G4VSolid* wrapInner
    = ConstructTrapazoidSolid( "WrapInner"
                             , _tilex + 2*wrapgap
                             , _tiley + 2*wrapgap
                             , _tilez + 2*wrapgap
                             , 0, 0 );
  G4VSolid* wrapbox = new G4SubtractionSolid( "WrapBox"
                                            , wrapOuter, wrapInner );

  return wrapbox;
}

G4VSolid*
LYSimDetectorConstruction::ConstructHollowWrapCladSolid() const
{
  G4VSolid* wrapOuter
    = ConstructTrapazoidSolid( "WrapOuter"
                             , _tilex  + 2*wrapthickness
                             , _tiley  + 2*wrapthickness
                             , _tilez
                             , 0, 0 );
  G4VSolid* wrapInner
    = ConstructTrapazoidSolid( "WrapInner"
                             , _tilex 
                             , _tiley
                             , _tilez*1.5
                             , 0, 0 );
  G4VSolid* wrapbox = new G4SubtractionSolid( "WrapBox"
                                            , wrapOuter, wrapInner );
  return wrapbox;
}

G4VSolid*
LYSimDetectorConstruction::ConstructDimpleSubtract() const
{
  switch( _dimple_type ){
  case SPHERICAL:
    return ConstructSphereDimpleSolid();
  case FLAT_DOME:
    return ConstructFlatDomeDimpleSolid();
  case CYLINDRICAL:
    return ConstructCylindricalDimpleSolid();
  case ELLIPSOID:
    return ConstructEllipsoidDimpleSolid();
  default:
    return ConstructSphereDimpleSolid();
  }
}

G4VSolid*
LYSimDetectorConstruction::ConstructSphereDimpleSolid() const
{
  return new G4Orb( "DimpleSolid"
                  , GetDimpleSphereRadius() );
}

G4VSolid*
LYSimDetectorConstruction::ConstructFlatDomeDimpleSolid() const
{
  G4VSolid* torus = new G4Torus( "DimpleTorus"
                               , 0, _dimple_indent
                               , _dimple_radius - _dimple_indent
                               , 0, 2*pi );
  G4VSolid* tub = new G4Tubs( "DimpleTub"
                            , 0, _dimple_radius-_dimple_indent
                            , _dimple_indent
                            , 0, 2*pi );

  return new G4UnionSolid( "DimpleSolid", torus, tub );
}

G4VSolid*
LYSimDetectorConstruction::ConstructCylindricalDimpleSolid() const
{
  return new G4Tubs( "DimpleSolid"
                   , 0.0, _dimple_radius
                   , _dimple_indent
                   , 0, 2*pi  );
}

G4VSolid*
LYSimDetectorConstruction::ConstructEllipsoidDimpleSolid() const
{
  return new G4Ellipsoid( "DimpleSolid"
                        , _dimple_radius
                        , _dimple_radius
                        , _dimple_indent );
}

G4ThreeVector
LYSimDetectorConstruction::CalcDimpleOffset() const
{
  switch( _dimple_type ){
  case FLAT_DOME:
    return G4ThreeVector( 0, 0, 0.5* _tilez );
  case CYLINDRICAL:
    return G4ThreeVector( 0, 0, 0.5*_tilez );
  case ELLIPSOID:
    return G4ThreeVector( 0, 0, 0.5*_tilez );
  default:
    return G4ThreeVector( 0, 0,
      0.5*_tilez + GetDimpleSphereRadius() - _dimple_indent );
  }
}

// Additional geometry factors
double
LYSimDetectorConstruction::WorldHalfX() const
{
  return 1*m;
}

double
LYSimDetectorConstruction::WorldHalfY() const
{
  return 1*m;
}

double
LYSimDetectorConstruction::WorldHalfZ() const
{
  return 10*m;
}

double
LYSimDetectorConstruction::LocalTileZ( const double x, const double y ) const
{
  /*
  const double r = sqrt( x*x+y*y );
  if( r < _dimple_radius ){
    switch( _dimple_type ){
    case FLAT_DOME:
      return LocalTileZFlatDome( r );
    case CYLINDRICAL:
      return LocalTileZCylindrical( r );
    case ELLIPSOID:
      return LocalTileZEllipsoid( r );
    default:
      return LocalTileZSpherical( r );
    }
  } else {
    return _tilez;
  }
  */
  if(x<_tile_x1+_hole_radius && x>_tile_x1-_hole_radius){
     return _tilez - 2*sqrt(pow(_hole_radius,2)-pow(x-_tile_x1,2));
  } else if(x<_tile_x2+_hole_radius && x>_tile_x2-_hole_radius){
    return _tilez - 2*sqrt(pow(_hole_radius,2)-pow(x-_tile_x2,2));
  }else{
    return _tilez;
  }
}


double
LYSimDetectorConstruction::LocalTileZFlatDome( const double r ) const
{
  if( r < _dimple_radius - _dimple_indent ){
    return _tilez - _dimple_indent;
  } else {
    const double dif = r- _dimple_radius + _dimple_indent;
    return _tilez - sqrt( _dimple_indent * _dimple_indent - dif * dif );
  }
}

double
LYSimDetectorConstruction::LocalTileZCylindrical( const double r ) const
{
  return _tilez - _dimple_indent;
}

double
LYSimDetectorConstruction::LocalTileZEllipsoid( const double r ) const
{
  return _tilez - _dimple_indent * sqrt(
    1 - r * r / _dimple_radius / _dimple_radius
    );
}

double
LYSimDetectorConstruction::LocalTileZSpherical( const double r ) const
{
  const double big_r = GetDimpleSphereRadius();
  return _tilez  - ( sqrt( big_r*big_r - r*r ) - ( big_r - _dimple_indent ) );
}


double
LYSimDetectorConstruction::GetDimpleSphereRadius() const
{
  return 0.5 * (
    ( _dimple_radius*_dimple_radius )/_dimple_indent
    + _dimple_indent );
}




void
LYSimDetectorConstruction::SetTileAbsMult( const double mult )
{
  _absmult = mult;
  Update_EJ200_AbsLength( fEJ200, _absmult );
}

void
LYSimDetectorConstruction::SetTileScintillation( const double mult )
{
  _ScintiN = mult;
  Update_EJ200_Scinti( fEJ200, _ScintiN );
}
void
LYSimDetectorConstruction::SetWrapReflect( const double r )
{
  // Add entries into properties table
  _wrap_reflect = r;
  static const unsigned nentries = 2;
  static double phoE[nentries]   = {1.0*eV, 6.0*eV};
  double reflectivity[nentries]  = {_wrap_reflect, _wrap_reflect};

  G4MaterialPropertiesTable* table = fESROpSurface->GetMaterialPropertiesTable();
  if( table ){
    table->RemoveProperty( "REFLECTIVITY" );
    table->AddProperty( "REFLECTIVITY", phoE, reflectivity, nentries );
  } else {
    table = new G4MaterialPropertiesTable();
    table->AddProperty( "REFLECTIVITY", phoE, reflectivity, nentries );
    fESROpSurface->SetMaterialPropertiesTable( table );
  }
  G4MaterialPropertiesTable* table2 = fTiO2Surface->GetMaterialPropertiesTable();
  if( table2 ){
    table2->RemoveProperty( "REFLECTIVITY" );
    table2->AddProperty( "REFLECTIVITY", phoE, reflectivity, nentries );
  } else {
    table2 = new G4MaterialPropertiesTable();
    table2->AddProperty( "REFLECTIVITY", phoE, reflectivity, nentries );
    fTiO2Surface->SetMaterialPropertiesTable( table );
  }
  
}

void
LYSimDetectorConstruction::SetPCBReflect( const double r )
{
  // Add entries into properties table
  _pcb_reflectivity = r;
  static const unsigned nentries = 2;
  static double phoE[nentries]   = {1.0*eV, 6.0*eV};
  double reflectivity[nentries]  = {r, r};

  G4MaterialPropertiesTable* table = fPCBSurface->GetMaterialPropertiesTable();
  if( table ){
    table->RemoveProperty( "REFLECTIVITY" );
    table->AddProperty( "REFLECTIVITY", phoE, reflectivity, nentries );
  } else {
    table = new G4MaterialPropertiesTable();
    table->AddProperty( "REFLECTIVITY", phoE, reflectivity, nentries );
    fPCBSurface->SetMaterialPropertiesTable( table );
  }
}


void
LYSimDetectorConstruction::SetTileAlpha( const double a )
{
  fTileBulkSurface->SetSigmaAlpha( a );
  _tile_alpha = a;
}

void
LYSimDetectorConstruction::SetDimpleAlpha( const double a )
{
  fTileDimpleSurface->SetSigmaAlpha( a );
  _dimple_alpha = a;
}


