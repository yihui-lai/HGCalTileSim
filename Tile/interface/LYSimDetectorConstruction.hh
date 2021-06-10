#ifndef LYSimDetectorConstruction_h
#define LYSimDetectorConstruction_h 1

#ifdef CMSSW_GIT_HASH
#include "HGCalTileSim/Tile/interface/LYSimPMTSD.hh"
#else
#include "LYSimPMTSD.hh"
#endif

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class LYSimDetectorMessenger;
class G4Box;
class G4LogicalVolume;
class G4Material;
class G4OpticalSurface;
class LYSimPMTSD;

class LYSimDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  LYSimDetectorConstruction();
  virtual ~LYSimDetectorConstruction();
  virtual G4VPhysicalVolume* Construct();
  void                       UpdateGeometry();

  // Get/Set functions of tile geometry parameters
  inline void
  SetTileZ ( const G4double z ){ _tilez = z;}
  inline G4double
  GetTileZ() const {return _tilez + (tile_number-1)*(tilegap+_tilez);}
  inline G4double
  Gettgap() const {return tilegap;}
  inline G4int
  Gettn() const {return tile_number;}
  inline void
  SetTileX ( const G4double x ){_tilex = x;}
  inline G4double
  GetTileX() const {return _tilex;}
  inline void
  SetTileY ( const G4double y ){_tiley = y;}
  inline G4double
  GetTileY() const {return _tiley;}

  inline void
  SetTileN ( const G4int n ){tile_number = n;}
  inline G4int
  GetTileN() const {return tile_number;}
  inline void
  SetWrapESR ( const G4bool tof ){is_ESR = tof;}
  inline G4bool
  GetWrapESR() const {return is_ESR;}

  inline void
  SetMaterial ( const int m ){material = m;}
  inline int
  GetMaterial () const {return material;}

  inline void
  SetPDE ( const int ee ){sipm_pde = ee;}
  inline int
  GetPDE () const {return sipm_pde;}

  //borrowed for gap study
  inline void
  SetDimpleSA ( const double d ){_dimplesa = d;wrapgap = d;}
  inline double
  GetDimpleSA () const {return _dimplesa;}

  // Get/Set functions of dimple geometry parameters
  inline double
  GetDimpleRadius() const { return _dimple_radius; }
  inline void
  SetDimpleRadius( const double x ){ _dimple_radius = x; }
  inline double
  GetDimpleIndent() const { return _dimple_indent; }
  inline void
  SetDimpleIndent( const double x ){ _dimple_indent = x; }
  int
  GetDimpleType() const {return _dimple_type; }
  inline void
  SetDimpleType( const int x ){ _dimple_type = x; }

  // Get/Set functions for SiPM geometry parameters
  inline double
  GetSiPMX() const {return _sipm_x; }
  inline void
  SetSiPMX( const double x ){ _sipm_x = x; }
  inline double
  GetSiPMY() const { return _sipm_y; }
  inline void
  SetSiPMY( const double x ){ _sipm_y = x; }
  inline double
  GetSiPMRim() const {return _sipm_rimwidth; }
  inline void
  SetSiPMRim( const double x ){ _sipm_rimwidth = x; }
  inline double
  GetSiPMStand() const { return _sipm_standz; }
  inline void
  SetSiPMStand( const double x ){ _sipm_standz = x; }
  inline double
  Getgap_pcb_wrap() const { return _gap_pcb_wrap; }
  inline void
  Setgap_pcb_wrap( const double x ){ _gap_pcb_wrap = x; }

  enum  DIMPLE_TYPE
  {
    SPHERICAL,
    PARABOLIC,
    PYRAMID
  };

  // Additional geometric factors.
  double WorldHalfX() const;
  double WorldHalfY() const;
  double WorldHalfZ() const;

  double GetDimpleSizeRadius() const;
  double LocalTileZ( const double x, const double y ) const;

  // Material updating functions
  void SetTileAbsMult( const double x );
  void SetSurfaceDimpleSA( const double x );

  inline double
  GetTileAbsMult() const { return _absmult; }

  void SetWrapReflect( const double x );
  inline double
  GetWrapReflect() const { return _wrap_reflect; }

  void
  SetPCBRadius( const double x ){ _pcb_radius = x; }
  inline double
  GetPCBRadius() const { return _pcb_radius; }

  void SetPCBReflect( const double x );
  inline double
  GetPCBReflect() const { return _pcb_reflectivity; }

  // Accesing Custom Materials
  inline G4Material*
  GetEJ200(){ return fTileMaterial; }
  inline G4OpticalSurface*
  GetSiPMSurface(){ return fSiPMSurface; }

private:
  // Subfunctions for cleaner code
  G4VSolid* ConstructTrapazoidSolid(
    const G4String& name,
    double          x,// Maximum Length along z direction.
    double          y,// Length along y direction
    double          z,// length along z direction
    double          indent_x1,// Indent in x direction
    double          indent_x2// indent of corner
    ) const;

  G4VSolid* ConstructHollowWrapSolid() const;
  G4VSolid* ConstructUnEvenHollowWrapSolid() const;
  G4VSolid* ConstructHollowWrapGapSolid() const;
  G4VSolid* ConstructUnEvenHollowWrapGapSolid() const;
  G4VSolid* ConstructSphereDimpleSolid() const;
  G4VSolid* ConstructPyramidDimpleSolid() const;
  G4VSolid* ConstructParabolicDimpleSolid() const;
  G4VSolid* ConstructArbiDimpleSolid() const;

  G4ThreeVector CalcDimpleOffset() const;
  G4ThreeVector CalcSiPMDimpleOffset() const;

  // Pointer to detector messenger class
  LYSimDetectorMessenger* fdetectorMessenger;

  // Pointers to materials
  G4Material* fAir;
  G4Material* fTileMaterial;
  G4Material* scsn81;
  G4Material* fBialkali;
  G4Material* fEpoxy;
  G4Material* fResin;
  G4Material* fBC_630_grease;

  // Pointers to surfaces
  G4OpticalSurface* fTyvekSurface;
  G4OpticalSurface* fESROpSurface;
  G4OpticalSurface* fPolishedOpSurface;
  G4OpticalSurface* fIdealPolishedOpSurface;
  G4OpticalSurface* fIdealWhiteOpSurface;
  G4OpticalSurface* fSiPMSurface;
  G4OpticalSurface* fPCBSurface;
  G4OpticalSurface* fdimpleSurface;
  G4OpticalSurface* fGrease;
  G4OpticalSurface* fDarkOpSurface;

  // Pointers for access to Sensitive Detector
  static LYSimPMTSD* fPMTSD;

  // Tile Geometric definitions
  double _tilex;
  double _tiley;
  double _tilez;
  double _tile_x1;
  double _tile_x2;

  // Wrapping Geometry definitions
  double wrapgap;// distance between tile and wrapping
  double wrapthickness;

  double tilegap;
  int tile_number;
  bool is_ESR;
  int material;
  int sipm_pde;
  double _dimplesa;
  // Absorption length multiplier
  double _absmult;
  double _wrap_reflect;

  // Dimple variables;
  double _dimple_indent;// Depth of dimple
  double _dimple_radius;// Dimple radius
  int _dimple_type;

  // Parameters for plate behind SiPM
  double _pcb_radius;
  double _pcb_reflectivity;

  double _sipm_deadwidth;
  double _sipm_x;
  double _sipm_y;
  double _sipm_z;
  double _sipm_standz;
  double _gap_pcb_wrap;
  double _sipm_rimwidth;
  double _sipm_glasswidth;


};

#endif
