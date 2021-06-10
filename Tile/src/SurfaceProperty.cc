#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"
#include "G4OpticalSurface.hh"

using namespace CLHEP;

static const unsigned flatentries = 2;
static const double minenergy = 1.0 * eV;
static const double maxenergy = 8.0 * eV;

G4OpticalSurface *
MakeS_dimple()
{
  static const unsigned nentries = 2;
  static double phoE[nentries] = {1.0 * eV, 6.0 * eV};
  static double specularlobe[nentries] = {1.0, 1.0};
  G4OpticalSurface *surface = new G4OpticalSurface("Dimple_Surface");
  surface->SetType(dielectric_dielectric);
  surface->SetModel(unified);
  surface->SetFinish(ground);
  surface->SetSigmaAlpha(1.3 * degree); // Janecek2010 (1.3 * degree)

  G4MaterialPropertiesTable *table = new G4MaterialPropertiesTable();

  table->AddProperty(
      "SPECULARLOBECONSTANT", phoE, specularlobe, nentries);
  surface->SetMaterialPropertiesTable(table);
  return surface;
}

void Update_Dimple_SigmaAlpha(G4OpticalSurface *surface, const double mult)
{
  std::cout << "------------------------------------" << std::endl;
  surface->DumpInfo();
  std::cout << "the original SA: " << surface->GetSigmaAlpha() << std::endl;
  surface->SetSigmaAlpha(mult * degree);
  std::cout << "Set Sigma Alpha to: " << mult << " degree" << std::endl;
}

G4OpticalSurface *
MakeS_TyvekCrystal()
{
  const unsigned num = 2;
  double Ephoton[num] = {1.5 * eV, 8.0 * eV};
  double Reflectivity[num] = {0.979, 0.979};

  //////////////////////////////////
  // Realistic Crystal-Tyvek surface
  //////////////////////////////////
  G4OpticalSurface *surface = new G4OpticalSurface("TyvekOpSurface");
  surface->SetType(dielectric_LUT);
  surface->SetModel(LUT);
  surface->SetFinish(polishedtyvekair);

  G4MaterialPropertiesTable *table = new G4MaterialPropertiesTable();

  table->AddProperty("REFLECTIVITY", Ephoton, Reflectivity, num);

  surface->SetMaterialPropertiesTable(table);

  return surface;
}

G4OpticalSurface *
MakeS_ESR()
{
  const unsigned num = 2;
  double Ephoton[num] = {1.0 * eV, 8.0 * eV};
  double Reflectivity[num] = {0.985, 0.985};
  // source: https://www.osti.gov/servlets/purl/1184400
  //////////////////////////////////
  // ESR surface
  //////////////////////////////////
  G4OpticalSurface *surface = new G4OpticalSurface("ESROpSurface");
  /*  surface->SetType( dielectric_LUT );
  surface->SetModel( LUT );
  surface->SetFinish( polishedvm2000air );
*/

  surface->SetType(dielectric_dielectric);
  surface->SetModel(unified);
  surface->SetFinish(polished);
  //  surface->SetSigmaAlpha( 1.3*degree );

  G4MaterialPropertiesTable *table = new G4MaterialPropertiesTable();
  table->AddProperty("REFLECTIVITY", Ephoton, Reflectivity, num);
  surface->SetMaterialPropertiesTable(table);

  return surface;
}

G4OpticalSurface *
MakeS_IdealTyvekCrystal()
{
  //////////////////////////////////
  // Ideal Crystal-Tyvek surface
  //////////////////////////////////
  G4OpticalSurface *surface = new G4OpticalSurface("IdealTyvekOpSurface");
  surface->SetType(dielectric_LUT);
  surface->SetModel(LUT);
  surface->SetFinish(polishedtyvekair);

  return surface;
}

G4OpticalSurface *
MakeS_UnifiedTyvek()
{
  const unsigned nentries = 2;
  double phoE[nentries] = {2.0 * eV, 6.0 * eV};
  double specularlobe[nentries] = {1.0, 1.0};
  double specularspike[nentries] = {0., 0.};
  double backscatter[nentries] = {0., 0.};
  double rindex[nentries] = {1.0, 1.0};
  double reflectivity[nentries] = {0.979, 0.979};
  double efficiency[nentries] = {0.0, 0.0};

  //////////////////////////////////
  // Unified Tyvek surface
  //////////////////////////////////
  G4OpticalSurface *surface = new G4OpticalSurface("UnifiedTyvekOpSurface");
  surface->SetType(dielectric_dielectric);
  surface->SetModel(unified);
  surface->SetFinish(groundbackpainted);
  surface->SetSigmaAlpha(1.3 * degree);

  G4MaterialPropertiesTable *table = new G4MaterialPropertiesTable();

  table->AddProperty("RINDEX", phoE, rindex, nentries);
  table->AddProperty("SPECULARLOBECONSTANT", phoE, specularlobe, nentries);
  table->AddProperty("SPECULARSPIKECONSTANT", phoE, specularspike, nentries);
  table->AddProperty("BACKSCATTERCONSTANT", phoE, backscatter, nentries);
  table->AddProperty("REFLECTIVITY", phoE, reflectivity, nentries);
  table->AddProperty("EFFICIENCY", phoE, efficiency, nentries);

  surface->SetMaterialPropertiesTable(table);

  return surface;
}

G4OpticalSurface *
MakeS_IdealTyvek()
{
  const unsigned nentries = 2;
  double phoE[nentries] = {2.0 * eV, 6.0 * eV};
  double specularlobe[nentries] = {1.0, 1.0};
  double specularspike[nentries] = {0., 0.};
  double backscatter[nentries] = {0., 0.};
  double rindex[nentries] = {1.0, 1.0};
  double reflectivity[nentries] = {1.0, 1.0};
  double efficiency[nentries] = {0.0, 0.0};

  //////////////////////////////////
  // Unified Ideal Tyvek surface
  //////////////////////////////////
  G4OpticalSurface *surface = new G4OpticalSurface(
      "UnifiedIdealTyvekOpSurface");
  surface->SetType(dielectric_dielectric);
  surface->SetModel(unified);
  surface->SetFinish(groundbackpainted);
  surface->SetSigmaAlpha(1.3 * degree);

  G4MaterialPropertiesTable *table = new G4MaterialPropertiesTable();

  table->AddProperty("RINDEX", phoE, rindex, nentries);
  table->AddProperty("SPECULARLOBECONSTANT", phoE, specularlobe, nentries);
  table->AddProperty("SPECULARSPIKECONSTANT", phoE, specularspike, nentries);
  table->AddProperty("BACKSCATTERCONSTANT", phoE, backscatter, nentries);
  table->AddProperty("REFLECTIVITY", phoE, reflectivity, nentries);
  table->AddProperty("EFFICIENCY", phoE, efficiency, nentries);

  surface->SetMaterialPropertiesTable(table);

  return surface;
}

G4OpticalSurface *
MakeS_Polished()
{
  static const unsigned nentries = 2;
  static double phoE[nentries] = {1.0 * eV, 6.0 * eV};
  static double specularlobe[nentries] = {1.0, 1.0};

  //////////////////////////////////
  // Realistic polished surface
  //////////////////////////////////
  G4OpticalSurface *surface = new G4OpticalSurface("PolishedOpSurface");
  surface->SetType(dielectric_dielectric);
  surface->SetModel(unified);
  surface->SetFinish(ground);
  // necessary even for polished surfaces to enable UNIFIED code
  surface->SetSigmaAlpha(1.3 * degree); // Janecek2010 (1.3 * degree)

  G4MaterialPropertiesTable *table = new G4MaterialPropertiesTable();

  table->AddProperty(
      "SPECULARLOBECONSTANT", phoE, specularlobe, nentries);
  surface->SetMaterialPropertiesTable(table);

  return surface;
}

G4OpticalSurface *
MakeS_IdealPolished()
{
  //////////////////////////////////
  // Ideal polished surface
  //////////////////////////////////
  G4OpticalSurface *surface = new G4OpticalSurface("IdealOpSurface");
  surface->SetType(dielectric_dielectric);
  surface->SetModel(glisur);
  surface->SetFinish(polished);

  return surface;
}

G4OpticalSurface *
MakeS_Mirror()
{
  const unsigned nentries = 2;
  double phoE[nentries] = {1.0 * eV, 6.0 * eV};
  double reflectivity[nentries] = {1.0, 1.0};
  //////////////////////////////////
  // Mirror surface
  //////////////////////////////////
  G4OpticalSurface *surface = new G4OpticalSurface("MirrorOpSurface");
  surface->SetType(dielectric_metal);
  surface->SetFinish(polished);
  surface->SetModel(unified);

  G4MaterialPropertiesTable *table = new G4MaterialPropertiesTable();
  table->AddProperty("REFLECTIVITY", phoE, reflectivity, nentries);
  surface->SetMaterialPropertiesTable(table);

  return surface;
}

G4OpticalSurface *
MakeS_IdealMirror()
{
  const unsigned nentries = 2;
  double phoE[nentries] = {1.0 * eV, 6.0 * eV};
  double reflectivity[nentries] = {1.0, 1.0};
  //////////////////////////////////
  // Ideal mirror surface
  //////////////////////////////////
  G4OpticalSurface *surface = new G4OpticalSurface("MirrorOpSurface");
  surface->SetType(dielectric_metal);
  surface->SetFinish(polished);
  surface->SetModel(unified);

  G4MaterialPropertiesTable *table = new G4MaterialPropertiesTable();
  table->AddProperty("REFLECTIVITY", phoE, reflectivity, nentries);
  surface->SetMaterialPropertiesTable(table);

  return surface;
}
G4OpticalSurface *
MakeS_DarkSurface()
{
  G4OpticalSurface *surface = new G4OpticalSurface("DarkOpSurface");
  surface->SetType(dielectric_metal);
  surface->SetFinish(ground);
  surface->SetModel(unified);

  double phoE[flatentries] = {minenergy, maxenergy};
  double reflectivity[flatentries] = {0.10, 0.10};

  G4MaterialPropertiesTable *table = new G4MaterialPropertiesTable();
  table->AddProperty("REFLECTIVITY", phoE, reflectivity, flatentries);
  surface->SetMaterialPropertiesTable(table);

  return surface;
}
G4OpticalSurface *
MakeS_IdealWhiteSurface()
{
  //////////////////////////////////
  // the case surface
  //////////////////////////////////
  G4OpticalSurface *surface = new G4OpticalSurface("WhiteOpSurface");
  surface->SetType(dielectric_metal);
  surface->SetFinish(ground);
  surface->SetModel(unified);

  double phoE[flatentries] = {minenergy, maxenergy};
  double reflectivity[flatentries] = {0.50, 0.50};

  G4MaterialPropertiesTable *table = new G4MaterialPropertiesTable();
  table->AddProperty("REFLECTIVITY", phoE, reflectivity, flatentries);
  surface->SetMaterialPropertiesTable(table);

  return surface;
}

G4OpticalSurface *
MakeS_Absorbing()
{
  const unsigned nentries = 2;
  double phoE[nentries] = {1.0 * eV, 8.0 * eV};
  double reflectivity[nentries] = {1, 1};
  //////////////////////////////////
  // Absorbing surface
  //////////////////////////////////
  G4OpticalSurface *surface = new G4OpticalSurface("AbsorbingOpSurface");
  surface->SetType(dielectric_dielectric);
  surface->SetFinish(polished);
  //  surface->SetFinish( groundfrontpainted );
  surface->SetModel(unified);

  G4MaterialPropertiesTable *table = new G4MaterialPropertiesTable();
  table->AddProperty("REFLECTIVITY", phoE, reflectivity, nentries);
  surface->SetMaterialPropertiesTable(table);

  return surface;
}

G4OpticalSurface *
MakeS_Rough()
{
  //////////////////////////////////
  // Rough surface used for esr
  //////////////////////////////////
  G4OpticalSurface *surface = new G4OpticalSurface("RoughOpSurface");
  surface->SetType(dielectric_metal);
  surface->SetFinish(ground);
  surface->SetModel(unified);

  return surface;
}

G4OpticalSurface *
MakeS_SiPM_S14160()
{
  //1.3mm SiPM, S14160-1315PS
  const unsigned nentries = 51;
  double phoE[nentries] = {4.456977245 * eV,
                           4.405663883 * eV, 4.371335649 * eV, 4.337993527 * eV, 4.305156186 * eV, 4.272812248 * eV,
                           4.120034327 * eV, 4.047934735 * eV, 3.885068561 * eV, 3.795976762 * eV, 3.710990512 * eV,
                           3.641025285 * eV, 3.562661791 * eV, 3.416389546 * eV, 3.254178294 * eV, 3.200252775 * eV,
                           3.130913965 * eV, 3.072870849 * eV, 3.016867241 * eV, 2.948003733 * eV, 2.868012792 * eV,
                           2.772393126 * eV, 2.622893865 * eV, 2.488842802 * eV, 2.467789116 * eV, 2.382432947 * eV,
                           2.325721122 * eV, 2.276066914 * eV, 2.224250888 * eV, 2.174741594 * eV, 2.115879533 * eV,
                           2.067333517 * eV, 2.024496146 * eV, 1.98339801 * eV, 1.943935293 * eV, 1.906041585 * eV,
                           1.869568783 * eV, 1.834465614 * eV, 1.80065635 * eV, 1.76540215 * eV, 1.731477711 * eV,
                           1.696391875 * eV, 1.66504429 * eV, 1.628029216 * eV, 1.594775069 * eV, 1.562852229 * eV,
                           1.546342471 * eV, 1.512389673 * eV, 1.477990547 * eV, 1.443354983 * eV, 1.40859115 * eV};

  double efficiency[nentries] = {1.1045,
                                 2.4287, 4.4945, 5.8187, 7.1429, 8.4672, 11.9633, 13.2347, 15.5127, 16.7841, 18.0556,
                                 19.327, 20.5984, 21.6053, 22.5064, 23.7778, 25.7379, 27.0093, 28.2808, 29.4993, 30.7179,
                                 31.7778, 32.0435, 31.7796, 31.6738, 30.7211, 29.5563, 28.3914, 27.1736, 26.0618, 24.6852,
                                 23.5734, 22.4086, 21.2438, 20.079, 18.9142, 17.7494, 16.5846, 15.4198, 14.255, 13.1432,
                                 12.0315, 10.9197, 9.755, 8.6962, 7.6375, 7.1611, 6.1554, 5.3086, 4.4618, 3.6681};
  for (int dd = 0; dd < 51; dd++)
    efficiency[dd] /= 100.;
  const unsigned ref_ent = 32;
  double phoE2[ref_ent] = {
      4.928613174 * eV, 4.760932071 * eV, 4.57675131 * eV, 4.416020551 * eV, 4.27782469 * eV, 4.132944198 * eV, 3.988297134 * eV,
      3.870877084 * eV, 3.751866883 * eV, 3.647774074 * eV, 3.538260695 * eV, 3.441997529 * eV, 3.347666946 * eV, 3.264460058 * eV,
      3.176475533 * eV, 3.098752668 * eV, 3.019365195 * eV, 2.751902006 * eV, 2.432636667 * eV, 2.254544997 * eV, 2.136148464 * eV,
      2.032627719 * eV, 1.907273067 * eV, 1.796976535 * eV, 1.722050515 * eV, 1.570692624 * eV, 1.511965476 * eV,
      1.459049532 * eV, 1.377862407 * eV, 1.305495288 * eV, 1.265726027 * eV, 1.240350474 * eV};
  double reflectivity[ref_ent] = {
      0.67091, 0.70309, 0.73346, 0.72132, 0.68334, 0.62105, 0.59036,
      0.57232, 0.5645, 0.56008, 0.56441, 0.58236, 0.5828, 0.54482, 0.51025,
      0.48443, 0.46541, 0.41855, 0.38169, 0.3657, 0.35705, 0.35096, 0.34397,
      0.33772, 0.33321, 0.32645, 0.32498, 0.32241, 0.31988, 0.31625, 0.31417, 0.31323};

  G4OpticalSurface *surface = new G4OpticalSurface("SiPM_Surface_S14");
  surface->SetType(dielectric_metal);
  surface->SetModel(unified);
    double flat_phoEs13[flatentries]   = {1*eV, 8*eV };
    double flat_reflectivitys13[flatentries] = {0, 0};

  G4MaterialPropertiesTable *table = new G4MaterialPropertiesTable();
  table->AddProperty("EFFICIENCY", phoE, efficiency, nentries);
//  table->AddProperty("REFLECTIVITY", phoE2, reflectivity, ref_ent);
  table->AddProperty("REFLECTIVITY", flat_phoEs13, flat_reflectivitys13, 2);
  surface->SetMaterialPropertiesTable(table);
  return surface;

}

G4OpticalSurface *
MakeS_SiPM_S13360()
{
  // From Meeting on HB SiPM selection (James Hirschauer)

  //1.3mm SiPM, S13360-1350PE
  const unsigned nSiPM13 = 63;
  double phoE_nSiPM13[nSiPM13] = {4.483036162 * eV, 4.419483603 * eV, 4.388386036 * eV,
                                  4.372985366 * eV, 4.357707729 * eV, 4.34253647 * eV, 4.31249367 * eV, 4.268208226 * eV, 4.253652707 * eV,
                                  4.210561468 * eV, 4.182311055 * eV, 4.11333664 * eV, 4.020487417 * eV, 3.94418248 * eV, 3.835017337 * eV,
                                  3.754184421 * eV, 3.687569314 * eV, 3.591974766 * eV, 3.530945075 * eV, 3.471954596 * eV, 3.414902801 * eV,
                                  3.341684828 * eV, 3.288801344 * eV, 3.229183902 * eV, 3.171689396 * eV, 3.116206415 * eV, 3.055122245 * eV,
                                  2.989205975 * eV, 2.912409692 * eV, 2.826584982 * eV, 2.763926011 * eV, 2.669260047 * eV, 2.586220633 * eV,
                                  2.523414636 * eV, 2.458730558 * eV, 2.383532172 * eV, 2.334444087 * eV, 2.291543551 * eV, 2.234062133 * eV,
                                  2.183212355 * eV, 2.116511629 * eV, 2.060547827 * eV, 2.010709904 * eV, 1.978798509 * eV, 1.938801027 * eV,
                                  1.891738959 * eV, 1.849646255 * eV, 1.814653335 * eV, 1.778420604 * eV, 1.736320087 * eV, 1.700797048 * eV,
                                  1.668926192 * eV, 1.638229971 * eV, 1.602438512 * eV, 1.57810373 * eV, 1.55643477 * eV, 1.529704014 * eV,
                                  1.493069504 * eV, 1.463283922 * eV, 1.428089882 * eV, 1.394549014 * eV, 1.376082338 * eV, 1.374564357 * eV};

  double efficiency_nSiPM13[nSiPM13] = {0, 0.030459, 0.044724,
                                        0.06395, 0.089997, 0.108603, 0.124729, 0.147677, 0.167523,
                                        0.179929, 0.195434, 0.214043, 0.227693, 0.241962, 0.259954,
                                        0.271123, 0.27919, 0.28664, 0.286646, 0.286651, 0.293478,
                                        0.300307, 0.315196, 0.328225, 0.346216, 0.358625, 0.369175,
                                        0.381586, 0.389657, 0.39587, 0.40022, 0.402715, 0.401488,
                                        0.393436, 0.384766, 0.368656, 0.355642, 0.342008, 0.325896,
                                        0.306683, 0.285613, 0.265162, 0.250912, 0.236037, 0.217444,
                                        0.201954, 0.184603, 0.16911, 0.154239, 0.141231, 0.128841,
                                        0.122652, 0.111501, 0.098493, 0.091061, 0.083009, 0.07558,
                                        0.064434, 0.055767, 0.045863, 0.0372, 0.03473, 0};


  // SiPM reflectivity paper: Using a flat 0.12 for now
  // https://reader.elsevier.com/reader/sd/pii/S016890021831297X?token=FADE142C43F500DEE5289D65B9D6573497F4AD220F3462EBA9239256CFD8E11FF98FF84E302016EAB89FA1B62B1EEDCE
    double flat_phoE2[flatentries]   = {1*eV, 8*eV };
    double flat_reflectivity[flatentries] = {0, 0};

  //reflectivity
  const unsigned ref_ent = 32;
  double phoE2[ref_ent] = {
      4.928613174 * eV, 4.760932071 * eV, 4.57675131 * eV, 4.416020551 * eV, 4.27782469 * eV, 4.132944198 * eV, 3.988297134 * eV,
      3.870877084 * eV, 3.751866883 * eV, 3.647774074 * eV, 3.538260695 * eV, 3.441997529 * eV, 3.347666946 * eV, 3.264460058 * eV,
      3.176475533 * eV, 3.098752668 * eV, 3.019365195 * eV, 2.751902006 * eV, 2.432636667 * eV, 2.254544997 * eV, 2.136148464 * eV,
      2.032627719 * eV, 1.907273067 * eV, 1.796976535 * eV, 1.722050515 * eV, 1.570692624 * eV, 1.511965476 * eV,
      1.459049532 * eV, 1.377862407 * eV, 1.305495288 * eV, 1.265726027 * eV, 1.240350474 * eV};
  double reflectivity[ref_ent] = {
      0.67091, 0.70309, 0.73346, 0.72132, 0.68334, 0.62105, 0.59036,
      0.57232, 0.5645, 0.56008, 0.56441, 0.58236, 0.5828, 0.54482, 0.51025,
      0.48443, 0.46541, 0.41855, 0.38169, 0.3657, 0.35705, 0.35096, 0.34397,
      0.33772, 0.33321, 0.32645, 0.32498, 0.32241, 0.31988, 0.31625, 0.31417, 0.31323};

  G4OpticalSurface *surface = new G4OpticalSurface("SiPM_Surface_S13");
  surface->SetType(dielectric_metal);
  surface->SetModel(unified);

  G4MaterialPropertiesTable *table = new G4MaterialPropertiesTable();


  table->AddProperty( "EFFICIENCY", phoE_nSiPM13, efficiency_nSiPM13,nSiPM13);
//  table->AddProperty("REFLECTIVITY", phoE2, reflectivity, ref_ent);
  table->AddProperty("REFLECTIVITY", flat_phoE2, flat_reflectivity, 2);
  surface->SetMaterialPropertiesTable(table);
  return surface;
}

G4OpticalSurface *
MakeS_PCBSurface()
{
  // PCB is a flat gray surface for now
  const int nentries = 43;
  double phoE2[nentries] = {360, 370, 380,
                            390, 400, 410, 420, 430, 440, 450, 460, 470, 480,
                            490, 500, 510, 520, 530, 540, 550, 560, 570, 580,
                            590, 600, 610, 620, 630, 640, 650, 660, 670, 680,
                            690, 700, 710, 720, 730, 740, 750, 760, 770, 780};
  double reflectivity[nentries] = {15.11, 15.4, 19.75,
                                   27.935, 38.715, 52.025, 59.16, 62.25, 63.19, 63.44, 63.325, 63.18, 62.705,
                                   62.295, 62.21, 62.035, 61.695, 61.26, 60.87, 60.505, 60.035, 59.485, 59.015,
                                   58.535, 57.89, 57.295, 57.07, 56.525, 56.09, 55.68, 55.07, 54.65, 54.18,
                                   53.72, 53.25, 52.83, 52.45, 52.195, 52.89, 54.16, 55.78, 57.59, 58.47};
  for (int i = 0; i < nentries; i++)
  {
    phoE2[i] = 1239.84193 * eV / phoE2[i];
    reflectivity[i] = reflectivity[i] / 100.0;
  }
  G4OpticalSurface *surface = new G4OpticalSurface("PCB_Surface");
  surface->SetType(dielectric_metal);
  surface->SetFinish(ground);
  surface->SetModel(unified);

  G4MaterialPropertiesTable *table = new G4MaterialPropertiesTable();
  table->AddProperty("REFLECTIVITY", phoE2, reflectivity, flatentries);
  surface->SetMaterialPropertiesTable(table);
  return surface;
}
