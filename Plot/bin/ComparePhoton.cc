#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/Common/interface/Format.hpp"
#include "UserUtils/Common/interface/STLUtils/Filesystem.hpp"
#include "UserUtils/Common/interface/STLUtils/VectorUtils.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

#include "HGCalTileSim/Plot/interface/LYvX_Common.hpp"

#include "TChain.h"
#include "TFile.h"
#include "TGraphErrors.h"

int
main( int argc, char** argv )
{
  usr::po::options_description desc( "Input options" );
  desc.add_options()
    ( "inputfile,f", usr::po::value<std::string>(),
    "File containing the generated files" )
    ( "output,o", usr::po::defvalue<std::string>(
    usr::resultpath( "HGCalTileSim", "Plot" ).string() ),
    "Output plot directory" )
  ;
  usr::ArgumentExtender args;
  args.AddOptions( GeometryOptions() );
  args.AddOptions( desc );
  args.ParseOptions( argc, argv );
  args.SetFilePrefix( args.Arg( "output" ) );
  args.AddNameScheme( {
    {"dimplerad", "r"}, {"dimpleind", "d"},
    {"tilewidth", "L"}, {"absmult", "a"}, {"wrapreflect", "m"},
    {"sipmwidth", "W"}, {"sipmstand", "S"}, {"PCBRadius", "b"}, {"PCBRef", "P"}
  } );

  const std::vector<std::string> options_list = {
    "dimplerad", "dimpleind",
    "tilewidth",
    "absmult",   "wrapreflect",
    "sipmwidth", "sipmstand",
    "PCBRadius", "PCBRef"
  };
  TChain tree( "LYSimRun", "LYSimRun" );
  LYSimRunFormat fmt;
  LYvXGraphContainer cont;
  int entry = -1;

  tree.Add( args.Arg<std::string>( "inputfile" ).c_str() );
  fmt.LoadBranches( &tree );
  cont.LoadBranches( &tree );

  for( int i = 0; i < tree.GetEntries(); ++i ){
    tree.GetEntry( i );
    bool pass = true;

    for( const auto& opt : options_list ){
      if( args.CheckArg( opt )
          && FormatOpt( fmt, opt ) != args.Arg<double>( opt ) ){
        pass = false;
        break;
      }
    }

    if( !pass ){
      continue;
    } else {
      entry = i;
      break;
    }
  }

  if( entry < 0 ){
    std::cout << "NO PLOTS SATISFIED" << std::endl;
    return 0;
  }

  tree.GetEntry( entry );
  std::cout << "Dimple Radius:" << fmt.dimple_rad << "[mm]\n"
            << "Dimple Indent:" << fmt.dimple_ind << "[mm]\n"
            << "Tile Width   :" << fmt.tile_x     << "[mm]\n"
            << "Absorption   :" << fmt.abs_mult   << "\n"
            << "Wrap Ref.    :" << fmt.wrap_ref   << "\n"
            << "SiPM Width   :" << fmt.sipm_width << "[mm]\n"
            << "SiPM Stand   :" << fmt.sipm_stand << "[mm]\n"
            << "PCB Radius   :" << fmt.pcb_rad    << "[mm]\n"
            << "PCB Ref      :" << fmt.pcb_ref    << std::endl;

  const std::vector<std::string> compare_list = {
    "OpticalLength",
    "NumWrapReflection",
    "NumPCBReflection"
  };

  for( const auto str : compare_list ){

    usr::plt::Simple1DCanvas c;

    TH1D* all_hist = MakePhotonHist( fmt, cont, str );
    TH1D* det_hist = MakePhotonHist( fmt, cont, str + "_Detected" );

    TH1D* all_acc = MakePhotonAccum( fmt, cont, str );
    TH1D* det_acc = MakePhotonAccum( fmt, cont, str + "_Detected" );

    c.PlotHist( det_acc,
      usr::plt::PlotType( usr::plt::hist ),
      usr::plt::TrackY( usr::plt::TrackY::min ),
      usr::plt::LineColor( usr::plt::col::red ),
      usr::plt::EntryText( "Det. Photons (Remaining)" ) );
    c.PlotHist( det_hist,
      usr::plt::PlotType( usr::plt::hist ),
      usr::plt::TrackY( usr::plt::TrackY::min ),
      usr::plt::LineColor( usr::plt::col::red ),
      usr::plt::LineStyle( usr::plt::sty::lindensedot ),
      usr::plt::EntryText( "Det. Photons" ) );
    c.PlotHist( all_acc,
      usr::plt::PlotType( usr::plt::hist ),
      usr::plt::TrackY( usr::plt::TrackY::min ),
      usr::plt::LineColor( usr::plt::col::black ),
      usr::plt::EntryText( "All Photons (Remaining)" ) );
    c.PlotHist( all_hist,
      usr::plt::PlotType( usr::plt::hist ),
      usr::plt::TrackY( usr::plt::TrackY::min ),
      usr::plt::LineColor( usr::plt::col::black ),
      usr::plt::LineStyle( usr::plt::sty::lindensedot ),
      usr::plt::EntryText( "All Photons" ) );

    const std::string xtitle
      = str == "OpticalLength" ?  "Optical Length [cm]" :
        str == "NumWrapReflection" ? "Number of wrap reflections" :
        "Number of PCB reflections";

    c.Pad().Yaxis().SetTitle( "A.U." );
    c.Pad().Xaxis().SetTitle( xtitle.c_str() );

    c.Pad().SetLogy( kTRUE );
    c.Pad().SetYaxisMin( 1e-4 );
    c.Pad().SetYaxisMax( 33 );

    c.Pad().DrawCMSLabel( "Simulation", "HGCal" );

    for( const auto opt : options_list ){
      if( args.CheckArg( opt ) ){
        c.Pad().WriteLine( FormatOptString( fmt, opt ) );
      }
    }

    c.SaveAsPDF( args.MakePDFFile( "ComparePhoton_" + str ) );
  }

  return 0;
}
