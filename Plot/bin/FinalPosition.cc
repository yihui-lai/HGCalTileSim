#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/Common/interface/Format.hpp"
#include "UserUtils/Common/interface/STLUtils/Filesystem.hpp"
#include "UserUtils/Common/interface/STLUtils/VectorUtils.hpp"
#include "UserUtils/PlotUtils/interface/Flat2DCanvas.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"

#include "HGCalTileSim/Plot/interface/LYvX_Common.hpp"

#include "TChain.h"
#include "TFile.h"
#include "TGraphErrors.h"

int
main( int argc, char** argv )
{
  usr::po::options_description desc( "Input output options" );
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
    {"dimplerad", "r"}, {"dimpleind", "d"}, {"dimpletype", "T"},
    {"tilewidth", "L"}, {"absmult", "a"},
    {"wrapreflect", "m"},
    {"sipmwidth", "W"}, {"sipmstand", "S"},
    {"tilealpha", "A"}, {"dimplealpha", "D"},
    {"PCBRadius", "b"}, {"PCBRef", "P"}
  } );

  const std::vector<std::string> options_list = {
    "dimplerad", "dimpleind",     "dimpletype",
    "tilewidth",
    "absmult",   "wrapreflect",
    "sipmwidth", "sipmstand",
    "PCBRadius", "PCBRef",
    "tilealpha", "dimplealpha"
  };
  TChain tree( "LYSimRun", "LYSimRun" );
  LYSimRunFormat fmt;
  LYSimRunFormat first_fmt;
  LYvXGraphContainer graph_container;
  std::vector<std::string> difflist;
  std::vector<int> entrylist;
  double xmax = 0;

  tree.Add( args.Arg<std::string>( "inputfile" ).c_str() );
  fmt.LoadBranches( &tree );

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

    if( !pass ){continue;}

    if( entrylist.size() == 0 ){
      first_fmt = fmt;
      xmax      = std::max( fmt.tile_x/2, xmax );
    }
    entrylist.push_back( i );

    for( const auto& opt : options_list ){
      if( entrylist.size() > 0 &&
          FormatOpt( fmt, opt ) != FormatOpt( first_fmt, opt ) ){
        difflist.push_back( opt );
      }
    }
  }

  if( entrylist.size() == 0 ){
    std::cout << "NO PLOTS SATISFIED" << std::endl;
    return 0;
  }

  difflist = usr::RemoveDuplicate( difflist );

  const std::vector<int> colorlist = {
    usr::plt::col::black,
    usr::plt::col::darkblue,
    usr::plt::col::darkred,
    usr::plt::col::darkgreen,
    usr::plt::col::darkorange,
    usr::plt::col::darkviolet
  };

  usr::plt::Simple1DCanvas c( usr::plt::len::a4textwidth_default(),
                              usr::plt::len::a4textwidth_default(),
                              usr::plt::FontSet( 12 ) );
  unsigned idx = 0;
  double width = 0;

  graph_container.LoadBranches( &tree );

  for( const auto& runentry : entrylist ){
    tree.GetEntry( runentry );

    TH1D* hist = MakeDetXHist( fmt, graph_container  );

    std::string entry = entrylist.size() > 0 ? "" : "Simulation Results";

    for( const auto opt : difflist ){
      entry += FormatOptString( fmt, opt ) + ";";
    }

    entry.pop_back();

    // if( idx >= 5 ){ continue; }
    const auto color = colorlist.at( idx % colorlist.size() );

    c.PlotHist( hist,
      usr::plt::PlotType( usr::plt::hist ),
      usr::plt::TrackY( usr::plt::TrackY::both ),
      usr::plt::LineColor( color ),
      usr::plt::EntryText( entry ) );

    ++idx;
  }

  c.Pad().SetDataMax( c.Pad().GetDataMax() * 1.3 );
  c.Pad().Xaxis().SetTitle( "Detected photon final x position X [mm]" );
  c.Pad().Yaxis().SetTitle( "Normalized num. of detected photons" );

  c.DrawLuminosity( usr::fstr( "Trigger Width: %.1lf[mm]", width ) );
  c.DrawCMSLabel( "Simulation", "HGCal" );

  for( const auto opt : options_list ){
    tree.GetEntry( entrylist.front() );
    if( args.CheckArg( opt ) ){
      c.Pad().WriteLine( FormatOptString( fmt, opt ) );
    }
  }

  std::cout << args.MakePDFFile( "FinalPositionX_Detected" ) << std::endl;
  c.SaveAsPDF( args.MakePDFFile( "FinalPositionX_Detected" ) );

  return 0;

}
