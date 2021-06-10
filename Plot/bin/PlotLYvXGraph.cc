#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/Common/interface/Format.hpp"
#include "UserUtils/Common/interface/STLUtils/Filesystem.hpp"
#include "UserUtils/Common/interface/STLUtils/VectorUtils.hpp"
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"
#include "UserUtils/PlotUtils/interface/Ratio1DCanvas.hpp"

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

  const std::vector<int> filllist = {
    usr::plt::col::gray,
    usr::plt::col::lightblue,
    usr::plt::col::pink,
    usr::plt::col::lightgreen,
    usr::plt::col::lightgoldenrodyellow,
    usr::plt::col::palevioletred
  };


  usr::plt::Ratio1DCanvas cr;
  usr::plt::Simple1DCanvas c( usr::plt::len::a4textwidth_default(),
                              usr::plt::len::a4textwidth_default(),
                              usr::plt::FontSet( 12 ) );
  usr::plt::Simple1DCanvas c_opt_length( usr::plt::len::a4textwidth_default(),
                                         usr::plt::len::a4textwidth_default(),
                                         usr::plt::FontSet( 12 ) );
  usr::plt::Simple1DCanvas c_num_ref( usr::plt::len::a4textwidth_default(),
                                      usr::plt::len::a4textwidth_default(),
                                      usr::plt::FontSet( 12 ) );
  usr::plt::Simple1DCanvas c_num_pcb( usr::plt::len::a4textwidth_default(),
                                      usr::plt::len::a4textwidth_default(),
                                      usr::plt::FontSet( 12 ) );
  TH1D dummyhist( "", "", 1, -xmax - 5, xmax + 5 );
  c.PlotHist( dummyhist, usr::plt::PlotType( usr::plt::hist ) );

  unsigned idx = 0;
  double width = 0;

  graph_container.LoadBranches( &tree );

  TGraph* front = nullptr;

  for( const auto& runentry : entrylist ){
    tree.GetEntry( runentry );

    TGraph* graph = MakeGraph( fmt, graph_container, "LYvX" );
    TH1D* opt_len = MakePhotonAccum( fmt, graph_container, "OpticalLength" );
    TH1D* num_ref = MakePhotonAccum( fmt, graph_container, "NumWrapReflection" );
    TH1D* num_pcb = MakePhotonAccum( fmt, graph_container, "NumPCBReflection" );

    std::string entry = entrylist.size() > 0 ? "" : "Simulation Results";
    if( front == nullptr ){
      front = graph ;
    }

    for( const auto opt : difflist ){
      entry += FormatOptString( fmt, opt ) + ";";
    }

    entry.pop_back();

    // if( idx >= 5 ){ continue; }
    const auto color = colorlist.at( idx % colorlist.size() );
    const auto fill  = filllist.at( idx % filllist.size() );

    c.PlotGraph( graph,
      usr::plt::PlotType( usr::plt::fittedfunc ),
      usr::plt::TrackY( usr::plt::TrackY::max ),
      usr::plt::LineColor( color ),
      usr::plt::FillColor( fill ),
      usr::plt::FillStyle( usr::plt::sty::fillsolid ),
      usr::plt::EntryText( entry ) );

    cr.PlotGraph( graph,
      usr::plt::PlotType( usr::plt::fittedfunc ),
      usr::plt::TrackY( usr::plt::TrackY::max ),
      usr::plt::LineColor( color ),
      usr::plt::FillColor( fill ),
      usr::plt::FillStyle( usr::plt::sty::fillsolid ),
      usr::plt::EntryText( entry ) );
    cr.PlotScale(graph,front,
      usr::plt::PlotType( usr::plt::fittedfunc ),
      usr::plt::TrackY( usr::plt::TrackY::max ),
      usr::plt::LineColor( color ),
      usr::plt::FillColor( fill ),
      usr::plt::FillStyle( usr::plt::sty::fillsolid ) );



    width = graph->GetErrorX( 0 );
    std::cout << "Plot" << std::endl;


    c_opt_length.PlotHist( opt_len,
      usr::plt::PlotType( usr::plt::hist ),
      usr::plt::TrackY( usr::plt::TrackY::min ),
      usr::plt::LineColor( color ),
      usr::plt::EntryText( entry ) );
    c_num_ref.PlotHist( num_ref,
      usr::plt::PlotType( usr::plt::hist ),
      usr::plt::TrackY( usr::plt::TrackY::min ),
      usr::plt::LineColor( color ),
      usr::plt::EntryText( entry ) );
    c_num_pcb.PlotHist( num_pcb,
      usr::plt::PlotType( usr::plt::hist ),
      usr::plt::TrackY( usr::plt::TrackY::min ),
      usr::plt::LineColor( color ),
      usr::plt::EntryText( entry ) );
    ++idx;
  }

  if( !usr::FindValue( difflist,  std::string( "dimplerad" ) ) ){
    tree.GetEntry( entrylist.front() );
    c.Pad().DrawVLine( fmt.dimple_rad,
      usr::plt::LineColor( usr::plt::col::gray ),
      usr::plt::LineStyle( usr::plt::sty::lindotted ) );
    c.Pad().DrawVLine( -fmt.dimple_rad,
      usr::plt::LineColor( usr::plt::col::gray ),
      usr::plt::LineStyle( usr::plt::sty::lindotted ) );
  }

  if( !usr::FindValue( difflist, std::string( "tilewidth" ) ) ){
    tree.GetEntry( entrylist.front() );
    c.Pad().DrawVLine( fmt.tile_x/2,
      usr::plt::LineColor( usr::plt::col::gray ),
      usr::plt::LineStyle( usr::plt::sty::lindotted ) );
    c.Pad().DrawVLine( -fmt.tile_x/2,
      usr::plt::LineColor( usr::plt::col::gray ),
      usr::plt::LineStyle( usr::plt::sty::lindotted ) );
  }

  c.Pad().SetDataMax( c.Pad().GetDataMax() * 1.3 );
  c.Pad().Xaxis().SetTitle( "Beam Center X [mm]" );
  c.Pad().Yaxis().SetTitle( "Detected Photons" );

  cr.BottomPad().Xaxis().SetTitle( "Beam Center X [mm]" );
  cr.BottomPad().Yaxis().SetTitle( "Ratio" );
  cr.TopPad().Yaxis().SetTitle( "Detected Photons" );
  cr.BottomPad().DrawHLine(0.5);
  cr.BottomPad().DrawHLine(1.5);

  c_opt_length.Pad().Xaxis().SetTitle( "Optical Length [cm]" );
  c_opt_length.Pad().Yaxis().SetTitle( "Remaining Photon fraction" );
  c_num_ref.Pad().Xaxis().SetTitle( "Number of wrap reflections" );
  c_num_ref.Pad().Yaxis().SetTitle( "Remaining Photon fraction" );
  c_num_pcb.Pad().Xaxis().SetTitle( "Number of PCB reflections" );
  c_num_pcb.Pad().Yaxis().SetTitle( "Remaining Photon fraction" );

  c_opt_length.Pad().SetDataMin( 1e-3 );
  c_opt_length.Pad().SetLogy( kTRUE );
  c_num_ref.Pad().SetDataMin( 1e-3 );
  c_num_ref.Pad().SetLogy( kTRUE );
  c_num_pcb.Pad().SetDataMin( 1e-3 );
  c_num_pcb.Pad().SetLogy( kTRUE );

  c.DrawLuminosity( usr::fstr( "Trigger Width: %.1lf[mm]", width ) );
  c.DrawCMSLabel( "Simulation", "HGCal" );
  c_opt_length.Pad().DrawCMSLabel( "Simulation", "HGCal" );
  c_num_ref.Pad().DrawCMSLabel( "Simulation", "HGCal" );
  c_num_pcb.Pad().DrawCMSLabel( "Simulation", "HGCal" );

  for( const auto opt : options_list ){
    tree.GetEntry( entrylist.front() );
    if( args.CheckArg( opt ) ){
      c.Pad().WriteLine( FormatOptString( fmt, opt ) );
      c_opt_length.Pad().WriteLine( FormatOptString( fmt, opt ) );
      c_num_ref.Pad().WriteLine( FormatOptString( fmt, opt ) );
      c_num_pcb.Pad().WriteLine( FormatOptString( fmt, opt ) );
    }
  }

  std::cout << args.MakePDFFile( "LYvX" ) << std::endl;
  c.SaveAsPDF( args.MakePDFFile( "LYvX" ) );
  cr.SaveAsPDF( args.MakePDFFile("LYvX_ratio") );
  c_opt_length.SaveAsPDF( args.MakePDFFile( "OptLen" ) );
  c_num_ref.SaveAsPDF( args.MakePDFFile( "NumWrap" ) );
  c_num_pcb.SaveAsPDF( args.MakePDFFile( "NumPCB" ) );
}
