#ifdef CMSSW_GIT_HASH
#include "UserUtils/PlotUtils/interface/PlotCommon.hpp"
#else
#include "UserUtils/PlotUtils/PlotCommon.hpp"
#endif

#include "CmdSetAttr.hpp"

#include <iostream>

extern void
SetTextAttr( TAttText&                   obj,
             const usr::RooArgContainer& args )
{
  if( args.Has<usr::plt::TextColor>() ){
    const auto cmd = args.Get<usr::plt::TextColor>();
    obj.SetTextColorAlpha( cmd.Col(), cmd.Alpha() );
  }

  if( args.Has<usr::plt::TextSize>() ){
    obj.SetTextSize( args.Get<usr::plt::TextSize>() );
  }
}

extern void
SetLineAttr( TAttLine&                   obj,
             const usr::RooArgContainer& args  )
{
  if( args.Has<usr::plt::LineColor>() ){
    const auto cmd = args.Get<usr::plt::LineColor>() ;
    obj.SetLineColorAlpha( cmd.Col(), cmd.Alpha()  );
  }

  if( args.Has<usr::plt::LineStyle>() ){
    obj.SetLineStyle( args.Get<usr::plt::LineStyle>() );
  }

  if( args.Has<usr::plt::LineWidth>() ){
    obj.SetLineWidth( args.Get<usr::plt::LineWidth>() );
  }
}

extern void
SetFillAttr( TAttFill&                   obj,
             const usr::RooArgContainer& args  )
{
  if( args.Has<usr::plt::FillColor>() ){
    const auto cmd = args.Get<usr::plt::FillColor>() ;
    obj.SetFillColorAlpha( cmd.Col(), cmd.Alpha() );
  }

  if( args.Has<usr::plt::FillStyle>() ){
    obj.SetFillStyle( args.Get<usr::plt::FillStyle>() );
  }
}

extern void
SetMarkAttr( TAttMarker&                 obj,
             const usr::RooArgContainer& args  )
{
  if( args.Has<usr::plt::MarkerColor>() ){
    const auto cmd = args.Get<usr::plt::MarkerColor>();
    obj.SetMarkerColorAlpha( cmd.Col(), cmd.Alpha() );
  }

  if( args.Has<usr::plt::MarkerStyle>() ){
    obj.SetMarkerStyle( args.Get<usr::plt::MarkerStyle>() );
  }

  if( args.Has( usr::plt::MarkerSize::CmdName ) ){
    obj.SetMarkerSize( args.Get<usr::plt::MarkerSize>() );
  }
}
