#ifdef CMSSW_GIT_HASH
#include "UserUtils/PlotUtils/interface/PlotCommon.hpp"
#else
#include "UserUtils/PlotUtils/PlotCommon.hpp"
#endif

#include "TProfile.h"

namespace usr {

namespace plt {

/**
 * Returning the maximum bin-value of a histogram object with the bin error
 * taken into consideration.
 */
double
GetYmax( const TH1D* hist )
{
  double ans = 0;

  for( int i = 1; i <= hist->GetNbinsX(); ++i ){
    const double bincont = hist->GetBinContent( i );
    // Special case for TProfile... (not sure why yet)
    const double binerr  = hist->InheritsFrom( TProfile::Class()) ? 0 :
                           hist->GetBinError( i );
    // Skipping over Zero bins
    if( bincont == 0 ) continue;
    ans = std::max( ans, bincont + binerr );
  }

  return ans;
}

/**
 * Returning the minimum bin value of a histogram object, with the bin error
 * taken into account.
 */
double
GetYmin( const TH1D* hist )
{
  double ans = 0.3;

  for( int i = 1; i <= hist->GetNbinsX(); ++i ){
    const double bincont = hist->GetBinContent( i );
    const double binerr  = hist->GetBinError( i );
    // Skipping over Zero bins
    if( bincont == 0 ) continue;
    ans = std::min( ans, bincont - binerr );
  }

  return ans;
}

/**
 * Returning maximum y value stored in a TGraphs object, with the Y error bar
 * taken into account (The GetErrorYhigh method is virtual, so TGraphs with
 * asymmetric errors would be handled properly.)
 */
double
GetYmax( const TGraph* x )
{
  double ans = -std::numeric_limits<double>::max();

  for( int i = 0; i < x->GetN(); ++i ){
    const double bin = x->GetY()[i] + std::max( x->GetErrorYhigh( i ), 0.0 );
    ans = std::max( ans, bin );
  }

  return ans;
}

/**
 * Returning the minimum y value stored in a TGraph object, with the y error
 * bar taken into account. (Again exploiting the fact that GetErrorYlow is
 * virtual.)
 */
double
GetYmin( const TGraph* x )
{
  double ans = std::numeric_limits<double>::max();

  for( int i = 0; i < x->GetN(); ++i ){
    const double bin = x->GetY()[i] - std::max( x->GetErrorYlow( i ), 0.0 );
    ans = std::min( ans, bin );
  }

  return ans;
}

/**
 * Returning the minimum value of x stored in the TGraph object with x errors
 * taken into account.
 */
double
GetXmin( const TGraph* x )
{
  double ans = std::numeric_limits<double>::max();

  for( int i = 0; i < x->GetN(); ++i ){
    const double bin = x->GetX()[i] - std::max( x->GetErrorXlow( i ), 0.0 );
    ans = std::min( ans, bin );
  }

  return ans;
}

/**
 * Returning the maximum value of x stored in the TGraph object with x errors
 * taken into account.
 */
double
GetXmax( const TGraph* x )
{
  double ans = -std::numeric_limits<double>::max();

  for( int i = 0; i < x->GetN(); ++i ){
    const double bin = x->GetX()[i] + std::max( x->GetErrorXhigh( i ), 0.0 );
    ans = std::max( ans, bin );
  }

  return ans;
}

/**
 * Retuning the maximum bin value of a THStack object, the bin errors of each
 * histogram would *not* be taken into account. This functions assumes
 * *identical* binning of the internal histograms, though this might not be the
 * case. Use at your own risk.
 */
double
GetYmax( const THStack* x )
{
  double ans = 0.3;
  if( x->GetNhists() == 0 ){ return ans; }

  TH1* first = (TH1*)x->GetHists()->At( 0 );

  for( int i = 1; i < first->GetNcells(); ++i ){
    double sum = 0;

    for( int j = 0; j < x->GetNhists(); ++j ){
      sum += ( (TH1*)x->GetHists()->At( j ) )->GetBinContent( i );
    }

    ans = std::max( ans, sum );
  }

  return ans;
}

/**
 * Getting the minimum bin value of a THStack object, the bin errors of each
 * internal histogram would *not* be taken into account. This function currently
 * assume *identical* binning for the histograms, though this is not necessary
 * the case. Use at your own risk.
 */
double
GetYmin( const THStack* x )
{
  double ans = 0.3;
  if( x->GetNhists() == 0 ){ return ans; }

  TH1* first = (TH1*)x->GetHists()->At( 0 );

  for( int i = 1; i < first->GetNcells(); ++i ){
    double sum = 0;

    for( int j = 0; j < x->GetNhists(); ++j ){
      sum += ( (TH1*)x->GetHists()->At( j ) )->GetBinContent( i );
    }

    ans = std::min( ans, sum );
  }

  return ans;
}

}

}
