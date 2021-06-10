/**
 * @file RooFitExt_KSText.cc
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 * @brief   Implementing extended RooFit functions for KSTest computation
 */

#ifdef CMSSW_GIT_HASH
#include "UserUtils/MathUtils/interface/RooFitExt.hpp"
#else
#include "UserUtils/MathUtils/RooFitExt.hpp"
#endif

#include <RooAbsData.h>
#include <RooAbsPdf.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <TMath.h>

#include <boost/format.hpp>
#include <memory>

namespace usr {

/**
 * @brief Helper structure for expressing a 1D data set (RooDataSet) simply
 *        as a list of numbers with weights.
 *
 * Additional helper routines are used for caching anc calculating the numbers
 * required for performing the Kolmogorov-Smirnov test. One in particular is the
 * "EffectiveNum" method, as the number of data points is required for the
 * calculation of the KS test, the "number of data points" for weighted data
 * sets would follow the "sum squared over sum of squared" approach:
 *
 * \f[
 *  N_\mathrm{eff} = \frac{\left(\sum w_i\right)^2 }{\sum w_i^2}
 * \f]
 *
 * The data set is also simplified to a list of doubles (with weights), with
 * the list ordered to make it simpler to calculated the KS test.
 */
struct SimplifiedData
{
  typedef std::pair<double, double> Evt;
  double           sum;   /**< @brief Sum of weights */
  double           sumsq; /**< @brief Sum of squared weights */
  std::vector<Evt> dataset;

  inline static bool
  Compare( const Evt& x, const Evt& y )
  { return x.first < y.first; }

  double EffectiveNum() const;

  SimplifiedData(
    RooDataSet&      data,
    RooRealVar&      x,
    const RooCmdArg& cut1,
    const RooCmdArg& cut2
    );
};

/**
 * @brief Helper structure for expressing a CDF object, and helping with
 *        adjusting the normalization of the CDF function when the Kolmogorov--
 *        Smirnov test is asked to be performed on certain regions of the
 *        variable range.
 *
 * The construction of the CDF functions as a RooFit RooAbsRel object handles
 * the normalization of the CDF function across the maximum range of the
 * variable used to constructed the CDF function. In the case that a range
 * is requested when calculating the KS test would require the effective CDF is
 * 0/1 at the edges of the range, so effectively:
 * \f[
 *     \mathrm{CDF}_{\mathrm{eff}}(x)
 *      = \frac{\mathrm{CDF}(x) -\mathrm{CDF}(x_\mathrm{min})}
 *             {\mathrm{CDF}(x_\mathrm{max})-\mathrm{CDF}(x_\mathrm{max})}
 * \f]
 *
 * In the case that two ranges is requested (say to test the goodness-of-fit
 * for a side band fit result), the effective CDF is construted such that CDF
 * is continuous for all values of \f$x\f$, an returns values 0/1 at the far
 * edges of the ranges. So effectively, given the two ranges \f$[a,b]\f$ and
 * \f$[c,d]\f$:
 *
 * \f[
 *  \mathrm{CDF}_\mathrm{eff}(x < b)
 *      = \frac{\mathrm{CDF}(x) - \mathrm{CDF}(a)}
 *      {\mathrm{CDF}(b)-\mathrm{CDF}(a) + \mathrm{CDF}(d)-\mathrm{CDF}(c)};
 *
 *  \mathrm{CDF}_\mathrm{eff}(x > c)
 *      = \frac
 *      {\mathrm{CDF}(x)-\mathrm{CDF}(c) + \mathrm{CDF}(b)-\mathrm{CDF}(a)}
 *      {\mathrm{CDF}(b)-\mathrm{CDF}(a) + \mathrm{CDF}(d)-\mathrm{CDF}(c)}
 * \f]
 *
 *
 * While this scheme readily extends to arbitrarily many ranges, our
 * implementation, would only allow for two ranges.
 */
struct SimplifiedCDF
{
  std::unique_ptr<RooAbsReal>             rawcdf;
  RooRealVar&                             var;
  std::vector<std::pair<double, double> > rangelist;
  double                                  norm;

  double operator()( const double x );
  double rawval( const double x );

  inline double
  min( unsigned x ){ return rangelist.at( x ).first; }
  inline double
  max( unsigned x ){ return rangelist.at( x ).second; }

  SimplifiedCDF(
    RooAbsPdf&       pdf,
    RooRealVar&      x,
    const RooCmdArg& cut1,
    const RooCmdArg& cut2
    );
  void AddRange( const RooCmdArg& cmd );
};

bool CheckCutCmd( const RooCmdArg& cmd, RooRealVar& x );


/**
 * @brief Calculating the Kolmogorov--Smirov distance of a RooDataSet and a
 *        RooAbsPdf in term of a variable.
 *
 * The calculation routine supports:
 *  - weighted data sets. To see how the number of data is calculated for
 *    weighted data sets, see the helper class SimplifiedData.
 *  - Up to two cuts in the data set. Currently you can only cut by specifying
 *    a range already declared in the variable of comparison. To see how
 *    the cumulative PDF is handled, see the helper class SimplifiedCDF.
 */
double
KSDistance(
  RooDataSet&      dataset,
  RooAbsPdf&       pdf,
  RooRealVar&      var,
  const RooCmdArg& cut1,
  const RooCmdArg& cut2
  )
{
  SimplifiedData sim( dataset, var, cut1, cut2 );
  SimplifiedCDF cdf( pdf, var, cut1, cut2 );

  double empcdf  = 0;
  double maxdist = 0;

  for( const auto& evt : sim.dataset ){
    // emperical cdf is the sum of weights
    // (normalization handled latter for stability )
    const double val = evt.first;
    const double wgt = evt.second;
    empcdf += wgt;
    maxdist = std::max( maxdist, fabs( cdf( val ) - ( empcdf/sim.sum ) ) );
  }

  return sqrt( sim.EffectiveNum() ) * maxdist;
}

/**
 * @brief Returning the corresponding KS probability after calculating the KS
 *        distance of between a data set and PDF.
 *
 * The probability function is already implemented in the
 * TMath::KolmogorovProb() function, will not attempt to recreate.
 */
double
KSProb(
  RooDataSet&      dataset,
  RooAbsPdf&       pdf,
  RooRealVar&      var,
  const RooCmdArg& cut1,
  const RooCmdArg& cut2
  )
{
  const double dist = KSDistance( dataset, pdf, var, cut1, cut2 );
  return TMath::KolmogorovProb( dist );
}

/**
 * @brief Calculating the Kolmogorov-Smirnov distance between two data sets.
 *
 * The code for looping though the two data set is based on the ROOT
 * implementation found in the TMath::KolmogorovTest() functions.
 * Supporting up to two range cuts and weighted data sets. For the calculation
 * of the effective number of events in a weighted data set, see the helper
 * class SimplifiedData.
 */
double
KSDistance(
  RooDataSet&      set1,
  RooDataSet&      set2,
  RooRealVar&      var,
  const RooCmdArg& cut1,
  const RooCmdArg& cut2
  )
{
  SimplifiedData sim1( set1, var, cut1, cut2 );
  SimplifiedData sim2( set2, var, cut1, cut2 );

  double empcdf1 = 0;
  double empcdf2 = 0;
  double maxdist = 0;

  for( auto evt1 = sim1.dataset.begin(), evt2 = sim2.dataset.begin();
       evt1 != sim1.dataset.end() && evt2 != sim2.dataset.end(); ){
    if( evt1->first < evt2->first ){
      empcdf1 += evt1->second;
      evt1++;
    } else if( evt1->first > evt2->first ){
      empcdf2 += evt2->second;
      evt2++;
    } else {// special case for identical value
      double x = evt1->first;

      while( evt1->first == x && evt1 != sim1.dataset.end() ){
        empcdf1 += evt1->second;
        evt1++;
      }

      while( evt2->first == x && evt2 != sim2.dataset.end() ){
        empcdf2 += evt2->second;
        evt2++;
      }
    }
    maxdist = TMath::Max( maxdist,
      fabs( ( empcdf1/sim1.sum ) - ( empcdf2/sim2.sum ) ) );
  }

  const double num1 = std::max( sim1.EffectiveNum(), sim2.EffectiveNum() );
  const double num2 = std::min( sim1.EffectiveNum(), sim2.EffectiveNum() );
  return maxdist * sqrt( ( num1/( num1 + num2 ) )*num2 );

}

/**
 * @brief Returning the KS probability of two data sets
 * @details this takes the KS distance of the two datasets and pass it directly
 * to the TMath::KolmogorovProb() functions.
 */
double
KSProb(
  RooDataSet&      set1,
  RooDataSet&      set2,
  RooRealVar&      var,
  const RooCmdArg& cut1,
  const RooCmdArg& cut2 )
{
  const double dist = KSDistance( set1, set2, var, cut1, cut2 );
  return TMath::KolmogorovProb( dist );
}

/**
 * @brief An alternative, effective p-value of the KS test for two data set.
 * @details the alternative calculation is described in
  [this]( https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test#Two-sample_Kolmogorov%E2%80%93Smirnov_test) page.
  Hard limiting the calculation output to 1 for a more "sane" output value
 */
double
KSProbAlt(
  RooDataSet&      set1,
  RooDataSet&      set2,
  RooRealVar&      var,
  const RooCmdArg& cut1,
  const RooCmdArg& cut2 )
{
  const double dist = KSDistance( set1, set2, var, cut1, cut2 );

  // Alternative calculation of p-value as described in wiki page
  return std::min( exp( -2. * dist * dist ) * 2., 1. );
}

/*-----------------------------------------------------------------------------
 *  SimplifiedData function implementations
   --------------------------------------------------------------------------*/
SimplifiedData::SimplifiedData(
  RooDataSet&      data,
  RooRealVar&      x,
  const RooCmdArg& cut1,
  const RooCmdArg& cut2 ) :
  sum( 0 ),
  sumsq( 0 )
{
  RooDataSet* op_set  = 0;
  RooDataSet* tmp_set = 0;
  if( !CheckCutCmd( cut1, x ) && !CheckCutCmd( cut2, x ) ){
    op_set = (RooDataSet*)( data.reduce( RooFit::SelectVars(x) ));
  } else if( CheckCutCmd( cut1, x ) && !CheckCutCmd( cut2, x ) ){
    op_set = (RooDataSet*)( data.reduce( RooFit::SelectVars( x ), cut1 ) );
  } else if( CheckCutCmd( cut2, x ) && !CheckCutCmd( cut1, x ) ){
    op_set = (RooDataSet*)( data.reduce( RooFit::SelectVars( x ), cut2 ) );
  } else {
    op_set  = (RooDataSet*)( data.reduce( RooFit::SelectVars( x ), cut1 ) );
    tmp_set = (RooDataSet*)( data.reduce( RooFit::SelectVars( x ), cut2 ) );
    op_set->append( *tmp_set );
  }

  if( !op_set->isWeighted() ){
    sum = sumsq = op_set->numEntries();
  }

  for( int i = 0; i < op_set->numEntries(); ++i ){
    const double d = op_set->get( i )->getRealValue( x.GetName() );
    const double w = op_set->weight();
    if( !op_set->isWeighted() ){
      dataset.push_back( std::pair<double, double>( d, 1 ) );
    } else {
      dataset.push_back( std::pair<double, double>( d, w ) );
      sum   += w;
      sumsq += w*w;
    }
  }

  delete op_set;
  if( tmp_set ){ delete tmp_set; }

  std::sort( dataset.begin(), dataset.end(), Compare );
}

/*----------------------------------------------------------------------------*/

double
SimplifiedData::EffectiveNum() const
{
  if( sum == sumsq ){// unweighted data set
    return sum;
  } else {
    // Effective number of events for weighted data set found here:
    // https://stats.stackexchange.com/questions/193439/two-sample-kolmogorov-smirnov-test-with-weights?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
    return sum * sum / sumsq;
  }
}

/*-----------------------------------------------------------------------------
 *  SimplifiedCDF implementations
   --------------------------------------------------------------------------*/
SimplifiedCDF::SimplifiedCDF(
  RooAbsPdf&       pdf,
  RooRealVar&      x,
  const RooCmdArg& cut1,
  const RooCmdArg& cut2 ) :
  rawcdf( pdf.createCdf( RooArgSet( x ) ) ),
  var( x )
{
  if( CheckCutCmd( cut1, var ) ){ AddRange( cut1 ); }
  if( CheckCutCmd( cut2, var ) ){ AddRange( cut2 ); }

  if( rangelist.size() == 2 ){
    // Making sure the lower edge is the first range entry.
    if( min( 0 ) > min( 1 ) ){
      auto temp = rangelist.at( 0 );
      rangelist.at( 0 ) = rangelist.at( 1 );
      rangelist.at( 1 ) = temp;
    }

    // reducing the second range if overlaps
    if( max( 0 ) >= min( 1 ) ){
      rangelist.at( 0 ).second = rangelist.at( 1 ).second;
      rangelist.erase( rangelist.begin() + 1 );
    }
  }

  // Saving normalization denominator
  if( rangelist.size() == 0 ){
    norm = rawval( var.getMax() )  - rawval( var.getMin() );// Should be 1
  } else if( rangelist.size() == 1 ){
    norm = rawval( max( 0 ) ) - rawval( min( 0 ) );
  } else {
    norm = ( rawval( max( 0 ) ) - rawval( min( 0 ) ) )
           + ( rawval( max( 1 ) ) - rawval( min( 1 ) ) );
  }
}

/**
 * @brief Returning the effective CDF value with the cut ranges taken into
 *        account.
 */
double
SimplifiedCDF::operator()( const double x )
{
  double ans = rawval( x );
  if( rangelist.size() == 0 ){
    return ans;
  } else if( rangelist.size() == 1 ){
    if( min( 0 ) < x && x < max( 0 ) ){
      return ( ans - rawval( min( 0 ) ) )/norm;
    } else {
      return ans;
    }
  } else {
    if( min( 0 ) < x && x < max( 0 ) ){
      return ( ans - rawval( min( 0 ) ) )/norm;
    } else if( min( 1 ) < x && x < max( 1 ) ){
      return ( ans - rawval( min( 1 ) ) + ( rawval( max( 0 ) )-rawval( min( 0 ) ) ) )/norm;
    } else {
      return ans;
    }
  }
}

/**
 * @brief returning the raw CDF functions value as constructed by RooFit.
 */
double
SimplifiedCDF::rawval( const double x )
{
  var = x;
  return rawcdf->getVal();
}

/*----------------------------------------------------------------------------*/

void
SimplifiedCDF::AddRange( const RooCmdArg& cmd )
{
  rangelist.push_back( std::pair<double, double>(
    var.getMin( cmd.getString( 0 ) ),
    var.getMax( cmd.getString( 0 ) ) )
    );
}

/*----------------------------------------------------------------------------*/

bool
CheckCutCmd( const RooCmdArg& cmd, RooRealVar& var )
{
  if( std::string( cmd.GetName() ) != "CutRange" ){
    return false;
  } else if( !var.hasRange( cmd.getString( 0 ) ) ){
    return false;
  }
  return true;
}


}/* usr */
