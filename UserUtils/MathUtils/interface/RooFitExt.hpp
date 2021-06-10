/**
 * @file   RooFitExt.hpp
 *  @brief  Extension to some mathematical computations for RooFit objects
 *  @author Yi-Mu "Enoch" Chen [ ensc@hep1.phys.ntu.edu.tw ]
 */
#ifndef USERUTILS_MATHUTILS_ROOFITEXT_HPP
#define USERUTILS_MATHUTILS_ROOFITEXT_HPP

#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/STLUtils/VectorUtils.hpp"
#else
#include "UserUtils/Common/STLUtils/VectorUtils.hpp"
#endif

#include <RooAbsPdf.h>
#include <RooCmdArg.h>
#include <RooDataSet.h>

namespace usr {

/**
 * @defgroup KSRooFit KSRooFit
 * @brief    Kolmogorov--Smirnov test routines for RooFit objects
 * @ingroup MathUtils
 * @details
 * For some strange reason, there doesn't exist the [Kolmogorov--Smirnov test](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test)
 * routines doesn't exit for RooFit objects. The functions provides includes
 * the test for both data--model comparison and data--data comparison, with
 * up to 2 variable range cut for regional and side-band goodness-of-fit tests.
 * Test includes both the maximum cumulative distribution function distances
 * and the resulting p-value when comparing the distance to the Kolmogorov
 * distribution
 * @{
 */
extern double KSDistance(
  RooDataSet&      dataset,
  RooAbsPdf&       pdf,
  RooRealVar&      var,
  const RooCmdArg& cut1 = RooCmdArg::none(),
  const RooCmdArg& cut2 = RooCmdArg::none()
  );

extern double KSProb(
  RooDataSet&      dataset,
  RooAbsPdf&       pdf,
  RooRealVar&      var,
  const RooCmdArg& cut1 = RooCmdArg::none(),
  const RooCmdArg& cut2 = RooCmdArg::none()
  );

extern double KSDistance(
  RooDataSet&      set1,
  RooDataSet&      set2,
  RooRealVar&      var,
  const RooCmdArg& cut1 = RooCmdArg::none(),
  const RooCmdArg& cut2 = RooCmdArg::none()
  );

extern double KSProb(
  RooDataSet&      set1,
  RooDataSet&      set2,
  RooRealVar&      var,
  const RooCmdArg& cut1 = RooCmdArg::none(),
  const RooCmdArg& cut2 = RooCmdArg::none()
  );

extern double KSProbAlt(
  RooDataSet&      set1,
  RooDataSet&      set2,
  RooRealVar&      var,
  const RooCmdArg& cut1 = RooCmdArg::none(),
  const RooCmdArg& cut2 = RooCmdArg::none()
  );
/* @} */

/**
 * @defgroup RooFitEnhance RooFitEnhance
 * @ingroup  MathUtils
 * @brief    Nicer routines for running RooFit routines functions in
 * @details
 * RooFit Routines are sometimes excessively verbose to type out, or require too
 * much intervention. These functions aim to trim down the verbosity of writing
 * certain functions.
 * @{
 */

/**
 * @{
 * @brief Fitting routing with arbitrary many arguments and disables screen
 * output by default.
 */
extern RooFitResult* FitPDFToData(
  RooAbsPdf&, RooAbsData&,
  const std::vector<RooCmdArg>& cmdargs
  );

inline RooFitResult*
FitPDFToData(
  RooAbsPdf& pdf, RooAbsData& data ){ return FitPDFToData( pdf, data, {} ); }

template<typename ... Args>
inline RooFitResult*
FitPDFToData(
  RooAbsPdf& pdf, RooAbsData& data, const RooCmdArg & arg1, Args ... args )
{
  return FitPDFToData( pdf, data, MakeVector<RooCmdArg>( arg1, args ... ) );
}
/* @} */

class MaxFitIteration : public RooCmdArg
{
public:
  static const std::string CmdName;
  MaxFitIteration( const RooCmdArg* x ) : RooCmdArg( *x ){}
  MaxFitIteration( unsigned x );
  virtual
  ~MaxFitIteration(){}
  inline
  operator unsigned() const { return getInt( 0 ); }
};

/**
 * @{
 * @brief Running the Fit Routine until converge or pass certain amount of fit
 * calls.
 */
extern RooFitResult* ConvergeFitPDFToData(
  RooAbsPdf&, RooAbsData&,
  const std::vector<RooCmdArg>& cmdargs
  );

inline RooFitResult*
ConvergeFitPDFToData(
  RooAbsPdf& pdf, RooAbsData& data )
{
  return ConvergeFitPDFToData( pdf, data, {} );
}

template<typename ... Args>
inline RooFitResult*
ConvergeFitPDFToData(
  RooAbsPdf& pdf, RooAbsData& data, const RooCmdArg & arg1, Args ... args )
{
  return ConvergeFitPDFToData( pdf, data,
    MakeVector<RooCmdArg>( arg1, args ... ) );
}
/** @} */

/** @}  */

}/* usr */

#endif/* end of include guard: USERUTILS_MATHUTILS_ROOFITEXT_HPP */
