#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/RootUtils/RooArgContainer.hpp"
#include "UserUtils/MathUtils/interface/RooFitExt.hpp"
#else
#include "UserUtils/Common/RootUtils/RooArgContainer.hpp"
#include "UserUtils/MathUtils/RooFitExt.hpp"
#endif

#include <RooAbsPdf.h>
#include <RooFitResult.h>

namespace usr
{

/**
 * @details
 * This uses the variadic interface to allow for arbitrarily many RooCmdArgs to
 * be added to the fit (instead of being hard capped at just 10 RooCmdArgs). By
 * default, it also adds all of the variables required to suppress the outputs.
 * Users can still enable the RooFit messages by explicitly adding back the
 * commands such as RooFit::Verbose( kTRUE ).
 */
extern RooFitResult*
FitPDFToData(
  RooAbsPdf&                    pdf,
  RooAbsData&                   data,
  const std::vector<RooCmdArg>& cmdargs )
{
  const usr::RooArgContainer args(
    cmdargs, {// Explicitly adding options to reduce verbosity
      RooFit::Verbose( false ),
      RooFit::PrintLevel( -100 ),
      RooFit::Warnings( false ),
      RooFit::PrintEvalErrors( -100 )  } );
  RooLinkedList roolist = args.MakeRooList();

  return pdf.fitTo( data, roolist );
}

USRUTILS_COMMON_REGISTERCMD( MaxFitIteration );
MaxFitIteration::MaxFitIteration( unsigned x ) :
  RooCmdArg( MaxFitIteration::CmdName.c_str(), x ){}

/**
 * @details
 * Running the PDF fit routine multiple times until the resulting RooFitResults
 * returns a nominal exit code (status() == 0). You can use the New RooCmdArg
 * `usr::MaxFitIteration` fo specify the maximum number of fit interation to
 * perform (defaults to 10).
 */
extern RooFitResult*
ConvergeFitPDFToData(
  RooAbsPdf& pdf, RooAbsData& data,
  const std::vector<RooCmdArg>& cmdargs
  )
{
  RooFitResult* status = nullptr;
  const usr::RooArgContainer args(
    cmdargs, {// Explicitly adding options to reduce verbosity
      RooFit::Save(),
      MaxFitIteration( 10 )// Default to at most 10
    } );
  const usr::RooArgContainer og_args( cmdargs );

  RooLinkedList roolist = args.MakeRooList();

  for( unsigned i = 0; i < args.Get<MaxFitIteration>(); ++i ){
    status = FitPDFToData( pdf, data, args );
    if( status->status() != 0 ){
      if( !og_args.Has( "Save" ) ){
        delete status;
        status = nullptr;
      }
      return status;
    }
    delete status;
  }

  if( !og_args.Has( "Save" ) ){
    delete status;
    status = nullptr;
  }
  return status;
}

}/* namespace usr */
