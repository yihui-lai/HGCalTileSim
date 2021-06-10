/**
 * @file    Pad1D.hpp
 * @brief   A pad designed for 1D data plotting.
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */
#ifndef USERUTILS_PLOTUTILS_PAD1D_HPP
#define USERUTILS_PLOTUTILS_PAD1D_HPP

#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/STLUtils/VectorUtils.hpp"
#include "UserUtils/Common/interface/RootUtils/RooArgContainer.hpp"
#include "UserUtils/PlotUtils/interface/Canvas.hpp"
#include "UserUtils/PlotUtils/interface/Constants.hpp"
#include "UserUtils/PlotUtils/interface/PlotCommon.hpp"
#else
#include "UserUtils/Common/STLUtils/VectorUtils.hpp"
#include "UserUtils/Common/RootUtils/RooArgContainer.hpp"
#include "UserUtils/PlotUtils/Canvas.hpp"
#include "UserUtils/PlotUtils/Constants.hpp"
#include "UserUtils/PlotUtils/PlotCommon.hpp"
#endif

#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooCmdArg.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooRealVar.h"

#include "TF1.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLine.h"
#include "TProfile.h"


#include <vector>

namespace usr  {

namespace plt {

/**
 * @brief Specialized class for 1D data plotting.
 */
class Pad1D : public PadBase
{
protected:
  Pad1D( const PadSize&, const RangeByVar& );
  Pad1D( const PadSize&, const double min, const double max );
  Pad1D( const PadSize& );

  virtual void InitDraw();
  virtual void Finalize();

  friend class Canvas;

public:
  virtual ~Pad1D ();
  Pad1D()               = delete;
  Pad1D( const Pad1D& ) = delete;

#define DECLARE_PLOT_FUNCTIONS( FUNC_NAME, TYPE, RET_TYPE )                    \
  RET_TYPE FUNC_NAME( TYPE&, const std::vector<RooCmdArg> & );                 \
  inline RET_TYPE FUNC_NAME( TYPE* x, const std::vector<RooCmdArg> &list ){    \
    return FUNC_NAME( *x, list );                                              \
  }                                                                            \
  inline RET_TYPE FUNC_NAME( TYPE& x ){ return FUNC_NAME( x, {} ); }           \
  inline RET_TYPE FUNC_NAME( TYPE* x ){ return FUNC_NAME( x, {} ); }           \
  template<typename ... Args>                                                  \
  inline RET_TYPE FUNC_NAME( TYPE& x, const RooCmdArg & arg1, Args ... args ){ \
    return FUNC_NAME( x, MakeVector<RooCmdArg>( arg1, args ... ) );            \
  }                                                                            \
  template<typename ... Args>                                                  \
  inline RET_TYPE FUNC_NAME( TYPE* x, const RooCmdArg & arg1, Args ... args ){ \
    return FUNC_NAME( x, MakeVector<RooCmdArg>( arg1, args ... ) );            \
  }

  /**
   * @{
   * @brief Plotting 1D histogram objects
   */
  DECLARE_PLOT_FUNCTIONS( PlotHist, TH1D, TH1D& );
  /** @} */

  /**
   * @{
   * @brief Special case for TProfile plotting
   */
  DECLARE_PLOT_FUNCTIONS( PlotProfile, TProfile, TH1D& );
  /** @} */

  /**
   * @{
   * @brief Plotting x-y scatter graph objects
   */
  DECLARE_PLOT_FUNCTIONS( PlotGraph, TGraph, TGraph& );
  /** @} */

  /**
   * @{
   * @brief Plotting ROOT flavoured 1D function
   */
  DECLARE_PLOT_FUNCTIONS( PlotFunc, TF1, TGraph& );
  /** @} */

  /**
   * @{
   * @brief Plotting RooFit's @ROOT{RooAbsData}s
   */
  DECLARE_PLOT_FUNCTIONS( PlotData, RooAbsData, TGraphAsymmErrors& );
  /** @} */

  /**
   * @{
   * @brief
   * Plotting RooFit's @ROOT{RooAbsPdf}s
   */
  DECLARE_PLOT_FUNCTIONS( PlotPdf, RooAbsPdf, TGraph& );
  /** @} */


  /**
   * @{
   *  @brief function for drawing a horizontal line across pad
   */
  DECLARE_PLOT_FUNCTIONS( DrawHLine, const double, TLine& );
  /** @} */

  /**
   * @{
   * @brief Function for drawing a vertical along in the pad
   */
  DECLARE_PLOT_FUNCTIONS( DrawVLine, const double, TLine& );
  /** @} */


#undef DECLARE_PLOT_FUNCTIONS

  void DrawCMSLabel( const std::string&      = cap::prelim,
                     const std::string& main = "CMS" );
  void DrawLuminosity( const double luminosity );
  void DrawLuminosity( const std::string& string );

  // Common text settings
  float InnerTextLeft() const;
  float InnerTextRight() const;
  float InnerTextHCenter() const;
  float InnerTextTop() const;
  float InnerTextBottom() const;
  float InnerTextVCenter() const;

  // Axis and Range setting function
  TObject*     GetAxisObject() const;
  virtual void SetAxisFont();

  TAxis&       Xaxis();
  TAxis&       Yaxis();
  const TAxis& Xaxis() const;
  const TAxis& Yaxis() const;
  TAxis*       GetXaxis() const;
  TAxis*       GetYaxis() const;

  double GetDataMax() const;
  double GetDataMin() const;
  double GetXaxisMin() const;
  double GetXaxisMax() const;
  double GetYaxisMax() const;
  double GetYaxisMin() const;
  void   SetYaxisMax( const double );
  void   SetYaxisMin( const double );
  void   SetDataMax( const double );
  void   SetDataMin( const double );

  // Histogram axis title setting
  static const double autobinwidth;
  static const double forcevarbinwidth;
  void SetHistAxisTitles(
    const std::string& title,
    const std::string& unit          = "",
    const double       forcebinwidth = autobinwidth
    );

  /**
   * @brief enum for specifying y range adjustment type.
   */
  enum rangetype { hist, graph, ratio, pull, aut };
  void AutoSetYRange( const rangetype = aut );
  void SetLogy( int );
  void SetLogx( int );

  void AddLegendEntry( TObject&,
                       const std::string& title,
                       const std::string& format );

  void FinalizeLegend( const align x );

  /** @brief returning reference to internal legend object. */
  inline TLegend&
  Legend(){ return _legend; }

  /** @brief returning reference to internal Roo1DFrame object. */
  inline RooFrame&
  FrameObj(){ return _frame; }

  /** @brief returning reference to internal Roo1DFrame object. */
  inline const RooFrame&
  FrameObj() const { return _frame; }

  /** @brief returning the stored y range adjustment type */
  inline rangetype&
  RangeType(){ return _prevrangetype; }

protected:

  /**
   * @brief RooPlot object used for generating graphs for ROOFIT objects
   */
  RooFrame _frame;

  /**
   * @brief Internal @ROOT{THStack} object for assisting with stacked histogram
   * plots.
   */
  THStack* _workingstack;

  /** @{ @brief  tracking the min/max values of the data plotted on the pad  */
  double _datamax;
  double _datamin;
  /** @} */

  /** @brief
   * The range type used to adjust the Y axis range on the previous adjustment
   * call
   */
  rangetype _prevrangetype;

  /** @brief the Legend object */
  TLegend _legend;

  /** @brief position of the legend object in the pad */
  align _legendposition;

  /** @brief list of lines that need finalizing just before plot is finalized */
  std::vector<TLine*> _linelist;

  static
  void SetAxisTitle(
    TAxis&,
    const std::string& title,
    const std::string& unit = "" );

  void AutoSetYRangeHist();
  void AutoSetYRangeGraph();
  void AutoSetYRangeRatio();
  void AutoSetYRangePull();
  void FixVLines();

  void _init_legend();
  void AddLegendEntry( TH1D&, const std::string&, const RooCmdArg& );
  void AddLegendEntry( TGraph&, const std::string&, const RooCmdArg& );

  // Helper function for Plot<> Functions
  void    TrackObjectY( const TObject& obj, const int tracky );
  TGraph& MakePdfGraph( RooAbsPdf&             pdf,
                        const RooArgContainer& arglist );
  TGraph& MakeTF1Graph( TF1&,
                        const RooArgContainer& arglist );
  TGraphAsymmErrors& MakeDataGraph( RooAbsData&            data,
                                    const RooArgContainer& arglist );
  TGraphAsymmErrors& GenGraph( RooAbsData& data, RooLinkedList& arglist );
  TGraph&            GenGraph( RooAbsPdf& pdf, RooLinkedList& arglist );

  bool CheckLogy() const;
  bool CheckLogx() const;
};

}/* plt */

}/* usr  */

#endif/* end of include guard: USERUTILS_PLOTUTILS_PADFRAME1D_HPP */
