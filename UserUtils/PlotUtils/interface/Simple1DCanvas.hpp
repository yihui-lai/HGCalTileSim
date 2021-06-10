/**
 * @file    Simple1DCanvas.hpp
 * @brief   A canvas with single 1D pad. Good for fast plotting!
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */
#ifndef USERUTILS_PLOTUTILS_SIMPLE1DCANVAS_HPP
#define USERUTILS_PLOTUTILS_SIMPLE1DCANVAS_HPP

#ifdef CMSSW_GIT_HASH
#include "UserUtils/PlotUtils/interface/Canvas.hpp"
#include "UserUtils/PlotUtils/interface/Pad1D.hpp"
#else
#include "UserUtils/PlotUtils/Canvas.hpp"
#include "UserUtils/PlotUtils/Pad1D.hpp"
#endif

namespace usr  {

namespace plt  {

class Simple1DCanvas : public Canvas
{
public:
  static length_t default_width;
  static length_t default_height;
  static FontSet  default_font;

  Simple1DCanvas(
    const RangeByVar&,
    const length_t width  = default_width,
    const length_t height = default_height,
    const FontSet& = default_font
    );

  Simple1DCanvas(
    const length_t width  = default_width,
    const length_t height = default_height,
    const FontSet& = default_font
    );

  virtual
  ~Simple1DCanvas ();

  /**
   * @brief Returning reference to the created Pad1D.
   */
  inline Pad1D&
  Pad(){ return GetPad<Pad1D>( 0 ); }

#define PASSTHROUGH_PLOTFUNC( FUNC_NAME, PLOT_TYPE, RET_TYPE ) \
  template<typename ... Args>                                  \
  inline RET_TYPE &FUNC_NAME( PLOT_TYPE& obj, Args ... args )  \
  { return this->Pad().FUNC_NAME( obj, args ... ); }           \
  template<typename ... Args>                                  \
  inline RET_TYPE &FUNC_NAME( PLOT_TYPE* obj, Args ... args )  \
  { return this->Pad().FUNC_NAME( obj, args ... );}

  /**
   * @{
   * @brief passing through for Pad1D plotting function so that user can 'plot'
   * on this canvas type.
   */
  PASSTHROUGH_PLOTFUNC( PlotHist,  TH1D,       TH1D );
  PASSTHROUGH_PLOTFUNC( PlotFunc,  TF1,        TGraph );
  PASSTHROUGH_PLOTFUNC( PlotGraph, TGraph,     TGraph );
  PASSTHROUGH_PLOTFUNC( PlotData,  RooAbsData, TGraphAsymmErrors );
  PASSTHROUGH_PLOTFUNC( PlotPdf,   RooAbsPdf,  TGraph );
  /** @} */
#undef PASSTHROUGH_PLOTFUNC

#define PASSTHROUGH_VOIDFUNC( FUNC_NAME ) \
  template<typename ... Args>            \
  inline void FUNC_NAME( Args ... args ) \
  { this->Pad().FUNC_NAME( args ... ); }

  /**
   * @{
   * @brief passing through addition Pad1D functions
   */
  PASSTHROUGH_VOIDFUNC( DrawCMSLabel );
  PASSTHROUGH_VOIDFUNC( DrawLuminosity );
  PASSTHROUGH_VOIDFUNC( DrawHLine );
  PASSTHROUGH_VOIDFUNC( DrawVLine );
  PASSTHROUGH_VOIDFUNC( SetLogy );
  PASSTHROUGH_VOIDFUNC( SetHistAxisTitles );
  /** @} */
#undef PASSTHROUGH_VOIDFUNC

#define PASSTHROUGH_ACCESSFUNC( FUNC_NAME, RET_TYPE )\
  inline RET_TYPE FUNC_NAME() { return Pad().FUNC_NAME(); }

  /**
   * @{
   * @brief Pass through for axis access
   */
  PASSTHROUGH_ACCESSFUNC( Xaxis, TAxis& );
  PASSTHROUGH_ACCESSFUNC( Yaxis, TAxis& );
  /** @} */

#undef PASSTHROUGH_ACCESSFUNC

protected:
  void _init_margin();
};

}/* plt  */

}/* usr  */


#endif/* end of include guard: USERUTILS_PLOTUTILS_SIMPLE1DCANVAS_HPP  */
