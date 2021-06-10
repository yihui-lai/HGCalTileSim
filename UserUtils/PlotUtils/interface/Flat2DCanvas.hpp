#ifndef USERUTILS_PLOTUTILS_FLAT2DCANVAS_HPP
#define USERUTILS_PLOTUTILS_FLAT2DCANVAS_HPP

#ifdef CMSSW_GIT_HASH
#include "UserUtils/PlotUtils/interface/Canvas.hpp"
#include "UserUtils/PlotUtils/interface/Pad2DFlat.hpp"
#else
#include "UserUtils/PlotUtils/Canvas.hpp"
#include "UserUtils/PlotUtils/Pad2DFlat.hpp"
#endif

namespace usr {

namespace plt {

/**
 * @brief Canvas containing a single 2D pad
 */
class Flat2DCanvas : public Canvas
{
public:
  static length_t default_width;
  static length_t default_height;
  static FontSet default_font;

  Flat2DCanvas(
    const length_t width  = default_width,
    const length_t height = default_height,
    const FontSet&        = default_font
    );

  virtual
  ~Flat2DCanvas ();

  /**
   * @brief Returning reference to the underlying Pad2DFlat instance
   */
  inline Pad2DFlat&
  Pad(){ return GetPad<Pad2DFlat>( 0 ); }

#define CANVAS2D_PASSTHROUGH_FUNC( FUNC_NAME, PLOT_TYPE, RET_TYPE ) \
  template<typename ... Args>                                       \
  inline RET_TYPE& FUNC_NAME( PLOT_TYPE& obj,    Args ... args )    \
  { return this->Pad().FUNC_NAME( obj, args ... ); \
  }                \
  template<typename ... Args>                                       \
  inline RET_TYPE& FUNC_NAME( PLOT_TYPE* obj,    Args ... args )    \
  { return this->Pad().FUNC_NAME( obj, args ... ); }                \

  /**
   * @{
   * @brief passing through for Pad2DFlat plotting function so that user can
   * 'plot' on this canvas type.
   * @details See Pad2DFlats plot functions for details.
   */
  CANVAS2D_PASSTHROUGH_FUNC( PlotHist,    TH2D,     TH2D     );
  CANVAS2D_PASSTHROUGH_FUNC( Plot1DGraph, TGraph,   TGraph   );
  CANVAS2D_PASSTHROUGH_FUNC( PlotGraph,   TGraph2D, TGraph2D );
  CANVAS2D_PASSTHROUGH_FUNC( PlotFunc,    TF2,      TGraph2D );
  /** @} */
#undef CANVAS2D_PASSTHROUGH_FUNC

#define CANVAS2D_ACCESS_PASSTHROUGH( FUNC_NAME, RET_TYPE )              \
  inline RET_TYPE FUNC_NAME(){ return Pad().FUNC_NAME(); }             \

  /**
   * @{
   * @brief Pass through for axis access
   */
  CANVAS2D_ACCESS_PASSTHROUGH( Xaxis, TAxis& );
  CANVAS2D_ACCESS_PASSTHROUGH( Yaxis, TAxis& );
  CANVAS2D_ACCESS_PASSTHROUGH( Zaxis, TAxis& );
  /** @} */

#undef CANVAS2D_ACCESS_PASSTHROUGH

#define CANVAS2D_VOID_PASSTHROUGH( FUNC_NAME ) \
  template<typename ... Args>            \
  inline void FUNC_NAME( Args ... args ) \
  { this->Pad().FUNC_NAME( args ... ); }

  /**
   * @{
   * @brief passing through addition Pad2DFlat functions
   */
  CANVAS2D_VOID_PASSTHROUGH( DrawCMSLabel );
  CANVAS2D_VOID_PASSTHROUGH( DrawLuminosity );
  /** @} */
#undef CANVAS2D_VOID_PASSTHROUGH

  inline void
  SetLogx( int x = 1  ){ Pad().SetLogx( x ); }
  inline void
  SetLogy( int x = 1  ){ Pad().SetLogy( x ); }
  inline void
  SetLogz( int x = 1  ){ Pad().SetLogz( x ); }

protected:
  void init_margin();
};

}

}

#endif
