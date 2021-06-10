/**
 * @file Pad2DFlat.hpp
 * @author Yi-Mu "Enoch" Chen
 */
#ifndef USERUTILS_PLOTUTILS_PAD2DFLAT_HPP
#define USERUTILS_PLOTUTILS_PAD2DFLAT_HPP

#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/STLUtils/VectorUtils.hpp"
#include "UserUtils/PlotUtils/interface/Canvas.hpp"
#include "UserUtils/PlotUtils/interface/Constants.hpp"
#include "UserUtils/PlotUtils/interface/PlotCommon.hpp"
#else
#include "UserUtils/Common/STLUtils/VectorUtils.hpp"
#include "UserUtils/PlotUtils/Canvas.hpp"
#include "UserUtils/PlotUtils/Constants.hpp"
#include "UserUtils/PlotUtils/PlotCommon.hpp"
#endif

#include "RooPlot.h"

#include "TH2D.h"
#include "TF2.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TLegend.h"

namespace usr {

namespace plt {

/**
 * @brief Specialized class for 2D data object representation in x-y plane
 * projections
 */
class Pad2DFlat : public PadBase
{
protected:
  Pad2DFlat( const PadSize& );
  Pad2DFlat( const PadSize&,
             const double xmin, const double xmax,
             const double ymin, const double ymax );

  virtual void InitDraw();
  virtual void Finalize();

  friend class Canvas;

public:
  virtual ~Pad2DFlat();
  Pad2DFlat()                   = delete;
  Pad2DFlat( const Pad2DFlat& ) = delete;

#define DECLARE_PLOT_FUNCTIONS( FUNC_NAME, TYPE, RET_TYPE ) \
  RET_TYPE  FUNC_NAME( TYPE&, const std::vector<RooCmdArg> & );                \
  inline RET_TYPE  FUNC_NAME( TYPE* x, const std::vector<RooCmdArg> &list ){   \
    return FUNC_NAME( *x, list );                                              \
  }                                                                            \
  inline RET_TYPE  FUNC_NAME( TYPE& x ){ return FUNC_NAME( x, {} ); }          \
  inline RET_TYPE  FUNC_NAME( TYPE* x ){ return FUNC_NAME( x, {} ); }          \
  template<typename ... Args>                                                  \
  inline RET_TYPE  FUNC_NAME( TYPE& x, const RooCmdArg & arg1, Args ... args ){\
    return FUNC_NAME( x, MakeVector<RooCmdArg>( arg1, args ... ) );            \
  }                                                                            \
  template<typename ... Args>                                                  \
  inline RET_TYPE  FUNC_NAME( TYPE* x, const RooCmdArg & arg1, Args ... args ){\
    return FUNC_NAME( x, MakeVector<RooCmdArg>( arg1, args ... ) );            \
  }

  /**
   * @{
   * @brief Plotting 2D histogram object
   */
  DECLARE_PLOT_FUNCTIONS( PlotHist, TH2D, TH2D& );
  /** @} */

  /**
   * @{
   * @brief Plotting 2D function object
   */
  DECLARE_PLOT_FUNCTIONS( PlotFunc, TF2, TGraph2D& );
  /** @} */

  /**
   * @{
   * @brief Plotting 2D graph object
   */
  DECLARE_PLOT_FUNCTIONS( PlotGraph, TGraph2D, TGraph2D& );
  /** @} */

  /**
   * @{
   * @brief Plotting 1D graph object
   */
  DECLARE_PLOT_FUNCTIONS( Plot1DGraph, TGraph, TGraph& );
  /** @} */

#undef DECLARE_PLOTFUNCTIONS

  TObject* GetAxisObject() const ;
  TAxis&       Xaxis();
  TAxis&       Yaxis();
  TAxis&       Zaxis();
  const TAxis& Xaxis() const;
  const TAxis& Yaxis() const;
  const TAxis& Zaxis() const;
  TAxis*       GetXaxis() const;
  TAxis*       GetYaxis() const;
  TAxis*       GetZaxis() const;

  void SetAxisFont();

  /**
   * @brief returning reference to internal legend object.
   */
  inline TLegend& Legend(){ return _legend; }


  void DrawCMSLabel( const std::string&      = cap::prelim,
                     const std::string& main = "CMS" );
  void DrawLuminosity( const double luminosity );
  void DrawLuminosity( const std::string& string );

  inline void SetLogx( int x = 1  ){ TPad::SetLogx(x) ; }
  inline void SetLogy( int x = 1  ){ TPad::SetLogy(x) ; }
  inline void SetLogz( int x = 1  ){ TPad::SetLogz(x) ; }

private:

  TLegend _legend;
  RooFrame _frame;
  void init_legend();
  void MakeLegend();
};




}

}



#endif
