#ifdef CMSSW_GIT_HASH
#include "UserUtils/PlotUtils/interface/Pad2DFlat.hpp"
#else
#include "UserUtils/PlotUtils/Pad2DFlat.hpp"
#endif

namespace usr {

namespace plt {

/**
 * @brief Getting the object that is responsible for plotting the axis.
 *
 * Assuming the first object that has axis settings (This should be true for
 * typical use cases)
 */
TObject*
Pad2DFlat::GetAxisObject() const
{
  if( !GetListOfPrimitives() ){ return nullptr; }

  for( const auto&& obj : *GetListOfPrimitives() ){
    if( obj->InheritsFrom( TH2::Class() )
        || obj->InheritsFrom( TGraph2D::Class() ) ){
      return obj;
    }
  }

  return nullptr;
}

/**
 * @brief common settings for setting the axis fonts
 */
void
Pad2DFlat::SetAxisFont()
{
  Xaxis().SetLabelFont( FontFace() );
  Xaxis().SetTitleFont( FontFace() );
  Xaxis().SetTitleSize( FontSize() );
  Xaxis().SetLabelSize( Font().small() );
  Xaxis().SetTitleOffset( 1.2 );

  Yaxis().SetLabelFont( FontFace() );
  Yaxis().SetTitleFont( FontFace() );
  Yaxis().SetTitleSize( FontSize() );
  Yaxis().SetLabelSize( Font().small() );
  Yaxis().SetTitleOffset( 1.2 );

  Zaxis().SetLabelFont( FontFace() );
  Zaxis().SetTitleFont( FontFace() );
  Zaxis().SetTitleSize( FontSize() );
  Zaxis().SetLabelSize( Font().small() );
  Zaxis().SetTitleOffset( 1.5 );
}

#define AXISOBJ_ACCESS( ACTION, DEFAULT )                  \
  TObject* axisobj = GetAxisObject();                      \
  if( !axisobj ){                                          \
    return DEFAULT;                                        \
  } else if( axisobj->InheritsFrom( TH2::Class() ) ){      \
    return dynamic_cast<TH1*>( axisobj )->ACTION;          \
  } else if( axisobj->InheritsFrom( TGraph2D::Class() ) ){ \
    return dynamic_cast<TGraph2D*>( axisobj )->ACTION;     \
  } else {                                                 \
    return DEFAULT;                                        \
  }

#define AXISOBJ_ADJUST( ACTION )                          \
  TObject* axisobj = GetAxisObject();                     \
  if( !axisobj ){                                         \
    return;                                               \
  } else if( axisobj->InheritsFrom( TH2::Class() ) ){     \
    return dynamic_cast<TH1*>( axisobj )->ACTION;         \
  } else if( axisobj->InheritsFrom( TGraph2D::Class() ) ){ \
    return dynamic_cast<TGraph*>( axisobj )->ACTION;      \
  }

/**
 * @{
 * @brief returning pointer to axis object.
 * @details
 * ROOT style pointer interface for getting the axis object used for plotting
 * the axis frame.
 */
TAxis*
Pad2DFlat::GetXaxis() const { AXISOBJ_ACCESS( GetXaxis(), nullptr ); }

TAxis*
Pad2DFlat::GetYaxis() const { AXISOBJ_ACCESS( GetYaxis(), nullptr ); }

TAxis*
Pad2DFlat::GetZaxis() const { AXISOBJ_ACCESS( GetZaxis(), nullptr ); }
/** @} */

/**
 * @{
 * @brief   reference interface for getting axis (No null pointer handling)
 */
TAxis&
Pad2DFlat::Xaxis(){ return *GetXaxis(); }

TAxis&
Pad2DFlat::Yaxis(){ return *GetYaxis(); }

TAxis&
Pad2DFlat::Zaxis(){ return *GetZaxis(); }

const TAxis&
Pad2DFlat::Xaxis() const { return *GetXaxis(); }

const TAxis&
Pad2DFlat::Yaxis() const { return *GetYaxis(); }

const TAxis&
Pad2DFlat::Zaxis() const { return *GetZaxis(); }
/** @} */

}

}// usr
