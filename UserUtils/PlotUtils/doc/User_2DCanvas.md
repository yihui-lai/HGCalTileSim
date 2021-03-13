@class usr::plt::Flat2DCanvas
@ingroup PlotUtilsUser

The Flat2DCanvas class is a canvas containing exactly Pad2DFlat instance. Its
member functions are generated such that it looks like a users plotting calls
are directly performed on the canvas, greatly simplifying the plotting process.

The pad assumes that the plotted objects would be fill up the main axis frame,
hence the choice of the CMS labels and luminosity labels being placed both
outside the main axis frame. The legend and other text, should therefore be
placed on the left of the main canvas.

The following are a few examples of plotting using the Flat2DCanvas (with
the style/color settings and data filling code removed)

```cpp
usr::plt::Flat2DCanvas c;

auto& gh = c.PlotHist( hist,
  usr::plt::Plot2DF( usr::plt::heat ) );
auto& gf1 = c.PlotFunc( xyg,
  usr::plt::Plot2DF( usr::plt::cont ),
  usr::plt::EntryText( "Model 1" )
  );
auto& gf2 = c.PlotFunc( xyg2,
  usr::plt::Plot2DF( usr::plt::cont ),
  usr::plt::EntryText( "Model 2" )
  );

c.DrawCMSLabel("2D plotting", "CWS" );
c.DrawLuminosity( 133.7 );

c.Xaxis().SetTitle( "My X value" );
c.Yaxis().SetTitle( "My y value" );
c.Zaxis().SetTitle( "Readout" );

c.SaveAsPNG( "image/flat2dcanvas.png" , 72 );
c.SaveAsPNG( "image/flat2dcanvas_highres.png" , 300 );
c.SaveAsPDF( "image/flat2dcanvas.pdf" );
```

@image html flat2dcanvas.png
