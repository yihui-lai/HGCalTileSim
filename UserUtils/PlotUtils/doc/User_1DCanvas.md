@defgroup PlotUtilsUser PlotUtilsUser
@ingroup  PlotUtils
@brief    User level classes for fast plotting with decent outputs.

While it is unambiguous to defined a rigid pad--canvas relation, it could be
inconvenient to the user who just wants to write a plot an figure.
This package provides classes that provides default layouts that are commonly
used in CMS analysis, as well as providing interface short hands that might be
more familiar to ROOT uses (e.x. plotting on a canvas).

@class usr::plt::Simple1DCanvas
@ingroup PlotUtilsUser
@brief Canvas with exactly 1 Pad1D instance

The Simple1DCanvas class, is a canvas containing exactly 1 Pad1D instance. Its
member functions are generated such that it looks like a users plotting calls
are directly performed on the canvas, greatly simplifying the plotting process.

This class is designed as a preset for simple data plotting, with the margin,
dimension and font defaults tuned to produce decent, publication ready plots
with minimal amount of user code.

The following are a few examples of plotting using the Simple1DCanvas (with the
data filling code removed)

# Plotting style with individual histograms

```cpp
    plt::Simple1DCanvas c;
    c.PlotHist( h1,
      plt::PlotType( plt::scatter ),
      plt::EntryText( "Bkg. 1" ),
      plt::LineColor( usr::plt::col::black ),
      plt::MarkerStyle( usr::plt::sty::mkrcircle ),
      plt::MarkerColor( usr::plt::col::black ),
      plt::MarkerSize( 0.2 ) );
    c.PlotHist( h2,
      plt::PlotType( plt::hist ),
      plt::EntryText( "Bkg. 2" ),
      plt::LineColor( usr::plt::col::blue ) );
    c.PlotHist( h3,
      plt::PlotType( plt::histerr ),
      plt::EntryText( "Signal" ),
      plt::FillColor( usr::plt::col::darkgray ),
      plt::FillStyle( plt::sty::filldotdense ) );

    c.DrawCMSLabel( "ROOTOBJ TEST", "CWS" );
    c.DrawLuminosity( 133.7 );
    c.SetHistAxisTitles( "P_{t}", plt::unit::GeVc );

    c.SetLogy( 1 );
    c.SaveAsPNG( "image/simple1dcanvas_log.png", 72 );
```

<div class="plot_example">
<img src="image/simple1dcanvas_log.png"/>

<ul>
  <li><a href="image/simple1dcanvas_log.png">72dpi output</a></li>
  <li><a href="image/simple1dcanvas_log_highres.png">300dpi output</a></li>
  <li><a href="image/simple1dcanvas_log.pdf">PDF output</a></li>
  </ul>
</div>

# Plotting with stacked histograms

```cpp
    plt::Simple1DCanvas c;
    c.PlotHist( h1,
      plt::PlotType( plt::histstack ),
      plt::EntryText( "Bkg. 1" ),
      plt::FillStyle( plt::sty::fillsolid ),
      plt::FillColor( plt::col::cyan ),
      plt::LineColor( plt::col::cyan ) );
    c.PlotHist( h2,
      plt::PlotType( plt::histstack ),
      plt::EntryText( "Bkg. 2" ),
      plt::FillStyle( plt::sty::filldotdense ),
      plt::FillColor( plt::col::pink ),
      plt::LineColor( plt::col::pink ) );
    c.PlotHist( h3,
      plt::PlotType( plt::hist ),
      plt::EntryText( "Sig" ),
      plt::LineColor( usr::plt::col::blue ) );

    c.DrawCMSLabel( "STACK TEST", "CWS" );
    c.DrawLuminosity( 133.7 );
    c.SetHistAxisTitles( "P_{t}", plt::unit::GeVc );

    c.SaveAsPNG( "image/simple1dcanvas_stack.png",         72 );
```

<div class="plot_example">
<img src="image/simple1dcanvas_stack.png"/>

<ul>
  <li><a href="image/simple1dcanvas_stack.png">72dpi output</a></li>
  <li><a href="image/simple1dcanvas_stack_highres.png">300dpi output</a></li>
  <li><a href="image/simple1dcanvas_stack.pdf">PDF output</a></li>
  </ul>
</div>

# Plotting with RooFit

```cpp
    plt::Simple1DCanvas c2( x );

    c2.PlotData( d,
      plt::EntryText("Data?"),
      plt::MarkerStyle( plt::sty::mkrcircle ),
      plt::MarkerSize( 0.5 ),
      plt::MarkerColor( plt::col::black ),
      plt::LineColor( plt::col::black ) );

    c2.PlotPdf( g,
      RooFit::VisualizeError( *fit, 1, false ),
      plt::EntryText( "My Fit" ),
      plt::PlotUnder( dgraph ),
      plt::LineColor( plt::col::blue ),
      plt::FillColor( plt::col::cyan )
      );

    c2.DrawCMSLabel( "ROOFIT test", "CWS" );
    c2.DrawLuminosity( 133.7 );
    c2.SetHistAxisTitles( "M_{HHHHH}", plt::unit::GeVcc );

    c2.SaveAsPNG( "image/simple1dcanvas_roofit.png",         72 );
```

<div class="plot_example">
<img src="image/simple1dcanvas_roofit.png"/>

<ul>
  <li><a href="image/simple1dcanvas_roofit.png">72dpi output</a></li>
  <li><a href="image/simple1dcanvas_roofit_highres.png">300dpi output</a></li>
  <li><a href="image/simple1dcanvas_roofit.pdf">PDF output</a></li>
  </ul>
</div>

@class   usr::plt::Ratio1DCanvas
@ingroup PlotUtilsUser
@brief   Canvas with exactly 2 Pad1D instances for comparison plots.

A canvas where the top pad displays the data, and the bottom pad is used to
display the ratio/pull between data sets is a common method of presenting
numbers in HEP. The ratio canvas, apart from setting up the pads for such
plotting, also provides tools to generate the required plots.

The `PlotObj` functions of the top pad is pass to the canvas, while additional
functions are used for plotting the ratio/pull plots onto the bottom pad. The
generated ratio/pull plot will be owned by the canvas and be destroyed when the
canvas goes out of scope.

The following is some examples of plots created with the Ratio1DCanvas class:

# Stack plots with error and signal comparison

```cpp
    plt::Ratio1DCanvas c;

    c.PlotHist( hist3,
      plt::PlotType( plt::histstack ),
      plt::EntryText( "Bkg_{1}: VV" ),
      plt::FillColor( plt::col::blue ),
      plt::LineColor( plt::col::blue ) );
    c.PlotHist( hist2,
      plt::PlotType( plt::histstack ),
      plt::EntryText( "Bkg_{2}: HH" ),
      plt::FillColor( plt::col::red ),
      plt::LineColor( plt::col::red ) );
    c.PlotHist( hist1,
      plt::PlotType( plt::histstack ),
      plt::EntryText( "Bkg_{3}: TT" ),
      plt::FillColor( plt::col::cyan ),
      plt::LineColor( plt::col::cyan ) );
    c.PlotHist( histsum,
      plt::PlotType( plt::histerr ),
      plt::EntryText( "Bkg_{1} Unc." ),
      plt::FillStyle( plt::sty::filldotdense ),
      plt::FillColor( plt::col::darkgray ) );
    c.PlotHist( data,
      plt::PlotType( plt::scatter ),
      plt::EntryText( "Data" ),
      plt::MarkerStyle( 20 ),
      plt::MarkerSize( 0.2 ),
      plt::LineColor( plt::col::black ) );

    c.PlotScale( histsum, histsum,
      plt::PlotType( plt::histerr ),
      plt::TrackY( plt::TrackY::none ) );
    c.PlotScale( data, histsum,
      plt::PlotType( plt::scatter ) );

    c.TopPad().SetHistAxisTitles( "M_{THTHTH}", plt::unit::GeVcc );
    c.BottomPad().SetHistAxisTitles( "M_{THTHTH}", plt::unit::GeVcc );
    c.BottomPad().Yaxis().SetTitle( "Data/Sim." );

    c.DrawCMSLabel( "Ratio1DCanvas", "CWS" );
    c.DrawLuminosity( 133.7 );
    c.SaveAsPNG(         "image/ratio1dcanvas.png",  72 );
```

<div class="plot_example">
<img src="image/ratio1dcanvas.png"/>

<ul>
  <li><a href="image/ratio1dcanvas.png">72dpi output</a></li>
  <li><a href="image/ratio1dcanvas_highres.png">300dpi output</a></li>
  <li><a href="image/ratio1dcanvas.pdf">PDF output</a></li>
  </ul>
</div>


# Plotting with RooFit and fitting uncertainties

```cpp
    plt::Ratio1DCanvas c( x );
    auto& datgraph = c.PlotData( d,
      plt::EntryText( "Fake Data" ),
      plt::MarkerStyle( 20 ),
      plt::MarkerSize( 0.2 ),
      plt::LineColor( kBlack ) );

    auto& altgraph = c.PlotPdf( gx,
      plt::EntryText( "Alt. Bkg. Model" ),
      plt::PlotUnder( datgraph ),
      plt::LineColor( plt::col::red ),
      plt::FillColor( plt::col::pink ) );

    auto& fitgraph = c.PlotPdf( g,
      plt::VisualizeError( *fit, 1, false ),
      plt::EntryText( "Fitted Bkg." ),
      plt::PlotUnder( datgraph ),
      plt::LineColor( usr::plt::col::blue ),
      plt::FillColor( usr::plt::col::cyan )  );

    c.PlotScale( fitgraph, fitgraph,
      plt::TrackY( plt::TrackY::none ), plt::PlotType( plt::fittedfunc ) );
    c.PlotScale( altgraph, fitgraph,
      plt::TrackY( plt::TrackY::none ), plt::PlotType( plt::simplefunc ) );
    c.PlotScale( datgraph, fitgraph,
      plt::PlotType( plt::scatter ) );

    c.BottomPad().Yaxis().SetTitle( "Data/Fit" );

    c.DrawCMSLabel( "Ratio1DCanvas", "CWS" );
    c.DrawLuminosity( 133.7 );
    c.SaveAsPNG( "image/ratio1dcanvas_roofit.png",  72 );
```

<div class="plot_example">
<img src="image/ratio1dcanvas_roofit.png"/>

<ul>
  <li><a href="image/ratio1dcanvas_roofit.png">72dpi output</a></li>
  <li><a href="image/ratio1dcanvas_roofit_highres.png">300dpi output</a></li>
  <li><a href="image/ratio1dcanvas_roofit.pdf">PDF output</a></li>
  </ul>
</div>
