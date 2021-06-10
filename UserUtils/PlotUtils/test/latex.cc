/**
 * @file    latex.cc
 * @brief   Testing the Latex dimesion estimation class. Also responsible for
 *          generating the character width table in PlotCommon_Latex.cc
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */
#ifdef CMSSW_GIT_HASH
#include "UserUtils/PlotUtils/interface/PlotCommon.hpp"
#else
#include "UserUtils/PlotUtils/PlotCommon.hpp"
#endif

#include "TCanvas.h"
#include "TLatex.h"
#include <map>
#include <regex>

int
main( int argc, char const* argv[] )
{
  std::cout << "Generating width and height map" << std::endl;
  const double fontsize      = 20;
  const unsigned mult        = 1;
  const unsigned canvas_size = 1000;
  TCanvas c( "c", "c", canvas_size, canvas_size );
  TLatex t;
  t.SetTextFont( 43 );
  t.SetTextSize( fontsize );

  std::map<char, double> charwidth;
  std::map<char, double> charheight;

  for( char i = 'a'; i <= 'z'; ++i ){
    const std::string s( mult, i );
    charwidth[i]  = t.DrawLatexNDC( 0.5, 0.5, s.c_str() )->GetXsize() * canvas_size;
    charheight[i] = t.DrawLatexNDC( 0.5, 0.5, s.c_str() )->GetYsize() * canvas_size;
  }

  for( char i = 'A'; i <= 'Z'; ++i ){
    const std::string s( mult, i );
    charwidth[i]  = t.DrawLatexNDC( 0.5, 0.5, s.c_str() )->GetXsize() * canvas_size;
    charheight[i] = t.DrawLatexNDC( 0.5, 0.5, s.c_str() )->GetYsize() * canvas_size;
  }

  for( char i = '0'; i <= '9'; ++i ){
    const std::string s( mult, i );
    charwidth[i]  = t.DrawLatexNDC( 0.5, 0.5, s.c_str() )->GetXsize() * canvas_size;
    charheight[i] = t.DrawLatexNDC( 0.5, 0.5, s.c_str() )->GetYsize() * canvas_size;
  }

  for( const char i : " ,./<>!@$%%&*()[]-=+;" ){
    const std::string s( mult, i );
    charwidth[i]  = t.DrawLatexNDC( 0.5, 0.5, s.c_str() )->GetXsize() * canvas_size;
    charheight[i] = t.DrawLatexNDC( 0.5, 0.5, s.c_str() )->GetYsize() * canvas_size;
  }

  for( auto& p : charwidth ){
    p.second            /= mult * fontsize;
    charheight[p.first] /= mult * fontsize;
  }

  for( const auto& p : charwidth ){
    printf( "{'%c', %.8lf},", p.first, p.second  );
  }

  std::cout << "Done Generating!" << std::endl;

  std::cout << "Checking Kerning" << std::endl;
  const double compound = t.DrawLatexNDC( 0.5, 0.5, "AB" )->GetXsize() * canvas_size;
  std::cout << compound << " "
            << compound - (charwidth['A'] + charwidth['B'])* fontsize
            << std::endl;

  std::cout << "Checking Matches" << std::endl;
  std::string test1 = "This is a test";
  std::cout << usr::plt::EstimateLatexWidth( test1 ) * fontsize << " "
            << t.DrawLatex( 0.5, 0.5, test1.c_str() )->GetXsize() * canvas_size<< std::endl;

  // ------------- Height matching stuff
  std::cout << "\n\n\n" << "Height stuff" << std::endl ;
  for( const auto& p : charheight ){
    printf( "{'%c', %.8lf},", p.first, p.second  );
  }
  std::cout << std::endl;

  std::string test2 = "#frac{()}{()}";
  std::cout << usr::plt::EstimateLatexHeight( test2 ) * fontsize << " "
            << t.DrawLatex( 0.3, 0.3, test2.c_str() )->GetYsize() * canvas_size << std::endl;
  std::string test3 = "#int";
  std::cout << usr::plt::EstimateLatexHeight( test3 ) * fontsize << " "
            << t.DrawLatex( 0.3, 0.3, test3.c_str() )->GetYsize() * canvas_size << std::endl;
  std::string test4 = "a^{a}";
  std::cout << usr::plt::EstimateLatexHeight( test4 ) * fontsize << " "
            << t.DrawLatex( 0.3, 0.3, test4.c_str() )->GetYsize() * canvas_size << std::endl;
  std::string test5 = "A^{a}";
  std::cout << usr::plt::EstimateLatexHeight( test5 ) * fontsize << " "
            << t.DrawLatex( 0.3, 0.3, test5.c_str() )->GetYsize() * canvas_size << std::endl;
  std::string test6 = "()^{a}";
  std::cout << usr::plt::EstimateLatexHeight( test6 ) * fontsize << " "
            << t.DrawLatex( 0.3, 0.3, test6.c_str() )->GetYsize() * canvas_size << std::endl;

  c.SaveAs("test.pdf");

  return 0;
}
