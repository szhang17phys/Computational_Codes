//*************************************************************
// cpROOT: standard ROOT setup for Computational Physics
//
// Modifications
//
// 22-Sep-20  add cpPadSetup
// 03-Sep-20  cpAxisLimits: add nDec option, remove trap for min or max = 0
// 26-Sep-18  add cpCleanUp
// 03-Apr-10  add cpFuncStyle
// 02-Feb-10  add cpAxisLimits
// 08-Jan-10  add cpRoundUp, cpRoundDown
// 25-Sep-09  created
//*************************************************************
#include <cmath>
#include <string>
#include <utility>

#include "TDirectory.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1.h"
#include "TROOT.h"
#include "TStyle.h"

//------------------------------------------------------------
// NicePlot: avoids grey background on figures
//------------------------------------------------------------
void cpNicePlot()
{
  const Int_t bgrd_col = 10;
  gStyle->SetCanvasColor( bgrd_col );      // background color to white
  gStyle->SetPadColor( bgrd_col );
  gStyle->SetFrameFillColor( bgrd_col );

  gStyle->SetOptStat( 0 );       // no histo stat's by default
  gStyle->SetPadTickX( 1 );      // tick marks on right/upper axes
  gStyle->SetPadTickY( 1 );
} // cpNicePlot

//------------------------------------------------------------
// cpCleanUp: deletes canvases and histo's left over from previously running a macro
//------------------------------------------------------------
void cpCleanUp()
{
  // Histograms
  TIter nextH( gDirectory->GetList() );
  TObject* obj;
  while( (obj = (TObject*)nextH()) )
    delete obj;

  // Canvases
  TIter nextC(gROOT->GetListOfCanvases());
  TCanvas* c;
  while( (c = (TCanvas*)nextC()) )
    delete c;
} // cpCleanUp

//--------------------------------------------------------------
// cpGraphStyle, cpHistStyle, cpFuncStyle: 
// standard setup for TGraphs and THists
//   c      = color
//   w      = line width (2 is a good value)
//   m      = marker style (20 is a good value if you want markers)
//            unused in cpFuncStyle
//   title  = title for plot
//   asiz   = size of axis and title labels
//   xtitle = title for x-axis
//   ytitle = title for y-axis
//--------------------------------------------------------------
void cpGraphStyle( TGraph *g, int c=1, double w=0, int m=0, std::string title="",
		   double asiz=0, std::string xtitle="", std::string ytitle="" )
{
  g->SetMarkerColor( c );    // colors
  g->SetLineColor( c );

  if( w > 0 )                // line width
    g->SetLineWidth( w );

  if( m > 0 )                // marker style
    g->SetMarkerStyle( m );

  if( title != "" )          // Graph Title
    g->SetTitle( title.data() );

  if( asiz > 0 ){            // title & axis label size
    gStyle->SetTitleFontSize( asiz );
    g->GetXaxis()->SetLabelSize( asiz );
    g->GetXaxis()->SetTitleSize( asiz );
    g->GetYaxis()->SetLabelSize( asiz );
    g->GetYaxis()->SetTitleSize( asiz );
  } // asiz > 0

  if( xtitle != "" )         // X-axis title
    g->GetXaxis()->SetTitle( xtitle.data() );

  if( ytitle != "" )         // X-axis title
    g->GetYaxis()->SetTitle( ytitle.data() );
}

void cpHistStyle( TH1 *g, int c=1, double w=0, int m=0, std::string title="",
		  double asiz=0, std::string xtitle="", std::string ytitle="" )
{
  g->SetMarkerColor( c );    // colors
  g->SetLineColor( c );

  if( w > 0 )                // line width
    g->SetLineWidth( w );

  if( m > 0 )                // marker style
    g->SetMarkerStyle( m );

  if( title != "" )          // Graph Title
    g->SetTitle( title.data() );

  if( asiz > 0 ){            // title & axis label size
    gStyle->SetTitleFontSize( asiz );
    g->GetXaxis()->SetLabelSize( asiz );
    g->GetXaxis()->SetTitleSize( asiz );
    g->GetYaxis()->SetLabelSize( asiz );
    g->GetYaxis()->SetTitleSize( asiz );
  } // asiz > 0

  if( xtitle != "" )         // X-axis title
    g->GetXaxis()->SetTitle( xtitle.data() );

  if( ytitle != "" )         // X-axis title
    g->GetYaxis()->SetTitle( ytitle.data() );
}

void cpFuncStyle( TF1 *g, int c=1, double w=0, int m=0, std::string title="",
		  double asiz=0, std::string xtitle="", std::string ytitle="" )
{
  g->SetLineColor( c );      // line color

  if( w > 0 )                // line width
    g->SetLineWidth( w );

  if( title != "" )          // Graph Title
    g->SetTitle( title.data() );

  if( asiz > 0 ){            // title & axis label size
    gStyle->SetTitleFontSize( asiz );
    g->GetXaxis()->SetLabelSize( asiz );
    g->GetXaxis()->SetTitleSize( asiz );
    g->GetYaxis()->SetLabelSize( asiz );
    g->GetYaxis()->SetTitleSize( asiz );
  } // asiz > 0

  if( xtitle != "" )         // X-axis title
    g->GetXaxis()->SetTitle( xtitle.data() );

  if( ytitle != "" )         // X-axis title
    g->GetYaxis()->SetTitle( ytitle.data() );
}

//-------------------------------------------------------------
// cpSciNot: returns a,b with val = a * 10^b
//-------------------------------------------------------------
double cpSciNot( double val, int& b )
{
  b = 0;
  double absval = abs( val );
  double sgn    = val / absval;

  // abs(val) > 10: divide by ten until result is < 10
  if( absval > 10 ){
    while( absval > 10 ){
      absval = absval/10;
      b++;
    } // while
  } // absval > 10

  // abs(val) < 1: multiply by ten until result > 1
  else if( absval < 1 ){
    while( absval < 1 ){
      absval = absval*10;
      b--;
    } // while
  } // absval < 1

  return( sgn* absval );
} // cpSciNot

//-------------------------------------------------------------
// cpRound: round a double (value) to nearest nDec dec places
//   val  = value to round
//   dir  = direction 1=up; -1=down
//   nDec = number of decimal places to retain
//-------------------------------------------------------------
double cpRound( double val, int dir=1, int nDec=0 )
{
  if( val == 0.0 ) return val;

  int b;
  double mult = pow( 10.0,nDec );   // 10^nDec
  double a    = cpSciNot( val,b );  // val = a * 10^b

  // shift a nDec places, round to nearest int
  // then shift back
  double adja;
  if( dir > 0 ) adja = ceil( a * mult ) / mult;  // round up
  else          adja = floor( a * mult ) / mult; // round down

  return( adja * pow( 10.0,b ) );
} // cpRound

//=======================================================================
// AxisLimits: sets up limits for plot axis based on min and max values
//=======================================================================
std::pair<double,double> cpAxisLimits( double min, double max, int nDec=0 )
{
  std::pair<double,double> lim;

  // default axes (min==max) are -1 -- +1
  if( min >= max ){ 
    lim.first  = -1;
    lim.second =  1;
    return lim;
  } // min == max

  // wrong order - change them back
  if( min > max ){
    double tmp = min;
    min = max;
    max = tmp;
  } // min > max

  // otherwise round min/max down/up to nearest decimal place
  min = cpRound( min,-1,nDec );
  max = cpRound( max, 1,nDec );

  // deal with cases where a bound is zero (doesn't plot well)
  //if( min == 0 ) min = -max / 10;
  //if( max == 0 ) max = -min / 10;

  lim.first  = min;
  lim.second = max;
  return lim;
} // AxisLimits

//==================================================================
// cpPadSetup: adjusts pad margin parameters to allow bigger axis labels
//==================================================================
void cpPadSetup( double margin )                      // adjust left, bottom, right by same amount
{
  // increase margins to allow bigger axis labels
  gPad->SetLeftMargin( margin );
  gPad->SetBottomMargin( margin );
  gPad->SetTopMargin( margin );
} // cpPadSetup

void cpPadSetup( double lmarg, double rmarg, double bmarg, double tmarg )
{
  // increase margins to allow bigger axis labels
  if( lmarg > 0 ) gPad->SetLeftMargin( lmarg );
  if( rmarg > 0 ) gPad->SetRightMargin( rmarg );
  if( bmarg > 0 ) gPad->SetBottomMargin( bmarg );
  if( tmarg > 0 ) gPad->SetTopMargin( tmarg );
} // cpPadSetup
