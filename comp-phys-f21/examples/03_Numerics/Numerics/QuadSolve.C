//******************************************************************************
// QuadSolve: compares solutions to quadratic equations
//            reads output of tests
//
// Inputs
//  fName  file name to read comparison values from
//
// Modifications
// 04-Sep-20  make plots more readable
// 27-Jun-18  update to ROOT v6
// 01-Dec-09  created
//******************************************************************************
R__ADD_INCLUDE_PATH($CPCODE/library/macro)   // necessary to get CompPhys ROOT macros
#include "cpROOT.C"

#include <string>
#include <iostream>
#include <cmath>

const int MAXPTS = 50;

void PadSetup( double margin );

//==================================================================
void QuadSolve( string file="tests.out" )
{
  // setup standard functions and default plotting options
  cpNicePlot();

  //---------------------
  // Open file
  //---------------------
  ifstream ifl( file.data() );
  if( !ifl ){
    cout << "Could not find file: " << file << endl;
    return;
  } // !ifl

  //----------------------
  // Read 1st line: 
  // Format: "a: A  b: B"
  //----------------------
  string tmp,line;
  double a,b,c;
  ifl >> tmp >> a >> tmp >> b;
  getline( ifl,line );            // cleans up whitespace at end of line

  cout << a << " " << b << " |" << line << "|" << endl;

  //----------------------
  // Now read the quad eqn solution comparison results
  // Format: "c: C  x1: Dev11 Dev12  x2: Dev21 Dev22"
  //----------------------
  double x[MAXPTS];
  double y11[MAXPTS],y12[MAXPTS],y21[MAXPTS],y22[MAXPTS];
  double dev11,dev12,dev21,dev22;
  int n = 0;
  while( ifl >> tmp >> c >> tmp >> dev11 >> dev12 >> tmp >> dev21 >> dev22 ){
    x[n]   = c ;
    y11[n] = (dev11 != 0) ? abs(dev11) : 1.0e-20;  // protect dev <= 0 in log plot
    y12[n] = (dev12 != 0) ? abs(dev12) : 1.0e-20;
    y21[n] = (dev21 != 0) ? abs(dev21) : 1.0e-20;
    y22[n] = (dev22 != 0) ? abs(dev22) : 1.0e-20;

    cout << n << " " << c << " " << dev11 << " " << dev12
	 << " " << dev21 << " " << dev22 << endl;
    cout << n << " " << x[n] << " " << y11[n] << " " << y12[n] 
	 << " " << y21[n] << " " << y22[n] << endl;

    n++;
    if( n >= MAXPTS ) break;  // check for too many points
  } // while

  //-------------------------
  // Create graphs of the results
  //-------------------------
  TGraph *g11 = new TGraph( n,x,y11 );
  cpGraphStyle( g11,kRed,2,20 );         // graph setup: marker-style(20), color(red)
  TGraph *g12 = new TGraph( n,x,y12 );
  cpGraphStyle( g12,kBlue,2,20 );

  TGraph *g21 = new TGraph( n,x,y21 );
  cpGraphStyle( g21,kRed,2,20 );
  TGraph *g22 = new TGraph( n,x,y22 );
  cpGraphStyle( g22,kBlue,2,20 );


  //-------------------------
  // Create Dummy Histograms: useful for plotting graphs on standard axes
  //-------------------------
  double asiz = 0.07;    // size of axis labels
  string title1( "x_{+,-} = [ -b #pm #sqrt{b^{2}-4ac} ] / 2a" );
  TH2F *hDum1 = new TH2F( "hDum1",title1.data(), 10,x[n-1]/10,x[0]*10, 10,1.0e-21,10.0 );
  cpHistStyle( hDum1,1,2,0,"",asiz,"c","deviation from exact" );

  string title2( "x_{+,-} = -2c / [ b #pm #sqrt{b^{2}-4ac} ]" );
  TH2F *hDum2 = new TH2F( "hDum2",title2.data(), 10,x[n-1]/10,x[0]*10, 10,1.0e-21,10.0 );
  cpHistStyle( hDum2,1,2,0,"",asiz,"c","deviation from exact" );

  //-------------------------
  // Draw Graphs
  //-------------------------
  double marg = 0.15;
  // create canvas w/ 2 pads in y
  TCanvas* canv = new TCanvas( "canv","Quadratic" );
  canv->Divide( 1,2 );

  // 1st pad
  canv->cd(1);
  PadSetup( marg );   // set up the pad
  hDum1->Draw();      // draw histo to set axes
  g11->Draw( "LP" );  // draw graphs as points(P) connected by lines(L)
  g12->Draw( "LP" );

  TLegend* leg = new TLegend( 0.2,0.3,0.3,0.5 );  // create legend at x1,y1,x2,y2
  leg->AddEntry( g11,"x_{+}","LP" );              // add graph to legend w/ explanation
  leg->AddEntry( g12,"x_{-}","LP" );
  leg->Draw();

  // 2nd pad
  canv->cd(2);
  PadSetup( marg );   // set up the pad
  hDum2->Draw();      // draw histo to set axes
  g21->Draw( "LP" );  // draw graphs as points(P) connected by lines(L)
  g22->Draw( "LP" );

  return;
} // QuadSolve

//==================================================================
// PadSetup: adjusts pad parameters to make plots more readable
//==================================================================
void PadSetup( double margin )
{
  // increase margins to allow bigger axis labels
  gPad->SetLeftMargin( margin );
  gPad->SetBottomMargin( margin );
  gPad->SetTopMargin( margin );

  // set log-log display for this pad
  gPad->SetLogx();
  gPad->SetLogy();
  
  // draw a grid on the graphs
  gPad->SetGrid();
} // PadSetup
