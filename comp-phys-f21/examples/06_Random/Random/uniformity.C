//*******************************************************************
// uniformity: plots results of uniformity test
//
// Input
//   string fName: file name for input
//
// Modifications
// 22-Sep-20  increase axis label size
// 25-Sep-18  updated for Fall 2018
//*******************************************************************
R__ADD_INCLUDE_PATH($CPCODE/library/macro)   // necessary to get CompPhys ROOT macros
#include "cpROOT.C"

#include <iostream>
#include <fstream>
#include <string>

void uniformity( string fName )
{
  cpNicePlot();
  cpCleanUp();

  //open file
  ifstream ifl( fName.data() );

  // read header
  string title,line,tmp;
  int nBins,nTrials;
  getline( ifl,title );
  ifl >> tmp >> nBins >> tmp >> nTrials;

  // expected average
  double ave = (double)nTrials / (double)nBins;

  // read histogram data
  const int MAX_BINS = 200;
  double x[MAX_BINS] = { 0 };
  double y[MAX_BINS] = { 0 };
  for( int i = 0; i < nBins; ++i )
    ifl >> x[i] >> y[i];

  // book and fill histogram
  gStyle->SetOptStat();
  double dx  = x[1] - x[0];
  double xlo = x[0] - dx/2;        // histo limits
  double xhi = x[nBins-1] + dx/2;
  double asiz = 0.05;              // axis label size
  double marg = 0.15;              // margins
  TH1F* h = new TH1F( "h",title.data(),nBins,xlo,xhi );
  cpHistStyle( h, kBlack, 2, 20, title.data(), asiz, "PRN", "Frequency" );
  h->SetFillStyle( 1001 );
  h->SetFillColor( kYellow );
  h->SetMinimum( 0 );
  for( int i = 0; i < nBins; ++i ) h->Fill( x[i],y[i] );

  // Plot
  TCanvas* canv = new TCanvas( "canv","uniformity" );
  canv->SetGrid();
  cpPadSetup( marg );
  h->Draw( "hist" );

  TLine* lAve = new TLine( xlo,ave,xhi,ave );
  lAve->SetLineWidth( 4 );
  lAve->SetLineStyle( 1 );
  lAve->SetLineColor( kRed );
  lAve->Draw();

} // uniformity
