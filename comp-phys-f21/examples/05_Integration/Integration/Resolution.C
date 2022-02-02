//*********************************************************************
// Resolution: makes plots of Exponential conv w/ a Gaussian
//
// Modifications
// 26-Sep-18  add cpCleanUp
// 13-Aug-18  updated for Fall 2018
// 19-Feb-10  change passed parameters
// 02-Feb-10  created
//*********************************************************************
R__ADD_INCLUDE_PATH($CPCODE/library/macro)   // necessary to get CompPhys ROOT macros
#include "cpROOT.C"

#include <iostream>
#include <cmath>
#include <string>
#include <vector>

//=====================================================================
void Resolution( int nFiles=1, double lambda=1.0 )
{
  cpNicePlot();
  cpCleanUp();

  const int MAXFILES = 10;     // File information
  string file[MAXFILES];       //   filenames
  string sLim[MAXFILES];       //   string with limits of int

  const int MAXBINS = 250;  // bin information
  int nBins,nSteps;
  double tlo,thi;
  double x[MAXBINS];
  double y[MAXBINS];

  TGraph* gRes[MAXFILES];     // the graphs

  //--------------------
  // Loop over Files
  //--------------------
  string line,tmp;
  for( int fnum = 0; fnum < nFiles; ++fnum ){
    // query for filename, and label
    //   note infinity = #infty
    cout << "File:" << fnum << " filename label? ";
    cin >> file[fnum] >> sLim[fnum];
    ifstream ifl( file[fnum].data() );     // Open File
    if( !ifl ){
      cout << "Could not find file: " << file[fnum] << endl;
      return;
    } // !ifl

    // Read File header
    getline( ifl,line );                              // line 1
    ifl >> tmp >> tlo >> thi >> tmp >> nBins >> tmp;  // line 2

    // loop over bins
    for( int i = 0; i < nBins; ++i )
      ifl >> x[i] >> y[i] >> nSteps;

    // create graph for this file
    gRes[fnum] = new TGraph( nBins,x,y );
    cpGraphStyle( gRes[fnum],2+fnum,2,20 );

    ifl.close();
  } // fnum loop

  //----------------------
  // Exponential Function & Dummy histogram for plot
  //----------------------
  TF1* fExp = new TF1( "fExp","exp(-[0]*x)/[0]",0,thi );
  fExp->SetParameter(0,lambda);
  fExp->SetLineWidth( 2 );

  TH2F* hDum = new TH2F( "hDum","Exp conv w/ Gaussian",
  			 nBins,tlo,thi, 10,1.0e-5,1.2/lambda );
  cpHistStyle( hDum,1,1,0,"",0.10,"#lambda t","Convolution" );

  //-------------------
  // Draw the Graph(s)
  //-------------------
  TCanvas* canv = new TCanvas( "canv","Resolution" );
  canv->Divide( 1,2 );
  TLegend* leg = new TLegend( 0.70, 0.90-0.15*nFiles/2, 0.90, 0.90 );

  for( int i = 0; i < 2; ++i ){   // linear and log scales
    canv->cd( i+1 );
    gPad->SetLeftMargin( 0.2 );
    gPad->SetBottomMargin( 0.2 );
    gPad->SetGrid();
    if( i == 1 ) gPad->SetLogy();

    hDum->Draw();
    for( int fnum = 0; fnum < nFiles; ++fnum ){
      gRes[fnum]->Draw( "P" );
      if( i == 0 ) leg->AddEntry( gRes[fnum],sLim[fnum].data(),"P" );
    } // fnum loop
    fExp->Draw( "same" );

    if( i == 0 ) leg->Draw();
  } // i loop

} // Resolution
