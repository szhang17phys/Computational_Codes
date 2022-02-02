//*******************************************************************
// filling.C: plots empty fraction vs t for 2D filling test
//
// Input
//   nFiles  number of files to plot
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
#include <cmath>

void filling( int nFiles )
{
  cpNicePlot();
  cpCleanUp();

  double asiz = 0.05;  // axis label size
  double marg = 0.15;  // margin size

  TCanvas* canv = new TCanvas( "canv","filling" );
  cpPadSetup( marg );
  canv->SetGrid();
  TLegend* leg = new TLegend( 0.7,0.6,0.95,0.95 );

  //------------------------
  // Loop over files
  //------------------------
  const int MAXFILES = 10;
  const int MAXSTEPS = 100;
  string fName;
  string title[MAXFILES];
  string line,tmp;
  int L,nSteps;
  double t[MAXSTEPS];
  double y[MAXFILES][MAXSTEPS];
  TGraph* gFrac[MAXFILES];

  for( int ifile = 0; ifile < nFiles; ++ifile ){
    // get filename and title
    cout << "File " << ifile << " name: ";
    cin >> fName;
    cout << "       title: ";
    cin >> title[ifile];

    //open file
    ifstream ifl( fName.data() );

    // read header
    getline( ifl,line );
    ifl >> tmp >> L >> tmp >> nSteps;
    //cout << nSteps << endl;

    // read entries
    for( int i = 0; i < nSteps; ++i ){
      ifl >> tmp >> t[i] >> tmp >> tmp >> tmp >> y[ifile][i];
      //cout << t[i] << " " << y[ifile][i] << endl;
    }

    // create and draw graph
    gFrac[ifile] = new TGraph( nSteps,t,&(y[ifile][0]) );
    cpGraphStyle( gFrac[ifile],ifile+2,2,20,"Empty Fraction",asiz,"t","f(Empty)" );
    if( ifile == 0 ) gFrac[ifile]->Draw("APL");
    else             gFrac[ifile]->Draw("PL");
    leg->AddEntry( gFrac[ifile],title[ifile].data(),"PL" );

    ifl.close();
  }  // ifile loop

  // exponential decay curve
  double norm = 1;
  TF1* fExp = new TF1( "fExp","[0]*exp(-x)",t[0],t[nSteps-1] );
  fExp->SetParameter( 0,norm );
  fExp->SetLineWidth( 2 );
  fExp->SetLineColor( kBlack );
  fExp->Draw( "same" );

  leg->AddEntry( fExp,"e^{-t}","L" );
  leg->Draw();

} // filling
