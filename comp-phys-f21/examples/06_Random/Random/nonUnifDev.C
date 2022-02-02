//*******************************************************************
// nonUnifDev: makes plots of Non Uniform Deviates
//
// Input
//   nFiles  number of files to plot
//
// Modifications
// 22-Sep-20  increase axis label size
// 19-Oct-18  change histogram options for doPoiss
// 09-Oct-18  updated for Fall 2018
// 22-Apr-10  add fitting option
//*******************************************************************
R__ADD_INCLUDE_PATH($CPCODE/library/macro)   // necessary to get CompPhys ROOT macros
#include "cpROOT.C"

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>

void nonUnifDev( int nFiles, int doPoiss=0, int doFit=0 )
{
  cpNicePlot();
  cpCleanUp();

  double asiz = 0.05;  // axis label size
  double marg = 0.15;  // margin size

  //------------------------
  // Loop over files
  //------------------------
  const int MAXFILES = 5;
  string fName;
  string sPRNG[MAXFILES];
  string sModel;
  string line,tmp;

  const int MAXBINS = 200;
  int nPts,nBins;
  double xlo,xhi;
  double xgen[MAXBINS];
  double ygen[MAXBINS];
  double egen[MAXBINS];
  double yprd[MAXBINS];
  double dum[MAXBINS] = { 0 };  // 0 error bars in x

  TGraph* gPred;
  TH1F* hPred;     // for Poisson Distrib
  TGraphErrors* gGen[MAXFILES];

  for( int ifile = 0; ifile < nFiles; ++ifile ){
    // get filename and sPRNG
    cout << "File " << ifile << " name: ";
    cin >> fName;
    cout << "      PRNG: ";
    cin >> sPRNG[ifile];

    //open file
    ifstream ifl( fName.data() );

    //-----------------
    // read header
    //-----------------
    ifl >> sModel >> tmp >> nPts >> tmp >> nBins >> tmp >> xlo >> tmp >> xhi; 

    //------------------
    // read bin contents
    //------------------
    for( int i = 0; i < nBins; ++i ){
      ifl >> xgen[i] >> ygen[i] >> egen[i] >> yprd[i];
    }

    // generated data graph
    gGen[ifile] = new TGraphErrors( nBins, xgen,ygen,dum,egen );
    cpGraphStyle( gGen[ifile],ifile+2,2,20 );

    // only create prediction graph from first file
    if( ifile == 0 ){
      gPred = new TGraph( nBins,xgen,yprd );
      cpGraphStyle( gPred,1,2 );

      if( doPoiss ){
	hPred = new TH1F( "hPred","",nBins,xlo,xhi );
	cpHistStyle( hPred,1,2 );
	for( int i = 0; i < nBins; ++i )
	  hPred->Fill( xgen[i],yprd[i] );
      } // doPoiss
    } // ifile = 0

    ifl.close();
  }  // ifile loop

  //-----------------------
  // Draw Graphs
  //-----------------------
  TCanvas* canv = new TCanvas( "canv","NonUniform" );
  cpPadSetup( marg );
  canv->SetGrid();

  // limits for histograms
  double ylo,yhi;
  if( doPoiss ){           // Poisson Distribution
    ylo = 0;
    yhi = 0.3;
  } // doPoiss

  else{
    canv->SetLogy();
    ylo = 1.0e-5;
    yhi = 0.1;
  } // !doPoiss

  if( doFit ) gStyle->SetOptFit();

  TH2F* hDum = new TH2F( "hDum",sModel.data(), 10,xlo,xhi, 10,ylo,yhi );
  cpHistStyle( hDum,1,2,20,sModel,asiz,"x","N_{bin}/N" );
  hDum->Draw();

  for( int ifile = 0; ifile < nFiles; ++ifile ){
    if( doFit ){ 
      gGen[ifile]->Fit( "gaus" );
      gGen[ifile]->GetFunction( "gaus" )->SetLineColor( kBlue );
    } // doFit
    gGen[ifile]->Draw("P");
  } // ifile
  if( doPoiss ) hPred->Draw( "same hist" );
  else          gPred->Draw( "C" );

  TLegend* leg = new TLegend( 0.75,0.8,0.95,0.95 );
  leg->AddEntry( gPred,"PDF","L" );
  for( int ifile = 0; ifile < nFiles; ++ifile )
    leg->AddEntry( gGen[ifile],sPRNG[ifile].data(),"P" );
  if( doFit == 0 ) leg->Draw();

} // filling
