//*******************************************************************
// autocor.C: makes autocorrelation plots
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
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>

void autocor( int nFiles )
{
  cpNicePlot();
  cpCleanUp();

  //------------------------
  // Loop over files
  //------------------------
  const int MAXFILES = 10;
  const int MAXEXPS  = 500;
  const int MAXLAGS  = 50;
  string fName;
  string title[MAXFILES];
  string line,tmp;
  int nExp,nMeas,kLo,kHi,kPrint;

  double CkExp[MAXEXPS];
  double CkMin  = 999, CkMax  = -999;
  double CkMinE = 999, CkMaxE = -999;
  double Clo,Chi;
  TH1F* hCk;

  double asiz = 0.05;  // axis label size
  double marg = 0.15;  // margin size

  double k[MAXLAGS];
  double aveCk[MAXLAGS];
  TGraph* gAveCk[MAXFILES];

  for( int ifile = 0; ifile < nFiles; ++ifile ){
    // get filename and title
    cout << "File " << ifile << " name: ";
    cin >> fName;
    cout << "       title: ";
    cin >> title[ifile];

    //open file
    ifstream ifl( fName.data() );

    //-----------------
    // read header
    //-----------------
    getline( ifl,line );
    ifl >> tmp >> nExp >> tmp >> nMeas >> tmp >> kLo >> tmp >> kHi 
	>> tmp >> kPrint;

    //------------------
    // read exp-by-exp results if kPrint > 0
    //------------------
    if( kPrint > 0 ){
      for( int i = 0; i < nExp; ++i ){
	ifl >> tmp >> tmp >> tmp >> tmp >> CkExp[i];
	//cout << i << " " << CkExp[i] << endl;
	if( CkExp[i] < CkMinE ) CkMinE = CkExp[i];  // find max & min for histo
	if( CkExp[i] > CkMaxE ) CkMaxE = CkExp[i];
      } // i = exp's
      Clo = cpRound( CkMinE,-1,1 );
      Chi = cpRound( CkMaxE, 1,1 );

      //cout << CkMin << " " << CkMax << endl;

      // only book and fill histogram for the first file
      if( ifile == 0 ){
	ostringstream osCkExp;
	osCkExp << "C(" << kPrint << ")";
	//cout << Clo << " " << Chi << " " << tCkExp << endl;
	gStyle->SetOptStat();
	hCk = new TH1F( "hCk","",25,3*Clo,3*Chi );
	cpHistStyle( hCk,1,2,0,title[ifile],asiz,osCkExp.str(),"" );
	hCk->SetFillStyle( 1001 );
	hCk->SetFillColor( kYellow );
	for( int i = 0; i < nExp; ++i )
	  hCk->Fill( CkExp[i] );

	// draw histogram
	TCanvas* canvE = new TCanvas( "canvE","Exps" );
	cpPadSetup( 0,0,marg,0 );
	canvE->SetGrid();
	hCk->Draw();
      } // ifile = 0
    } // kPrint > 0

    //----------------------
    // read C(k) averages
    //----------------------
    int nLag = 0;

    while( ifl >> tmp >> k[nLag] >> tmp >> aveCk[nLag] ){
      getline( ifl,line ); // rest of line
      if( aveCk[nLag] < CkMin ) CkMin = aveCk[nLag];  // find max & min
      if( aveCk[nLag] > CkMax ) CkMax = aveCk[nLag];
      //cout << nLag << " " << k[nLag] << " " << aveCk[nLag] << endl;
      nLag++;
    } // while

    //--------------------------
    // create graph
    //--------------------------
    gAveCk[ifile] = new TGraph( nLag,k,aveCk );
    cpGraphStyle( gAveCk[ifile],ifile+2,2,20,"<C(k)> vs k",asiz,"k","<C(k)>" );

    ifl.close();
  }  // ifile loop

  Clo = cpRound( CkMin,-1,1 );
  Chi = cpRound( CkMax, 1,1 );


  //-----------------------
  // Draw Graphs
  //-----------------------
  TCanvas* canv = new TCanvas( "canv","Auto-Corr" );
  cpPadSetup( marg );
  canv->SetGrid();

  gStyle->SetOptStat( 0 );
  TH2F* hDum = new TH2F( "hDum","<C(k)> vs k", 10,kLo-1,kHi+1, 10,Clo,Chi );
  cpHistStyle( hDum,1,2,20,"<C(k)> vs k",asiz,"k","<C(k)>" );
  hDum->Draw();

  for( int ifile = 0; ifile < nFiles; ++ifile )
    gAveCk[ifile]->Draw("PL");

  TLegend* leg = new TLegend( 0.7,0.6,0.95,0.95 );
  for( int ifile = 0; ifile < nFiles; ++ifile )
    leg->AddEntry( gAveCk[ifile],title[ifile].data(),"PL" );
  leg->Draw();

} // filling
