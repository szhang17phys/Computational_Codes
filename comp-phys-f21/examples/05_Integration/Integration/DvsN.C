//*********************************************************************
// DvsN.C: Deviation vs N(steps)
//
// Modifications
// 28-Sep-20  increase axis label size
// 03-Aug-20  increase marker size for Romberg graph
// 13-Aug-18  updated for Fall 2018
// 26-Mar-10  plot log_10(N) instead of log_2(N)
// 02-Feb-10  created
//*********************************************************************
R__ADD_INCLUDE_PATH($CPCODE/library/macro)   // necessary to get CompPhys ROOT macros
#include "cpROOT.C"

#include <iostream>
#include <cmath>
#include <string>
#include <vector>

// Function prototypes
void updateLim( double x1, double x2, double x3, double& xmin, double& xmax );

//=====================================================================
void DvsN( string file, string title="" )
{
  cpNicePlot();
  cpCleanUp();

  double asiz = 0.05;  // axis label size
  double marg = 0.15;  // margin size

  //---------------------
  // Read Filename / Open File
  //---------------------
  ifstream ifl( file.data() );
  if( !ifl ){
    cout << "Could not find file: " << file << endl;
    return;
  } // !ifl

  //----------------------
  // Read file header
  // Title of plot: read from file if not passes
  //----------------------
  int log2nStep;
  string line,tmp;
  ifl >> tmp >> tmp;                // line 1
  getline( ifl,line );
  if( title == "" ) title = line;

  ifl >> log2nStep;                 // line 2
  getline( ifl,line );

  //cout << title << " - " << log2nStep << " intervals" << endl;

  //-----------------
  // Read Data
  //------------------
  const int MAXINTVL = 100;
  int N;
  double x[MAXINTVL];
  double qR[MAXINTVL], dR[MAXINTVL];
  double qT[MAXINTVL], dT[MAXINTVL];
  double qS[MAXINTVL], dS[MAXINTVL];
  double xRo,qRo,dRo;
  double ybase = -20;    // set log10(0) to this value
  double ymin  = 100;
  double ymax  = -100;

  // Rectangle, Trapezoid, Simpson
  for( int i = 0; i < log2nStep; ++i ){
    ifl >> x[i] >> N >> qR[i] >> dR[i]>> tmp >> qT[i] >> dT[i] 
	>> tmp >> qS[i] >> dS[i];
    //x[i] = log10( (double)N );
    //cout << x[i] << " " << N << " " << qR[i] << " " << dR[i] << " " 
    //	 << qT[i] << dT[i] << " " << qS[i] << " " << dS[i] << endl;
    dR[i] = ( dR[i] != 0 ) ? log10( abs(dR[i]) ) : ybase;
    dT[i] = ( dT[i] != 0 ) ? log10( abs(dT[i]) ) : ybase;
    dS[i] = ( dS[i] != 0 ) ? log10( abs(dS[i]) ) : ybase;
    updateLim( dR[i],dT[i],dS[i],ymin,ymax );
  } // i loop

  // Romberg
  ifl >> xRo >> N >> qRo >> dRo;
  dRo = ( dRo != 0 ) ? log10( abs(dRo) ) : ybase;
  if( dRo < ymin ) ymin = dRo;
  if( dRo > ymax ) ymax = dRo;

  ifl.close();


  //------------------
  // Create the Graphs for this file
  //------------------
  TGraph* gR = new TGraph( log2nStep,x,dR );
  cpGraphStyle( gR,2,2,20 );
  TGraph* gT = new TGraph( log2nStep,x,dT );
  cpGraphStyle( gT,3,2,20 );
  TGraph* gS = new TGraph( log2nStep-1,&(x[1]),&(dS[1]) ); // N=1 Simpson meaningless
  cpGraphStyle( gS,4,2,20 );
  TGraph* gRo = new TGraph( 1,&xRo,&dRo );
  cpGraphStyle( gRo,1,2,20 );
  gRo->SetMarkerSize( 2 );

  //----------------------
  // Create dummy histogram that spans min and max x and y values
  //----------------------

  // set axis limits to min and max values rounded to nearest 1 decimal place
  pair<double,double> ylim = cpAxisLimits( ymin,ymax );

  TH2F* hDum = new TH2F( "hDum",title.data(),
			 log2nStep,x[0],x[log2nStep-1],
			  10,ylim.first,ylim.second );
  cpHistStyle( hDum,1,1,0,"",asiz,"log_{2}(N)","log_{10}(Num/Ana-1)" );

  //-------------------
  // Draw the Graph(s)
  //-------------------
  TCanvas* canv = new TCanvas( "canv","1DQuadrature" );
  cpPadSetup( marg );
  canv->SetGrid();

  hDum->Draw();
  gR->Draw( "PL" );  gT->Draw( "PL" );  gS->Draw( "PL" );
  gRo->Draw( "P" );

  TLegend* leg = new TLegend( 0.75,0.80,0.99,0.99 );
  leg->AddEntry( gR,"Rectangle","PL" );
  leg->AddEntry( gT,"Trapezoid","PL" );
  leg->AddEntry( gS,"Simpson","PL" );
  leg->AddEntry( gRo,"Romberg","P" );
  leg->Draw();

} // DvsN



//=======================================================================
// updateLim: updates min and max values base on inputs
//=======================================================================
void updateLim( double x1, double x2, double x3, double& xmin, double& xmax )
{
  if( x1 < xmin ) xmin = x1;  // minimum
  if( x2 < xmin ) xmin = x2;
  if( x3 < xmin ) xmin = x3;

  if( x1 > xmax ) xmax = x1;  // maximum
  if( x2 > xmax ) xmax = x2;
  if( x3 > xmax ) xmax = x3;
} // updateLim
