//*********************************************************************
// YvsX.C: ROOT macro for plotting ODE test results
// Parameters
//   model   1 Decay
//           2 Falling
//           3 Simple Harmonic Oscillator
//           4 Falling w/ Air Resist
//   opt     1 (./tests 1)  X vs t
//           2 (./tests 2) (Nana-Nnum) vs h
//           3 (./tests 3) (Nana-Nnum)/Nana vs nStep
//           4 (./tests 2) (Nana-Nnum)/Nana vs h
//   nFiles  number of tests.cpp output files to use in plots
//
// Modifications
// 18-Sep-20  make sure analytic curve is on top of plot
// 03-Sep-20  use cpAxisLimits(), remove AxisLimits()
// 22-Sep-10  change to AxisLimits()
//            do not use AxisLimits for x-axes
// 22-Jan-10  add model=4 (air resist)
// 20-Jan-10  add opt=4
// 18-Jan-10  add AxisLimits()
//*********************************************************************
R__ADD_INCLUDE_PATH($CPCODE/library/macro)   // necessary to get CompPhys ROOT macros
#include "cpROOT.C"

#include <iostream>
#include <cmath>
#include <string>
#include <vector>

//-------------------------------------------------------------------
// Function templates: needed for all functions that do not correspond to the file name
//--------------------------------------------------------------------
bool GetDataLine( ifstream&, double &, double&, double&, double*, double*, double* );
pair<double,double> AxisLimits( double, double );
string PlotTitle( int, int, int );


//--------------------------------------------------------------------
// Constants
//--------------------------------------------------------------------
const int idDecay   = 1;     // radioactive decay
const int idFall    = 2;     // 1D falling w/out air resist
const int idSHO     = 3;     // simple harmonic oscillator
const int idFallAir = 4;     // 1D falling w/ air resist

//=====================================================================
void YvsX( int model, int opt, int nFiles = 1 )
//void YvsX(1, 1, 1)
{
  // standard CompPhys ROOT setup
  cpNicePlot();
  cpCleanUp();

  //-----------------------------
  // Plotting options
  //-----------------------------
  bool doLegend = ( nFiles > 1 );             // draw legend
  bool doVE     = ( opt != 3 &&               // plot v and E too
		    ( model == idFall || model == idSHO ) );
  bool doV      = ( opt == 1 && model == idFallAir ); // plot x,v(t)

  const int MAXVALS = 500;    // max allowed number of h-values
  const int MAXFILE = 10;     // max number of files

  string file[MAXFILE];
  TGraph* gYvX[MAXFILE];
  TGraph* gYanavX[MAXFILE];
  TGraph* gVvX[MAXFILE];
  TGraph* gVanavX[MAXFILE];
  TGraph* gEvX[MAXFILE];
  TGraph* gEanavX[MAXFILE];

  //------------------------------
  // Loop over input files
  //------------------------------
  string line;
  double n,t,h;
  double pos[3],vel[3],ene[3];
  double x[MAXVALS];
  double y[MAXVALS],yana[MAXVALS];
  double v[MAXVALS],vana[MAXVALS];
  double e[MAXVALS],eana[MAXVALS];
  double xmin=1.0e15,xmax=-999;
  double ymin=1.0e15,ymax=-999;
  double vmin=1.0e15,vmax=-999;
  double emin=1.0e15,emax=-999;
  for( int fnum = 0; fnum < nFiles; ++fnum ){

    //---------------------
    // Read Filename / Open File
    //---------------------
    cout << "File " << fnum << " name? ";
    cin >> file[fnum];
    ifstream ifl( file[fnum].data() );
    if( !ifl ){
      cout << "Could not find file: " << file[fnum] << endl;
      return;
    } // !ifl

    //-------------------------
    // Read information for each h-value
    //-------------------------
    int nVals = 0;
    while( GetDataLine( ifl,n,h,t,pos,vel,ene ) ){
      switch( opt ){
      case 1:         // State vs t (numeric and analytic)
	x[nVals] = t; 
	y[nVals] = pos[1];     yana[nVals] = pos[0];
	v[nVals] = vel[1];     vana[nVals] = vel[0];
	e[nVals] = ene[1];     eana[nVals] = ene[0];
	break;
      case 2:         // Ana-Num vs h
	x[nVals] = h; 
	y[nVals] = pos[0] - pos[1];
	v[nVals] = vel[0] - vel[1];
	e[nVals] = ene[0] - ene[1];
	break;
      case 3:         // log(devX) vs log(nStep)
	x[nVals] = log10( n );
	y[nVals] = ( abs(pos[2]) > 0.0 ) ? log10( abs(pos[2]) ) : -20;
	getline( ifl,line );  // read clock time line
	break;
      case 4:         // Ana-Num/Ana vs h
	x[nVals] = h; 
	y[nVals] = pos[2];
	v[nVals] = vel[2];
	e[nVals] = ene[2];
	break;
      } // switch opt

      // find maxima and minima over all graphs
      // allows to plot them all on the same axes
      //if( nVals > 0 ){  // first point is often an anomaly
	if( x[nVals] < xmin ) xmin = x[nVals];
	if( x[nVals] > xmax ) xmax = x[nVals];
	if( y[nVals] < ymin ) ymin = y[nVals];
	if( y[nVals] > ymax ) ymax = y[nVals];
	if( v[nVals] < vmin ) vmin = v[nVals];
	if( v[nVals] > vmax ) vmax = v[nVals];
	if( e[nVals] < emin ) emin = e[nVals];
	if( e[nVals] > emax ) emax = e[nVals];
	//} // nVals > 0

      nVals++;  // update number of values read
    } // while

    //------------------
    // Create the Graphs for this file
    //------------------
    gYvX[fnum] = new TGraph( nVals,x,y );
    cpGraphStyle( gYvX[fnum],2+fnum,2,20 );
    gYanavX[fnum] = new TGraph( nVals,x,yana );
    cpGraphStyle( gYanavX[fnum],1,2 );

    gVvX[fnum] = new TGraph( nVals,x,v );
    cpGraphStyle( gVvX[fnum],2+fnum,2,20 );
    gVanavX[fnum] = new TGraph( nVals,x,vana );
    cpGraphStyle( gVanavX[fnum],1,2 );

    gEvX[fnum] = new TGraph( nVals,x,e );
    cpGraphStyle( gEvX[fnum],2+fnum,2,20 );
    gEanavX[fnum] = new TGraph( nVals,x,eana );
    cpGraphStyle( gEanavX[fnum],1,2 );

    ifl.close();  // close file
  } // fnum loop

  //----------------------
  // Create dummy histogram that spans min and max x and y values
  //----------------------
  double asiz = 0;          // need to change default axis label sizes
  if( doVE || doV ) asiz = 0.08;   //   for 3 or 2 pads in canvas

  // set axis limits to min and max values rounded to nearest 0 or 1 decimal place
  pair<double,double> xlim = cpAxisLimits( xmin,xmax,1 );
  pair<double,double> ylim = cpAxisLimits( ymin,ymax );

  cout << "X: " << xmin << "," << xmax << "  Axis: " << xlim.first << "," << xlim.second << endl;
  cout << "Y: " << ymin << "," << ymax << "  Axis: " << ylim.first << "," << ylim.second << endl;

  TH2F* hDumY = new TH2F( "hDumY",PlotTitle(model,opt,0).data(),
			  10,xlim.first,xlim.second,
			  10,ylim.first,ylim.second );
  cpHistStyle( hDumY,1,1,0,"",asiz );

  pair<double,double> vlim = cpAxisLimits( vmin,vmax );
  TH2F* hDumV = new TH2F( "hDumV",PlotTitle(model,opt,1).data(),
			  10,xlim.first,xlim.second,
			  10,vlim.first,vlim.second );
  cpHistStyle( hDumV,1,1,0,"",asiz );

  pair<double,double> elim = cpAxisLimits( emin,emax );
  TH2F* hDumE = new TH2F( "hDumE",PlotTitle(model,opt,2).data(),
			  10,xlim.first,xlim.second,
			  10,elim.first,elim.second );
  cpHistStyle( hDumE,1,1,0,"",asiz );

  //-------------------
  // Draw the Graph(s)
  //-------------------
  TCanvas* canv = new TCanvas( "canv","ODE" );
  TLegend* leg = new TLegend( 0.75,0.70,0.99,0.99 );
  if( doVE ) canv->Divide( 1,3 );
  if( doV  ) canv->Divide( 1,2 );

  //----------
  // Position or N Plot
  //---------
  canv->cd(1);  gPad->SetGrid();
  hDumY->Draw();  // dummy histogram sets the scale
  for( int fnum = 0; fnum < nFiles; ++fnum ){
    gYvX[fnum]->Draw( "P" );        // numeric results
    if( opt == 1 && fnum == 0 ){    // analytic results (only once)
      leg->AddEntry( gYanavX[fnum],"Analytic","L" );    // legend
    }
    leg->AddEntry( gYvX[fnum],file[fnum].data(),"P" );  // legend
  } // fnum loop
  gYanavX[0]->Draw( "L" );         // make sure analytic curve is on top
  if( doLegend ) leg->Draw();      // draw legend

  if( doVE || doV ){
    //---------
    // Velocity Plot
    //---------
    canv->cd(2);  gPad->SetGrid();
    hDumV->Draw();
    for( int fnum = 0; fnum < nFiles; ++fnum )
      gVvX[fnum]->Draw( "P" );        // numeric results
    gVanavX[0]->Draw( "L" );          // make sure analytic curve is on top

    //---------
    // Energy Plot
    //---------
    if( doVE ){
      canv->cd(3);  gPad->SetGrid();
      hDumE->Draw();
      for( int fnum = 0; fnum < nFiles; ++fnum )
	gEvX[fnum]->Draw( "P" );        // numeric results
      gEanavX[0]->Draw( "L" );          // make sure analytic curve is on top
    } // doVE
  } // doVE || doV
} // YvsX


//=======================================================================
// GetDataLine: reads output of ./tests
//   fills variables
//   returns false if no more lines
//=======================================================================
bool GetDataLine( ifstream& ifl, double & n, double& h, double& t,
		  double* x, double* v, double* E )
{
  string tmp;
  int nEqn;

  // defaults
  n = -1; h = 0; t = 0;
  for( int i = 0; i < 3; ++i ){
    x[i] = 0;  v[i] = 0;  E[i] = 0;
  } // i loop

  // read nSteps, h, t, nEqn
  ifl >> n >> tmp >> h >> tmp >> t >> tmp >> nEqn;

  // nothing on this line = end of file
  if( n < 0 ) return false;

  // read state vector
  double ana,num,dev;
  for( int i = 0; i < nEqn; ++i ){
    ifl >> tmp >> ana >> num >> dev;  // the state

    if( i == 0 ){    // position
      x[0] = ana; x[1] = num;  x[2] = dev;
    } // i == 0

    if( i == 1 ){    // velocity (if present)
      v[0] = ana; v[1] = num;  v[2] = dev;
    }
  } // i = nEqn

  // read energy
  ifl >> tmp >> E[0] >> E[1] >> E[2];

  // clean up after last value read
  getline( ifl,tmp );

  return true;
} // GetDataLine


//=======================================================================
// PlotTitle: returns title for the plot base on model, opt and plot
//=======================================================================
string PlotTitle( int model, int opt, int var=0 )
{
  //----------------------------
  // Titles for each plotting option
  // Note: root doesn't deal well with arrays of strings
  //       so can't initialize them in declaration
  //-----------------------------
  const int MAXOPT  = 10;
  string titleOpt[MAXOPT];
  titleOpt[1]  = "State vs t  ";
  titleOpt[2]  = "Ana-Num vs h ";
  titleOpt[3]  = "log[(Ana-Num)/Ana] vs log[N_{steps}] ";
  titleOpt[4]  = "(Ana-Num)/Ana vs h ";

  //----------------------------
  // IDs and titles for various models
  //----------------------------
  const int MAXMODEL  = 10;
  string titleModel[MAXMODEL];
  titleModel[idDecay]   = "Decay: ";
  titleModel[idFall]    = "Falling: ";
  titleModel[idSHO]     = "SHO: ";
  titleModel[idFallAir] = "Falling w/ AirResist: ";

  //----------------------------
  // The plot variable type
  //----------------------------
  string titleVar[3];
  if( model == idDecay ) 
    titleVar[0] = "(N)";
  else
    titleVar[0] = "(x)";
  titleVar[1] = "(v)";
  titleVar[2] = "(E)";

  //----------------------------
  // Construct the title
  //----------------------------
  string title;
  title = titleModel[model];
  title.append( titleOpt[opt] );
  title.append( titleVar[var] );
  return title;
} // PlotTitle
