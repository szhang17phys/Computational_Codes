// Group Problem 2
// 9/15/2021

#include <iostream>
#include <iomanip>
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "cpROOT.hpp"

using namespace std;

int main ()
{
// Dummy data
//Int_t n = 3;
//Double_t x[n] = {1, 2, 3};
//Double_t y1[n] = {4, 5, 6};
//Double_t y2[n] = {7, 8, 9};

  Double_t x[101]  = {};
  Double_t y1[101] = {};
  Double_t y2[101] = {};

  ifstream infile;

  infile.open("data.txt");
//Check if the file opened
     if(infile.fail())
    {
      cout << "file did not open" << endl;
      return 1;
    }
//Put the data from the file into arrays
   for(int num = 0; num< 100; num++)
      {
         infile >> x[num];
         infile >> y1[num];
         infile >> y2[num];
         ++num;
        }

         infile.close();

///////////////////////////////////////////////////////

// Set Range for Data
Int_t ymin = 0;
Int_t ymax = 10;

// Create graphs in ROOT
Int_t marg = 0.15;
TCanvas *canvas = new TCanvas("c", "c", 600, 600);
canvas->Divide( 1,2 );

// Linear Pad
canvas->cd(1);
cpPadSetup( marg, marg, marg, marg );
TGraph *gr1 = new TGraph(n, x, y1);
cpGraphStyle( gr1,kRed,2,20 );
TGraph *gr2 = new TGraph(n, x, y2);
cpGraphStyle( gr2,kBlue,2,20 );
gr1->SetMaximum(ymax);
gr1->SetMinimum(ymin);
gr1->Draw("ALPE");
gr2->Draw("LPE");
gr1->GetXaxis()->SetTitle("X");
gr1->GetYaxis()->SetTitle("Y");
gr1->SetTitle("Linear Data");

// Logarithmic Pad
canvas->cd(2);
cpPadSetup( marg, marg, marg, marg );
gPad->SetLogx();
gPad->SetLogy();
TGraph *gr3 = new TGraph(n, x, y1);
cpGraphStyle( gr3,kRed,2,20 );
TGraph *gr4 = new TGraph(n, x, y2);
cpGraphStyle( gr4,kBlue,2,20 );
gr3->SetMaximum(ymax);
gr3->SetMinimum(ymin);
gr3->Draw("ALPE");
gr4->Draw("LPE");
gr3->GetXaxis()->SetTitle("log(X)");
gr3->GetYaxis()->SetTitle("log(Y)");
gr3->SetTitle("Logarithmic Data");

// Save
canvas->SaveAs("graphs.png");

return 0;
}
