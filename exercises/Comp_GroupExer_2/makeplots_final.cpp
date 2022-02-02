// Group Problem 2
// 9/15/2021

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "cpROOT.hpp"//very important for using some root variables in c++
#include <fstream> //read data from file && output data to file


using namespace std;

int main(int argc, char *argv[]){
   Double_t x[100]  = {};
   Double_t y1[100] = {};
   Double_t y2[100] = {};

   ifstream infile;
   infile.open("data.txt");

//Check if the file opened
   if(infile.fail()){
       cout << "file did not open" << endl;
       return 1;
   }
//Put the data from the file into arrays
    for(int num = 0; num < 100; num++){
        infile >> x[num];
        infile >> y1[num];
        infile >> y2[num];
      //   cout << "Y2 araay: "<< y2[num] <<endl; // prints out the y2 values from the array
       
    }

    infile.close();

//===Draw the plots================================
    Int_t n = 100;
// Create and divide canvas---
    TCanvas *canvas = new TCanvas("c", "c", 1200, 500);
    canvas->Divide(2, 1);
    canvas->SetGrid();
//Details about curves---
    TGraph *gr1 = new TGraph(n, x, y2);
    gr1->SetLineColor(kBlue);
    gr1->SetLineWidth(2);
    gr1->GetXaxis()->SetTitle("x");
    gr1->GetYaxis()->SetTitle("f(x)");    
//    gr1->SetTitle("Linear Axis");

    TGraph *gr2 = new TGraph(n, x, y1);
    gr2->SetLineColor(kRed);
    gr2->SetLineWidth(2);

    canvas->cd(1);
    canvas->SetGrid();
    gr1->Draw("AC");
    gr2->Draw("C");

// Logarithmic Pad
    canvas->cd(2);
    canvas->SetGrid();
    gPad->SetLogy();
//Draw---
//    gr1->SetTitle("Logarithmic Axis");
    gr1->Draw("AC");
    gr2->Draw("C");


// Save
    canvas->SaveAs("graphs.png");

return 0;
}
