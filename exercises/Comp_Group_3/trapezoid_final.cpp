// Group Exercise 3 of Computational Physics---
//Brooks, Hagen, Zhang---
//Oct 12, 2021---

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>

#include "TMultiGraph.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "cpROOT.hpp"//very important---

#define PI 3.1415926//accurary?----
using namespace std;

double poly(double x){
	return x*((x*x)-1); //First Integral---
}

double trig(double x){
    return 2*PI*cos(2*PI*x);
}

void trape(double (*func)(double), double *array, int size){
    double point = 0.0;//to store new points---
    double sum = 0.0;//to store function value of new points---
    
    array[0] = 0.5*(func(0.0)+func(1.0));

    for(int i=1; i<size; ++i){
	for(int j=0; j<pow(2, (i-1)); ++j){
	    point = (2*j+1)*1.0/pow(2, i);
	    sum += func(point);
        }
        array[i] = 0.5*array[i-1] + 1.0*sum/pow(2, i);
	sum = 0;
    }    

}


int main(int argc, char *argv[]){

    Double_t N[20];//to store N---
    Double_t y1[20];
    Double_t y2[20];

    for(int i=0; i<20; ++i){//initialize---
        N[i] = i;
    }

    trape(poly, y1, 20);
    trape(trig, y2, 20);

    for(int i=0; i<20; ++i){//to calculate the differences---
	y1[i] = fabs(y1[i]+0.25);
	y2[i] = fabs(y2[i]);
    }

//===test=====================
/*    for(int i=0; i<20; ++i){
	cout.setf(ios::scientific);
	cout.precision(12);
	cout<<y2[i]<<endl;
    }
*/

//===Draw the plots======
    TCanvas *canvas = new TCanvas("c", "c", 700, 500);

    gPad->SetLogy();

    TGraph *gr1 = new TGraph(20, N, y1);
    TGraph *gr2 = new TGraph(20, N, y2);
    gr2->SetLineColor(kBlue);
    gr2->SetLineWidth(2);
    gr2->SetMarkerStyle(21);
    gr2->SetMarkerColor(kBlue);
    gr2->GetXaxis()->SetTitle("Log(N)");
    gr2->GetYaxis()->SetTitle("Log(|Num-Ana|)");    
    gr1->SetLineColor(kRed);
    gr1->SetLineWidth(2);
    gr1->SetMarkerStyle(21);
    gr1->SetMarkerColor(kRed);

    canvas->SetGrid();
    gr2->Draw("ACP");
    gr1->Draw("CP");

/*---Legend does not work---
    TLegend *leg = new TLegend(.65, .9, .75, .9);
    leg->AddEntry(gr1, "Polynominal");
    leg->AddEntry(gr2, "Trigonometric");
    leg->Draw("same");
*/
    canvas->SaveAs("graphs.png");//save the figure---

    return 0;
}


