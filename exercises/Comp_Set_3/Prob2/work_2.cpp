//--szhang; Nov 22, 2021------

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <cstdlib>//random related---
#include <ctime>//time related---

#include "TROOT.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"//(x, y) graphs---
#include "TLegend.h" 
#include "TArrow.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TFile.h"

#include "cpROOT.hpp"//very important---
#include "cpRandom.hpp"//random---

#include "../comp_set3.hpp"//my own headfile---
using namespace std;
#define Anal 0.25*PI

void format_h(TH1F* h, int linecolor){
    h->SetLineWidth(3);
    h->SetLineColor(linecolor);
}

//===Main==================================
void work_2(){
    srand((unsigned)time(NULL));
    cpRandMT random(rand());

    double Num[8];//different points---
    double simu[8];//simulated values---
    double devi[8];//deviation bet simu and anal---
    double fluc[8];//statistical uncertainties---
    double temp = 0.0;//temporary sums of f/g---
    double tempSq = 0.0;//temporary sums of f^2/g^2---
    double tempfSq = 0.0;//temporary sums of f^2---
    double value_x = 0.0;//show distributions of x---
    double var_fp = 0.0;//variance of f/p---
    double var_f = 0.0;//variance of f---

    for(int i=0; i<8; ++i){
	Num[i] = 1.0*pow(10, (i+2));
    }

    for(int i=0; i<8; ++i){//different numbers---
	temp = 0.0;
        tempSq = 0.0;
	for(int j=0; j<Num[i]; ++j){
	    value_x = xDis(random.randDouble());
	    temp += 1.0*f(value_x)/p(value_x);
	    tempSq += 1.0*f(value_x)*f(value_x)/(p(value_x)*p(value_x));
//	    tempfSq += 1.0*f(value_x)*f(value_x)/p(value_x);//for part d---
        }
	simu[i] = temp*1.0/Num[i];
	devi[i] = simu[i] - Anal;
	temp = simu[i];
	tempSq = tempSq*1.0/Num[i];
//        tempfSq = tempfSq*1.0/Num[i];//for part d---
	fluc[i] = sqrt((tempSq - temp*temp)*1.0/Num[i]);

	cout<<"Num: "<<Num[i]<<"------------------"<<endl;
	cout<<"Simulated: "<<simu[i]<<endl;
	cout<<"Deviation: "<<devi[i]<<endl;    
	cout<<"Fluctuati: "<<fluc[i]<<endl;

/*	if(i==7){//for part d---
	    cout<<"========Part d========"<<endl;
	    cout<<"Total Points: "<<Num[7]<<endl;
	    var_fp = tempSq - temp*temp;
	    var_f = tempfSq - temp*temp;
	    cout<<"Variance of f/p: "<<var_fp<<endl;
	    cout<<"Variance of f  : "<<var_f<<endl;
	    cout<<"======================"<<endl;
	}
*/
    }



//---Graph------
    TFile *file_1 =TFile::Open("out.root", "new");
    TGraphErrors *value = new TGraphErrors(8, Num, devi, 0, fluc);
    TCanvas *canv = new TCanvas("", "", 700,500);
	
    value->SetMarkerColor(kRed);
    value->SetMarkerStyle(21);
    value->SetLineColor(kBlue);
    value->SetLineWidth(2);
    value->GetXaxis()->SetTitle("Num");
    value->GetYaxis()->SetTitle("(Simulated - Analytical");
    value->SetTitle("Deviation");  

    gPad->SetLogx();
    gPad->SetGrid();
//    value[i]->GetYaxis()->SetRangeUser(0.0, (anal[i]+1.5));//set range, must be here(after deciding the site)---
    value->Draw("APL");

    canv->SaveAs("results.png");
    canv->Write();
    file_1->Close();


}
 
#ifndef __CINT__
int main(){
    work_2();
}
#endif


