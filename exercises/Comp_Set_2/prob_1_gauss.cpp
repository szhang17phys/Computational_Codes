//--szhang; Oct 30, 2021------

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <time.h>//to timing---

#include "TROOT.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TCanvas.h"            //(x, y) graphs
#include "TLegend.h" 
#include "TArrow.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TFile.h"

#include "cpROOT.hpp"//very important---

#include "comp_set2.hpp"//my own headfile---
#include "/N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/library/include/cpQuadrature.hpp"
#include "/N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/library/include/cpUtils.hpp"
#include "/N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/library/include/cpFunctions.hpp"

using namespace std;

void format_h(TH1F* h, int linecolor){
    h->SetLineWidth(3);
    h->SetLineColor(linecolor);
}

//===Main==================================
void prob_1_gauss(){
    double startTime1 = clock();

    double num[8];//store number of points used in the simu---
    double dev1[8];//To store the deviation for first function---
    double dev2[8];
    double dev3[8];    
    vector<double> x1;//x1: initial empty vector---
    vector<double> w1;
    vector<double> x2;//x2 are for the Laguerre--- 
    vector<double> w2;
    double sum1 = 0.0;//to store temporary sum---
    double sum2 = 0.0;
    double sum3 = 0.0;

    for(int i=0; i<8; ++i){
        num[i] = 2*(i+1);
        x1.push_back(1.0);//make sure volume increased---
        x1.push_back(1.0);
        w1.push_back(1.0);
        w1.push_back(1.0);
        x2.push_back(1.0);
        x2.push_back(1.0);
        w2.push_back(1.0);
        w2.push_back(1.0);

        GaussHermite(1.0, 1.0, x1, w1);//key---
        GaussLaguerre(0.0, 1.0, x2, w2);//key---


        for(int j=0; j<num[i]; ++j){
    	    sum1 += 1.0/sqrt(PI)*w1[j];//outputing vector---
	    sum2 += 1.0*w2[j];
            sum3 += w2[j]*cos(1.0*x2[j]);
        }

        dev1[i] = fabs(sum1-1.0);
	dev2[i] = fabs(sum2-1.0);
        dev3[i] = fabs(sum3-0.5);

//        cout<<"Sum 1: "<<sum1<<endl;
//        cout<<"Sum 2: "<<sum2<<endl;
//        cout<<"Sum 3: "<<sum3<<endl;

        sum1 = 0.0;
	sum2 = 0.0;
	sum3 = 0.0;
    }

    double endTime1 = clock();

//---Graph------
    double startTime2 = clock();

    TGraph *dif_1 = new TGraph(8, num, dev1);//8 CANNOT variable---
    TGraph *dif_2 = new TGraph(8, num, dev2);
    TGraph *dif_3 = new TGraph(8, num, dev3);

    dif_1->SetMarkerColor(kRed);
    dif_1->SetMarkerStyle(21);
    dif_1->SetLineColor(kRed);
    dif_1->SetLineWidth(2);
    dif_1->GetXaxis()->SetTitle("N");
    dif_1->GetYaxis()->SetTitle("Log(|Num-Ana|)");
    dif_1->GetYaxis()->SetRangeUser(1.0*pow(10, -20), 0.1);//Set Range---
    dif_2->SetMarkerColor(kBlue);
    dif_2->SetMarkerStyle(22);
    dif_2->SetLineColor(kBlue);
    dif_2->SetLineWidth(2);
    dif_3->SetMarkerColor(kGreen);
    dif_3->SetMarkerStyle(23);
    dif_3->SetLineColor(kGreen);
    dif_3->SetLineWidth(2);


    TCanvas *c1 = new TCanvas("c1", "c1", 700, 500);
    gStyle->SetOptStat("nemr");
    c1->SetGrid();
    gPad->SetLogy();

    dif_1->Draw("APC");//P means marker---
    dif_2->Draw("PC");
    dif_3->Draw("PC");

    TLegend *leg = new TLegend(.6, .65, .9, .9);
    leg->AddEntry(dif_1, "exp(-x^2)");
    leg->AddEntry(dif_2, "exp(-x)");
    leg->AddEntry(dif_3, "exp(-x)*cos(x)");
    leg->Draw();

//    c1->SaveAs("graphs.png");

    TFile *file_1 =TFile::Open("out.root", "new");
    c1->Write();        
    file_1->Close();
 
   double endTime2 = clock();

    cout<<"Time for calculation is: "<<(double)(endTime1-startTime1)/CLOCKS_PER_SEC<<"s"<<endl; 
    cout<<"Time for drawing is:     "<<(double)(endTime2-startTime2)/CLOCKS_PER_SEC<<"s"<<endl;

}
 
#ifndef __CINT__
int main(){
    prob_1_gauss();
}
#endif


