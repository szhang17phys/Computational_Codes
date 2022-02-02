//--szhang; Oct 30, 2021------

#include "TCanvas.h"
#include "TROOT.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TCanvas.h"            //(x, y) graphs
#include "TLegend.h" 
#include "TArrow.h"
#include "TLatex.h"
#include "TH1F.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>
#include <time.h>//to use clock---

#include "cpROOT.hpp"

#include "../comp_set2.hpp"

#include "/N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/library/include/cpODEmodel.hpp"
#include "/N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/library/include/cpUtils.hpp"
#include "/N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/library/include/cpODEstep.hpp"
#include "/N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/library/include/cpInterpolate.hpp"

using namespace std;

void format_h(TH1F* h, int linecolor){
    h->SetLineWidth(3);
    h->SetLineColor(linecolor);
}

//===Main==================================
void prob_2(){
    Int_t N = 15;//number of integrals: 2^15---
    Int_t epNum = 15;//number of epsilon---
    Double_t ep[epNum];//different values of epsilon---
    Double_t epShow[epNum];//values of log_2(ep)---
    Double_t dev_1[epNum];//deviation for fixed step-size---
    Double_t dev_2[epNum];//for adaptive step method---
    Double_t sim1 = 0.0;//to record value of simulation---
    Double_t sim2 = 0.0;

//---initialize------
    for(int i=0; i<epNum; ++i){
	ep[i] = 1.0/pow(2, epNum-1-i);
        epShow[i] = -14+i;
        sim1 = simpson2(delta, ep[i], 100,  N);//[0,100] integral---
	dev_1[i] = fabs(sim1-1.0);
        cout<<"Epsilon: "<<ep[i]<<endl;
	cout<<"sim1: "<<sim1<<endl;
//	cout<<"fixed step: "<<dev_1[i]<<endl;
    }
  
    

//---Graph------

    TGraph *dif_1 = new TGraph(15, epShow, dev_1);//15 CANNOT variable---
//    TGraph *dif_2 = new TGraph(15, x, dev_2);

    dif_1->SetMarkerColor(kRed);
    dif_1->SetMarkerStyle(21);
    dif_1->SetLineColor(kRed);
    dif_1->SetLineWidth(2);
    dif_1->GetXaxis()->SetTitle("Log_2(epsilon)");
    dif_1->GetYaxis()->SetTitle("Log(|Num-Ana|)");
//    dif_2->SetMarkerColor(kBlue);
//    dif_2->SetMarkerStyle(22);
//    dif_2->SetLineColor(kBlue);
//    dif_2->SetLineWidth(2);


    TCanvas *c1 = new TCanvas("c1", "c1", 700, 500);
    gStyle->SetOptStat("nemr");
    c1->SetGrid();
    gPad->SetLogy();

    dif_1->Draw("APC");//P means marker---
//    dif_2->Draw("PC");

    TLegend *leg = new TLegend(.6, .6, .9, .9);
    leg->AddEntry(dif_1, "Fixed step-size");
//    leg->AddEntry(dif_2, "Adaptive step");
    leg->Draw();

    c1->SaveAs("graph.png");

}
 
#ifndef __CINT__
int main(){
    prob_2();
}
#endif


