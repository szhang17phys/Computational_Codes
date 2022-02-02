//--szhang; Oct 30, 2021------

#include "TCanvas.h"
#include "TROOT.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TCanvas.h"            //(x, y) graphs
#include "TLegend.h" 
#include "TArrow.h"
#include "TLatex.h"

#include <time.h>//to use clock---
#include "comp_set2.hpp"
using namespace std;

void format_h(TH1F* h, int linecolor){
    h->SetLineWidth(3);
    h->SetLineColor(linecolor);
}

//===Main==================================
void prob_1_trans(){
    double startTime1 = clock();//start timing---

    Int_t N = 15;
    Double_t x[N];//x is array for number of integrals---
    Double_t dev_1[N];//used to record the deviation---
    Double_t dev_2[N];
    Double_t dev_3[N];
    Double_t sim1 = 0.0;//to record value of simulation---
    Double_t sim2 = 0.0;
    Double_t sim3 = 0.0;

//---initialize------
    for(int i=0; i<N; ++i){
	x[i] = i+1;
        sim1 = triple(f1_trans, (i+1));
        sim2 = triple(f2_trans, (i+1));
	sim3 = triple(f3_trans, (i+1));
        dev_1[i] = fabs(sim1-1.0);//exp(-x^2)---
        dev_2[i] = fabs(sim2-1.0);//exp(-x))---
	dev_3[i] = fabs(sim3-0.5);//exp(-x)*cos(x)---
    }
    
    double endTime1 = clock();//end timing---
//    cout<<"test: "<<triple(f3_trans, 4)<<endl;
//---Graph------
    double startTime2 = clock();//start timing---    

    TGraph *dif_1 = new TGraph(15, x, dev_1);//15 CANNOT variable---
    TGraph *dif_2 = new TGraph(15, x, dev_2);
    TGraph *dif_3 = new TGraph(15, x, dev_3);

    dif_1->SetMarkerColor(kRed);
    dif_1->SetMarkerStyle(21);
    dif_1->SetLineColor(kRed);
    dif_1->SetLineWidth(2);
    dif_1->GetXaxis()->SetTitle("Log(N)");
    dif_1->GetYaxis()->SetTitle("Log(|Num-Ana|)");
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

    TLegend *leg = new TLegend(.6, .6, .9, .9);
    leg->AddEntry(dif_1, "exp(-x^2)");
    leg->AddEntry(dif_2, "exp(-x)");
    leg->AddEntry(dif_3, "exp(-x)*cos(x)");
    leg->Draw();
 
    double endTime2 = clock();//end timing---
     
    cout<<"Time for calculation is: "<<(double)(endTime1-startTime1)/CLOCKS_PER_SEC<<"s"<<endl;
    cout<<"Time for drawing is:     "<<(double)(endTime2-startTime2)/CLOCKS_PER_SEC<<"s"<<endl;

}
 
#ifndef __CINT__
int main(){
    prob_1_trans();
}
#endif


