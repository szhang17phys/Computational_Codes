//--szhang; Sep 29, 2021------
/*
#include "TCanvas.h"
#include "TROOT.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TCanvas.h"            //(x, y) graphs
#include "TLegend.h" 
#include "TArrow.h"
#include "TLatex.h"
*/
#include "comp_set1.hpp"
using namespace std;
const double step = 0.0002;
const int repeat = 40000;

void format_h(TH1F* h, int linecolor){
    h->SetLineWidth(3);
    h->SetLineColor(linecolor);
}

//===Main==================================
void prob_1(){

    pend mot;//record each step of the motion---
    double x[repeat];//x is time.---
    double y[repeat];//y is theta---
    double E[repeat];//E is the energy--- 
    double anal[repeat];//analytical angle evolution---   
    double Amp = 0.1;//amplitude---

    mot.s = Amp;
    mot.dsdt = 0;
    mot.ddsddt = f(mot.s);
    x[0] = 0.0;
    y[0] = mot.s;
    E[0] = mot.energy();
    anal[0] = Amp;

    for(int i=1; i<repeat; ++i){
        mot.s += mot.dsdt*step;
        mot.dsdt += mot.ddsddt*step;
        mot.ddsddt = f(mot.s);
        x[i] = i*step;
        y[i] = mot.s;
//        E[i] = Energy(mot.s, mot.dsdt);
        E[i] = mot.energy();
        anal[i] = theta(Amp, x[i]);
    }

//    cout<<"Final time is: "<<x[repeat-1]<<endl;//check time---

//---TEST----------------------------------------
//    double te[7] = {1, 2, 3, 4, 5, 6, 7};
//    cout<<"====================="<<endl;
//    cout<<"Energy Dev: "<<Dev(E, repeat)<<endl;
    cout<<"====================="<<endl;
//-----------------------------------------------

//to get the period------
    double y1 = 0.0;
    double y2 = 0.0;
    double y3 = 0.0;
    double period[10];//record all period points---
    int num = 0;//number of period points---
    double final_period = 0.0;

    for(int i=0; i<repeat; ++i){
        y1 = y[i];
        y2 = y[i+1];
        y3 = y[i+2];
   
        if(y2>y1 && y2>y3){
    	    period[num] = x[i+1];
            num += 1;
            cout<<"Peorid points: t =  "<<period[num-1]<<endl;
        }
    }
    final_period = 1.0*period[num-1]/num;
    cout<<"The period of the pendulum is: "<<final_period<<endl;
    cout<<"The difference to small-angle period is: "<<abs(final_period-1.0)<<endl;

//---Graph------
    TGraph *theta = new TGraph(repeat, x, y);
    TGraph *theta0 = new TGraph(repeat, x, anal);
    TGraph *energy = new TGraph(repeat, x, E);

    theta->SetMarkerColor(4);
    theta->SetMarkerSize(0.2);
    theta->SetMarkerStyle(21);
    theta->GetXaxis()->SetTitle("Time");
    theta->GetYaxis()->SetTitle("Theta");
    theta0->SetLineColor(2);
    theta0->SetLineWidth(2);
    energy->SetLineColor(2);
    energy->SetLineWidth(2);
    energy->SetMarkerStyle(21);
    energy->SetMarkerSize(0.2);
    energy->SetMarkerColor(4);
    energy->GetXaxis()->SetTitle("Time");
    energy->GetYaxis()->SetTitle("Energy");

    TCanvas *c1 = new TCanvas("c1", "c1", 1300, 500);
    c1->Divide(2, 1);

    c1->cd(1);
    gStyle->SetOptStat("nemr");
    theta->Draw("AP");//P means marker---
    theta0->Draw("C");//draw curve---

    TLegend *leg = new TLegend(.5120378, .7829703, .7736582, .9008746);
    leg->AddEntry(theta, "General Angle Solution"); 
    leg->AddEntry(theta0, "Small Angle Approximation");
    leg->Draw();

//    TCanvas *c2 = new TCanvas("c2", "c2", 200, 10, 700, 500);
    c1->cd(2);
    gStyle->SetOptStat("nemr");
    energy->Draw("AP");


}
       
#ifndef __CINT__
int main(){
    prob_1();
}
#endif


