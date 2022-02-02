//--szhang; Sep 30, 2021------

#include "TCanvas.h"
#include "TROOT.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TCanvas.h"            //(x, y) graphs
#include "TLegend.h" 
#include "TArrow.h"
#include "TLatex.h"

#include "comp_set1.hpp"
using namespace std;
const double step = 0.1;//time of each step---
const int repeat = 23;
const double x_0 = 0.0;//initial value of x---
const double y_0 = 0.0;//initial value of y---
const double v_x0 = 15*sin(PI/4.0);//initial value of vx---
const double v_y0 = 15*cos(PI/4.0);//initial value of vy---

void format_h(TH1F* h, int linecolor){
    h->SetLineWidth(3);
    h->SetLineColor(linecolor);
}

//===Main==================================
void prob_2(){

    ball point;//record each step of the motion---
    double x[repeat];//record x---
    double y[repeat];//record y---
    double R[repeat];//range--- 
    double x_a[repeat];//analytic evolution of x---
    double y_a[repeat];//analytic evolution of y---
    double dxdt_mid = 0.0;//(dxdt)_{i+1}---
    double dydt_mid = 0.0;//(dydt)_{i+1}---
    double stepNum[repeat];//to record step numbers
    double rangeDif[repeat];//store range differences---
    double rangeRK = 0.0;//tranfer range vaule of RK---
    double rangeAn = 0.0;//tranfer range value of analytic---

//---Initialize------
    point.x = x_0;
    point.y = y_0;
    point.dxdt = v_x0;
    point.dydt = v_y0;
    x[0] = x_0;
    y[0] = y_0;
    x_a[0] = x_0;
    y_a[0] = y_0;
    stepNum[0] = 1.0;
    rangeDif[0] = 0.0;

    for(int i=1; i<repeat; ++i){
        point.ddxddt = func_x(point.dxdt, point.dydt);
        point.ddyddt = func_y(point.dxdt, point.dydt);
        dxdt_mid = point.dxdt + step*point.ddxddt;
        dydt_mid = point.dydt + step*point.ddyddt; 
        point.x += step/2.0*(point.dxdt + dxdt_mid);
        point.y += step/2.0*(point.dydt + dydt_mid);
	point.dxdt = dxdt_mid;
	point.dydt = dydt_mid;
        x[i] = point.x;//store---
        y[i] = point.y;
        x_a[i] = v_x0*step*i;
        y_a[i] = v_y0*step*i - g*pow(step*i, 2)/2.0;
    }

//---TEST----------------------------------------
    cout<<"====================="<<endl;
//-----------------------------------------------

//---Graph------
    TGraph *curve = new TGraph(repeat, x, y);
    TGraph *curve0 = new TGraph(repeat, x_a, y_a);

    curve->SetMarkerColor(4);
    curve->SetMarkerSize(0.8);
    curve->SetMarkerStyle(21);
    curve0->GetXaxis()->SetTitle("x(t)");
    curve0->GetYaxis()->SetTitle("y(t)");
    curve0->SetLineColor(2);
    curve0->SetLineWidth(2);
    curve0->SetMarkerColor(6);
    curve0->SetMarkerStyle(24);
    curve0->SetMarkerSize(0.6);

    TCanvas *c1 = new TCanvas("c1", "c1", 700, 500);
    gStyle->SetOptStat("nemr");

    curve0->Draw("APC");
    curve->Draw("P");//P means marker---

//---Legend------
    TLegend *leg = new TLegend(.5120378, .7829703, .7736582, .9008746);
    leg->AddEntry(curve, "2nd Runge-Kutta"); 
    leg->AddEntry(curve0, "Analytical, 0 Resis");
    leg->Draw();

//---dashed line------
   TLine *line = new TLine(0,0,25,0);
   line->SetLineStyle(2);
   line->SetLineWidth(2);
   line->Draw();

//---remark------
   TLatex *   tex = new TLatex(2.479554,-0.3091959,"y=0");
   tex->SetTextFont(42);
   tex->SetTextSize(0.04201681);
   tex->SetLineWidth(2);
   tex->Draw();

}
       
#ifndef __CINT__
int main(){
    prob_2();
}
#endif


