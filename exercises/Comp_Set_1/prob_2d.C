//--szhang; Oct 2, 2021------
//For d part of problem 2------

#include "TCanvas.h"
#include "TROOT.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h" 
#include "TArrow.h"
#include "TLatex.h"

#include "comp_set1.hpp"
using namespace std;
const double step = 0.0004;//time of each step---
const int repeat = 5750;
const double x_0 = 0.0;//initial value of x---
const double y_0 = 0.0;//initial value of y---

void format_h(TH1F* h, int linecolor){
    h->SetLineWidth(3);
    h->SetLineColor(linecolor);
}

//===Main==================================
void prob_2d(){
   
    int anl_num = 5;//number of angles we consider---
    double angle[anl_num];//store angles---
    double v_x0[anl_num];//vx---
    double v_y0[anl_num];//vy---
    double v_x0_a = 0.0;//vx of analytical process---
    double v_y0_a = 0.0;//vy of analytical process---
    ball point;//record each step of the motion---
    double x[anl_num][repeat];//record x---
    double y[anl_num][repeat];//record y---
    double x_a[repeat];//analytic evolution of x---
    double y_a[repeat];//analytic evolution of y---
    double dxdt_mid = 0.0;//(dxdt)_{i+1}---
    double dydt_mid = 0.0;//(dydt)_{i+1}---

//---Initialize------
    x_a[0] = x_0;
    y_a[0] = y_0;
    for(int i=0; i<anl_num; ++i){
        angle[i] = (42.0 + 0.5*i)*PI/180;//set steps of angle---
	v_x0[i] = 15*cos(angle[i]);
        v_y0[i] = 15*sin(angle[i]);
        x[i][0] = x_0;
	y[i][0] = y_0;
        cout<<"Angles: "<<angle[i]<<endl;
    }
    v_x0_a = 15*cos(PI/4.0);
    v_y0_a = 15*sin(PI/4.0);

//--Runge-Kutta Algorithm------
    for(int j=0; j<anl_num; ++j){    
        point.x = x_0;//initialize of x---
        point.y = y_0;//initialize of y---
        point.dxdt = v_x0[j];//initialize of dxdt---
        point.dydt = v_y0[j];//initialize of dydt---
        for(int i=1; i<repeat; ++i){
            point.ddxddt = func_x(point.dxdt, point.dydt);
            point.ddyddt = func_y(point.dxdt, point.dydt);
            dxdt_mid = point.dxdt + step*point.ddxddt;
            dydt_mid = point.dydt + step*point.ddyddt; 
            point.x += step/2.0*(point.dxdt + dxdt_mid);
            point.y += step/2.0*(point.dydt + dydt_mid);
    	    point.dxdt = dxdt_mid;
       	    point.dydt = dydt_mid;
            x[j][i] = point.x;//store---
            y[j][i] = point.y;
        }
    }

//analytical result------
    for(int i=0; i< repeat; ++i){
        x_a[i] = v_x0_a*step*i;
        y_a[i] = v_y0_a*step*i - g*pow(step*i, 2)/2.0;
    }

//---TEST----------------------------------------
//    double te[7] = {1, 2, 3, 4, 5, 6, 7};
    cout<<"====================="<<endl;
    for(int i=0; i<anl_num; ++i){
        for(int j=0; j<repeat; ++j){
            cout<<"x"<<j<<": "<<x[i][j]<<"; y"<<j<<": "<<y[i][j]<<endl;
        }
        cout<<"-------------"<<endl;
    }
    cout<<"====================="<<endl;
//-----------------------------------------------

//---Graph------
    TMultiGraph *mg = new TMultiGraph();//draw multigraph---
    TGraph *curve[anl_num];
    for(int i=0; i<anl_num; ++i){ 
        curve[i] = new TGraph(repeat, &x[i][0], &y[i][0]);
    }//& here is of key important!!!---
    TGraph *curve_a = new TGraph(repeat, x_a, y_a);

    curve_a->GetXaxis()->SetTitle("x(t)");
    curve_a->GetYaxis()->SetTitle("y(t)");
    curve_a->SetLineColor(2);
    curve_a->SetLineWidth(2);
    curve_a->SetMarkerColor(6);
    curve_a->SetMarkerStyle(24);
    curve_a->SetMarkerSize(0.6);
    for(int i=0; i<anl_num; ++i){
        curve[i]->SetMarkerColor(i+3);
        curve[i]->SetMarkerSize(0.7);
        curve[i]->SetMarkerStyle(21);
    }
    for(int i=0; i<anl_num; ++i){
        mg->Add(curve[i]);
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 700, 500);
    gStyle->SetOptStat("nemr");

    curve_a->Draw("APC");
    for(int i=0; i<anl_num; ++i){
        curve[i]->Draw("P");//P means marker---
    }
//---Legend------
    TLegend *leg = new TLegend(.7, .5, .9, .9);
    leg->AddEntry(curve_a, "Analytical, 0 Resis");
    for(int i=0; i<anl_num; ++i){
        leg->AddEntry(curve[i], "RK, Angle= deg"); 
    }
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
    prob_2d();
}
#endif


