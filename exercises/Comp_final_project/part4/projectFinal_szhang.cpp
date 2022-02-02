//--szhang; Nov 29, 2021------
//this file is about part 4 of the project---
//the method used in this file to set x range is cumbersome---
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <cstdlib>//random related---
#include <ctime>//time related---

#include "TROOT.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"            //(x, y) graphs
#include "TLegend.h" 
#include "TArrow.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"

#include "cpROOT.hpp"//very important---
#include "cpRandom.hpp"//random---

#define PI 3.14159265359

using namespace std;

void format_h(TH1F* h, int linecolor){
    h->SetLineWidth(3);
    h->SetLineColor(linecolor);
}

//===Main==================================
void projectFinal(){
    srand((unsigned)time(NULL));
    cpRandMT randomC(rand());//for capture---
    cpRandMT randomD(rand());//for distance---
    cpRandMT randomT(rand());//for theta---
    cpRandMT randomP(rand());//for phi---

    int num = 10;//number of neutrons---
    int steps = 20;//expected maximal steps---
    double lambda = 0.5;//lambda---
    double thick = 1.0;//thickness---
    double pc = 0.1;//capture prob for each neutron---
    double length_x = 0.0;//propagation length along x---
    double length_y = 0.0;//propagation length along y---
    double length_z = 0.0;//propagation length along z---
    double pos_x[num][steps] = {0.0};//position x of each step---
    double pos_y[num][steps] = {0.0};//position y of each step---
    double pos_z[num][steps] = {0.0};//z position in each step---
    int valid_xyz[num] = {0};//valid steps of x y z---
    double dis = 0.0;//distance along theta for each step---
    double cosTheta = 0.0;//cos(theta)---
    double phi = 0.0;//phi---
    
    TGraph *rout_xz[num];
    TGraph *rout_yz[num];
    TMultiGraph *mg_xz = new TMultiGraph();
    TMultiGraph *mg_yz = new TMultiGraph();

    for(int i=0; i<num; ++i){
	length_x = 0.0;
	length_y = 0.0;
        length_z = 0.0; 
        for(int j=0; j<1000000; ++j){//100000: Maximal steps---
	    dis = -lambda*log(randomD.randDouble());
	    cosTheta = 1 - 2*randomT.randDouble();
	    phi = 2*PI*randomP.randDouble();
	    length_z += dis*cosTheta;
	    length_x += dis*sqrt(1-cosTheta*cosTheta)*cos(phi);
	    length_y += dis*sqrt(1-cosTheta*cosTheta)*sin(phi);
	   
            pos_x[i][j+1] = length_x;
            pos_y[i][j+1] = length_y;
            pos_z[i][j+1] = length_z;

	    if(randomC.randDouble() < pc){//Capture---
		cout<<"-----------------------"<<endl;
		cout<<"Total Steps: "<<(j+1)<<endl;
		cout<<"-----------------------"<<endl;
		break;
	    }
	    if(length_z < 0 || length_z > thick){//throw away those events first!---
		cout<<"-----------------------"<<endl;
		cout<<"Total Steps: "<<(j+1)<<endl;
		cout<<"-----------------------"<<endl;
		break;
	    }
	}
	for(int m=0; m<steps; ++m){   
	    cout<<"x: "<<pos_x[i][m]<<endl;
	    cout<<"y: "<<pos_y[i][m]<<endl;
	    cout<<"z: "<<pos_z[i][m]<<endl;
	    if(fabs(pos_x[i][m+1])<pow(10, -20) && pos_y[i][m+1]==0 && pos_z[i][m+1]==0){//reset range---
		cout<<"======================"<<endl;
		cout<<"Valid Steps: "<<m<<endl;
		cout<<"======================"<<endl;
		valid_xyz[i] = m;
		break;
	    }
	}
    }

//===test======   
    cout<<"===---Valid Steps---==="<<endl;
    for(int i=0; i<num; ++i){
	cout<<"------"<<(i+1)<<" : "<<valid_xyz[i]<<"------"<<endl;
        for(int j= valid_xyz[i]+1; j<steps; ++j){//make left points the same as the final vaild points---
            pos_x[i][j] = pos_x[i][valid_xyz[i]]; 
            pos_y[i][j] = pos_y[i][valid_xyz[i]];
            pos_z[i][j] = pos_z[i][valid_xyz[i]];
	}
   
        for(int k=0; k<steps; ++k){
	    cout<<"x: "<<pos_x[i][k]<<endl;
	    cout<<"y: "<<pos_y[i][k]<<endl;
	    cout<<"z: "<<pos_z[i][k]<<endl;
	}
    }	

    for(int i=0; i<num; ++i){
	cout<<"Valide Steps: "<<valid_xyz[i]<<endl;
    }

    double beauty_z[5] = {-0.2, -0.2, 1.2, 1.2, -0.2};
    double beauty_x[5] = {-1.0, 1.0, 1.0, -1.0, -1.0}; 
    TGraph *beauty = new TGraph(5, beauty_x, beauty_z);
    mg_xz->Add(beauty);
    mg_yz->Add(beauty);
    


//---fill the graph------
    for(int i=0; i<num; ++i){
        rout_xz[i] = new TGraph(num, &pos_x[i][0], &pos_z[i][0]);
        rout_xz[i]->GetXaxis()->SetTitle("x");
        rout_xz[i]->GetYaxis()->SetTitle("z");
	rout_xz[i]->SetLineWidth(2);
	rout_xz[i]->SetLineColor(i+1);
	mg_xz->Add(rout_xz[i]);
        mg_xz->GetXaxis()->SetTitle("x");
        mg_xz->GetYaxis()->SetTitle("z");
	mg_xz->GetXaxis()->SetRangeUser(-1.0, 1.0);
	mg_xz->GetYaxis()->SetRangeUser(-0.2, 1.2);

        rout_yz[i] = new TGraph(num, &pos_y[i][0], &pos_z[i][0]);
        rout_yz[i]->GetXaxis()->SetTitle("y");
        rout_yz[i]->GetYaxis()->SetTitle("z");
	rout_yz[i]->GetXaxis()->SetRangeUser(-1.0, 1.0);
	rout_yz[i]->GetYaxis()->SetRangeUser(-0.2, 1.2);
	rout_yz[i]->SetLineWidth(2);
	rout_yz[i]->SetLineColor(i+1);
	mg_yz->Add(rout_yz[i]);
	mg_yz->GetXaxis()->SetTitle("y");
        mg_yz->GetYaxis()->SetTitle("z");
	mg_yz->GetXaxis()->SetRangeUser(-1.0, 1.0);
	mg_yz->GetYaxis()->SetRangeUser(-0.2, 1.2);
    }
		
//---TFile------
    TFile *file_1 =TFile::Open("out.root", "new");

    TCanvas *canv1 = new TCanvas("", "", 700, 500);
    mg_xz->Draw("ALP");
    canv1->Write();

    TCanvas *canv2 = new TCanvas("", "", 700, 500);
    mg_yz->Draw("ALP");
    canv2->Write();


    file_1->Close();



}
 
#ifndef __CINT__
int main(){
    projectFinal();
}
#endif


