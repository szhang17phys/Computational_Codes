//--szhang; Nov 29, 2021------
//part 2 of the project---
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

    double lambda[6] = {0.01, 0.05, 0.1, 0.5, 1.0, 10.0};//lambda---
    double thick = 1.0;//thickness---
    double pc = 0.1;//capture prob---
    double length_x = 0.0;//propagation length along x---
    double length_y = 0.0;//propagation length along y---
    double length_r = 0.0;//length_r =sqrt(x^2+y^2)---
    double length_z = 0.0;//propagation length along z---
    double dis = 0.0;//distance along theta for each step---
    double cosTheta = 0.0;//cos(theta)---
    double phi = 0.0;//phi---
    int N_a[6] = {0};//number of absorbed neutrons---
    int N_r[6] = {0};//number of reflected neutrons---
    int N_t[6] = {0};//number of transmitted neutrons---
    int N_total = 100000000;//total neutrons---
    
    TH1F *dis_z[6];
    TH3D *dis_xyz = new TH3D("dis_xyz", "Distributions of absorbed neutrons", 201, -1.005, 1.005, 201, -1.005, 1.005, 101, -0.005, 1.005);//D: double---
    TH2D *dis_tran_xy = new TH2D("dis_tran_xy", "Distributions of transmitted neutrons", 201, -1.005, 1.005, 201, -1.005, 1.005);
    TH2D *speed_xy = new TH2D("speed_xy", "", 201, -0.005, 2.005, 101, -0.005, 1.005);

    for(int i=0; i<6; ++i){
        dis_z[i] = new TH1F("dis_z", "dis_z", 201, -0.0025, 1.0025);//F: float---
    }

    for(int s=0; s<6; ++s){
        for(int i=0; i<N_total; ++i){
	    length_x = 0.0;
	    length_y = 0.0;
	    length_z = 0.0;
            for(int j=0; j<1000000; ++j){//100000: Maximal steps---
	        dis = -lambda[s]*log(randomD.randDouble());
	        cosTheta = 1 - 2*randomT.randDouble();
	        phi = 2*PI*randomP.randDouble();
	        length_z += dis*cosTheta;
	        length_x += dis*sqrt(1-cosTheta*cosTheta)*cos(phi);
	        length_y += dis*sqrt(1-cosTheta*cosTheta)*sin(phi);

	        if(randomC.randDouble() < pc){//Capture---
		    N_a[s] += 1;
	            dis_z[s]->Fill(length_z);
		    break;
	        }
	        if(length_z < 0){//reflected neutrons---
		    N_r[s] += 1;
		    break;
	        }
	        if(length_z > thick){
	 	    N_t[s] += 1;
		    break;
	        }
       	    }
        }
        dis_z[s]->GetXaxis()->SetTitle("Path Length");
        dis_z[s]->GetYaxis()->SetTitle("Counts");
	dis_z[s]->SetLineColor(s+2);
	dis_z[s]->SetLineWidth(2);
    }


//---PRINT------
    for(int i=0; i<6; ++i){
        cout<<"========================="<<endl;
	cout<<"#lambda = "<<lambda[i]<<": "<<endl;
        cout<<"Absorbed   : "<<N_a[i]<<"; P_abs = "<<(1.0*N_a[i]/N_total)<<endl;
        cout<<"Reflected  : "<<N_r[i]<<"; P_ref = "<<(1.0*N_r[i]/N_total)<<endl;
        cout<<"Transmitted: "<<N_t[i]<<"; P_tra = "<<(1.0*N_t[i]/N_total)<<endl;
        cout<<"Total Num  : "<<(N_a[i]+N_r[i]+N_t[i])<<endl;
//    cout<<"========================="<<endl;
    }
		
//---TFile & Drawing------
    TFile *file_1 =TFile::Open("out.root", "new");
    TCanvas *canv = new TCanvas("","", 700, 500);
    
    dis_z[0]->Draw();
    for(int i=1; i<6; ++i){
	dis_z[i]->Draw("same");
    }

    TLegend *leg = new TLegend(.7, .5, .9, .9);
    leg->AddEntry(dis_z[0], "#lambda = 0.01");
    leg->AddEntry(dis_z[1], "#lambda = 0.05");
    leg->AddEntry(dis_z[2], "#lambda = 0.10");
    leg->AddEntry(dis_z[3], "#lambda = 0.50");
    leg->AddEntry(dis_z[4], "#lambda = 1.00");
    leg->AddEntry(dis_z[5], "#lambda =10.00");
    leg->Draw();

    canv->Write();
    file_1->Close();



}
 
#ifndef __CINT__
int main(){
    projectFinal();
}
#endif


