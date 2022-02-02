//--szhang; Nov 29, 2021------
//this file is about part 1 of the project---
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

    double lambda = 0.1;//lambda---
    double thick = 1.0;//thickness---
    double pc = 0.1;//capture prob for each neutron---
    double length_x = 0.0;//propagation length along x---
    double length_y = 0.0;//propagation length along y---
    double length_r = 0.0;//length_r =sqrt(x^2+y^2)---
    double length_z = 0.0;//propagation length along z---
    double dis = 0.0;//distance along theta for each step---
    double cosTheta = 0.0;//cos(theta)---
    double phi = 0.0;//phi---
    int N_a = 0;//number of absorbed neutrons---
    int N_r = 0;//number of reflected neutrons---
    int N_t = 0;//number of transmitted neutrons---
    int N_total = 100000000;//total neutrons---
    
    TH1F *dis_z = new TH1F("dis_z", "Distribution of Path Length", 101, -0.005, 1.005);//F: float---
    TH3D *dis_xyz = new TH3D("dis_xyz", "Distributions of absorbed neutrons", 201, -1.005, 1.005, 201, -1.005, 1.005, 101, -0.005, 1.005);//D: double---
    TH2D *dis_tran_xy = new TH2D("dis_tran_xy", "Distributions of transmitted neutrons", 201, -1.005, 1.005, 201, -1.005, 1.005);
    TH2D *speed_xy = new TH2D("speed_xy", "", 201, -0.005, 2.005, 101, -0.005, 1.005);

    for(int i=0; i<N_total; ++i){
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

	    if(randomC.randDouble() < pc){//Capture---
		N_a += 1;
		dis_xyz->Fill(length_x, length_y, length_z);
	        dis_z->Fill(length_z);
		break;
	    }
	    if(length_z < 0){//reflected neutrons---
		N_r += 1;
		break;
	    }
	    if(length_z > thick){
		N_t += 1;
		dis_tran_xy->Fill(length_x, length_y);
		length_r = sqrt(length_x*length_x + length_y*length_y);//calculate r here---
		speed_xy->Fill(length_r, cosTheta);
		break;
	    }
	}
    }



//---PRINT------
    cout<<"========================="<<endl;
    cout<<"Absorbed   : "<<N_a<<"; P_abs = "<<(1.0*N_a/N_total)<<endl;
    cout<<"Reflected  : "<<N_r<<"; P_ref = "<<(1.0*N_r/N_total)<<endl;
    cout<<"Transmitted: "<<N_t<<"; P_tra = "<<(1.0*N_t/N_total)<<endl;
    cout<<"Total Num  : "<<(N_a+N_r+N_t)<<endl;
    cout<<"========================="<<endl;

		
//---TFile------
    TFile *file_1 =TFile::Open("out.root", "new");

    dis_z->GetXaxis()->SetTitle("Path Length");
    dis_z->GetYaxis()->SetTitle("Counts");
//    dis_z->SetTitle("lambda = 0.01");

    dis_xyz->GetXaxis()->SetTitle("x");
    dis_xyz->GetYaxis()->SetTitle("y");
    dis_xyz->GetZaxis()->SetTitle("z");
    dis_tran_xy->GetXaxis()->SetTitle("x");
    dis_tran_xy->GetYaxis()->SetTitle("y");
    speed_xy->GetXaxis()->SetTitle("sqrt(x^2+y^2)");
    speed_xy->GetYaxis()->SetTitle("cos(theta)");
   
    dis_z->Write();
    dis_xyz->Write();
    dis_tran_xy->Write();
    speed_xy->Write();

    file_1->Close();



}
 
#ifndef __CINT__
int main(){
    projectFinal();
}
#endif


