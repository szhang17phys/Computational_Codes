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
#include "TCanvas.h"            //(x, y) graphs
#include "TLegend.h" 
#include "TArrow.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TFile.h"

#include "cpROOT.hpp"//very important---
#include "cpRandom.hpp"//random---

#include "comp_set3.hpp"//my own headfile---
using namespace std;

void format_h(TH1F* h, int linecolor){
    h->SetLineWidth(3);
    h->SetLineColor(linecolor);
}

//===Main==================================
void work_1(){
    srand((unsigned)time(NULL));
    cpRandMT random(rand());
    cpRandLCG lcgRand(rand(), 5, 3, 32);//for part c---

    int D[5] = {2, 3, 4, 5, 6};//dimensions---
    double Num[8];//number of points---
    double sqrtNum[8];//to record sqrt(Num)---
    double anal[5];//analytical value---
    double simu[5][8];//simulated value---
    double fluc[5][8];//errors of each simu---
    double accu[5][8];//differences bet simu and anal---
    double inver_Diff[5][8];//1/(simu-anal)---
    double axis[6];//x, y, z, w---
    int GoodNum = 0;//points inside boundaries---
    double tempSum = 0.0;//sum of squared values---

    for(int i=0; i<8; ++i){
	Num[i] = 1.0*pow(10, (i+2));
	sqrtNum[i] = 1.0*sqrt(Num[i]);
    }


    for(int i=0; i<5; ++i){//different dimensions---
	if(i > 2)//for part c---
	    break;
	if(i != 2)
	    continue;

	anal[i] = V_ana(i+2);
        for(int j=0; j<8; ++j){//different numbers---
	    GoodNum = 0;//set 0 each time---
	    for(int k=0; k<Num[j]; ++k){
		tempSum = 0.0;//set 0 for each cycle---
		for(int l=0; l<(i+2); ++l){//(i+2) axis---
//		    axis[l] = random.randDouble()*2-1; 
		    axis[l] = lcgRand.randDouble()*2-1;//part c---
                }
		for(int l=0; l<(i+2); ++l){//check boundary---
		    axis[l] = 1.0*axis[l]*axis[l];
		    tempSum += 1.0*axis[l];
    		}
		if(tempSum<1){
		    GoodNum += 1;
		}

            }
            simu[i][j] = 1.0*pow(2, (i+2))*GoodNum/Num[j];
	    accu[i][j] = simu[i][j] - anal[i];
	    fluc[i][j] = 1.0*GoodNum*(Num[j] - GoodNum);
	    fluc[i][j] = 1.0*sqrt(fluc[i][j])/Num[j];
	    fluc[i][j] = fluc[i][j]*pow(2, (i+2))*1.0/sqrt(Num[j]);
            inver_Diff[i][j] = 1.0/fabs(accu[i][j]);

	    cout<<"Dim: "<<(i+2)<<", num: "<<Num[j]<<"------"<<endl;
	    cout<<"Analytical: "<<anal[i]<<endl;
	    cout<<"Difference: "<<accu[i][j]<<endl;
	    cout<<"Fluctuatio: "<<fluc[i][j]<<endl;
	    cout<<"sqrt num  : "<<sqrtNum[j]<<endl;
	    cout<<"InverseDif: "<<inver_Diff[i][j]<<endl;
        }
        cout<<"==============================="<<endl;
    }



//---Graph------
    TGraphErrors *value[5];//simulated values for d=2,3,4...---
    TGraph *diff[5];//error for d=2,3,4...---
    TCanvas *canv[5];
    TFile *file_1 =TFile::Open("out.root", "new");

    for(int i=0; i<5; ++i){
	value[i] = new TGraphErrors(8, Num, &accu[i][0], 0, &fluc[i][0]);
        diff[i] = new TGraph(8, sqrtNum, &inver_Diff[i][0]);
        canv[i] = new TCanvas("", "", 1200, 500);
	canv[i]->Divide(2, 1);
	
        value[i]->SetMarkerColor(kRed);
        value[i]->SetMarkerStyle(21);
        value[i]->SetLineColor(kBlack);
        value[i]->SetLineWidth(2);
        value[i]->GetXaxis()->SetTitle("Num");
        value[i]->GetYaxis()->SetTitle("(Simulated - Analytical");
        value[i]->SetTitle("Accuracies");
        
        diff[i]->SetMarkerColor(kBlue);
        diff[i]->SetMarkerStyle(3);
        diff[i]->SetLineColor(kBlack);
        diff[i]->SetLineWidth(2);
        diff[i]->GetXaxis()->SetTitle("sqrt(Num)");
        diff[i]->GetYaxis()->SetTitle("1/|Simulated - Analytical|");
        diff[i]->SetTitle("Inversed Accuracies");         

	canv[i]->cd(1);
	gPad->SetLogx();
        gPad->SetGrid();
//	value[i]->GetYaxis()->SetRangeUser(0.0, (anal[i]+1.5));//set range, must be here(after deciding the site)---
        value[i]->Draw("APL");

        canv[i]->cd(2);
	gPad->SetGrid();
	diff[i]->Draw("APL");

//	canv[i]->SaveAs("%d-D_results.png", (i+2));
        canv[i]->Write();
    }

    file_1->Close();



}
 
#ifndef __CINT__
int main(){
    work_1();
}
#endif


