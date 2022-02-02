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
#include "TCanvas.h"//(x, y) graphs---
#include "TLegend.h" 
#include "TArrow.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TFile.h"

#include "cpROOT.hpp"//very important---
#include "cpRandom.hpp"//random---

#include "../../comp_set3.hpp"//my own headfile---
using namespace std;
#define Anal 0.25*PI

void format_h(TH1F* h, int linecolor){
    h->SetLineWidth(3);
    h->SetLineColor(linecolor);
}

//===Main==================================
void work_2e(){
    srand((unsigned)time(NULL));
    cpRandMT random(rand());

    double Num[8];//different points---
    double simu_1[8];//simulated values, importance---
    double devi_1[8];//deviation bet simu and anal, importance---
    double simu_2[8];//simulated values, sample mean---
    double devi_2[8];//deviation, sample mean---

    double temp_1 = 0.0;//temporary sums of f/g, importance---
    double temp_2 = 0.0;//sample mean---
    double value_x = 0.0;//show distributions of x---

    for(int i=0; i<8; ++i){
	Num[i] = 1.0*pow(10, (i+2));
    }

    for(int i=0; i<7; ++i){//different numbers---
	temp_1 = 0.0;
	temp_2 = 0.0;

	if(i<5){//just calculate case of Num=10^7---
	    continue;
	}
	if(i>5){
	    break;
	}

	double startTimeIS = clock();

	for(int j=0; j<Num[i]; ++j){//importance sampling---
	    value_x = xDis(random.randDouble());
	    temp_1 += 1.0*f(value_x)/p(value_x);
        }
	simu_1[i] = temp_1*1.0/Num[i];
	devi_1[i] = fabs(simu_1[i] - Anal);
	
	double endTimeIS = clock();
	double startTimeSM = clock();

        for(int j=0; j<Num[i]; ++j){//sample mean---
	    value_x = random.randDouble();
	    temp_2 += 1.0*f(value_x);
        }
	simu_2[i] = temp_2*1.0/Num[i];
	devi_2[i] = fabs(simu_2[i] - Anal);

	double endTimeSM = clock();

	cout<<"Num: "<<Num[i]<<"------------------"<<endl;
	cout<<"Deviation 1 (IS)  : "<<devi_1[i]<<endl;    
	cout<<"Time for IS Method: "<<(endTimeIS - startTimeIS)<<"us"<<endl;
	cout<<"Deviation 2 (SM)  : "<<devi_2[i]<<endl;        
	cout<<"Time for SM Method: "<<(endTimeSM - startTimeSM)<<"us"<<endl;

/*	if(devi_1[i]<0.0001){
	    cout<<"Num for IS < 10^{-4}: "<<Num[i]<<" ======"<<endl;
	}
	if(devi_2[i]<0.0001){
	    cout<<"Num for SM < 10^{-4}: "<<Num[i]<<" ======"<<endl;
	}
*/



    }


}
 
#ifndef __CINT__
int main(){
    work_2e();
}
#endif


