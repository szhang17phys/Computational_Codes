//Ashely, Katy, Shuaixiang------
//Nov 5, 2021------

#include <iostream>
#include<iomanip>
#include<cmath>
#include<string>

#include "TROOT.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TCanvas.h"            //(x, y) graphs
#include "TLegend.h" 
#include "TArrow.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TLine.h"

#include "cpROOT.hpp"//very important---
#include "cpRandom.hpp"
#include <cstdlib>
#include <ctime>
#include <chrono>

using namespace std::chrono;

void format_h(TH1F* h, int linecolor){
    h->SetLineWidth(3);
    h->SetLineColor(linecolor);
}

//===Main==================================
void random_final(){
    double lo = -1.0;
    double hi = 1.0;
    double value = 0.0;//set variable---
   
    //Generating a Random Seed value with respect to time---
    Ullong seed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    std::cout << seed << std::endl;//output seed---

    cpRandMT Check(seed);//Messener PRNG Function---

    TH1F *h1 = new TH1F("h1", "Random Number", 100, -1, 1);//The histogram---

    for(int i=0; i<= 10000; i++){//Messener PRNG with Range [-1, 1]---
	value = Check.randDouble()*2 -1;
    	h1->Fill(value);	
    }	

//---Graph------
    h1->SetLineColor(kRed);
    h1->SetLineWidth(2);
    h1->GetXaxis()->SetTitle("Value");
    h1->GetYaxis()->SetTitle("Numbers");
    h1->GetYaxis()->SetRangeUser(0, 150);//set range of y axis---

    TLine *line = new TLine(-1, 100, 1, 100);
    line->SetLineColor(kBlue);
    line->SetLineStyle(2);
    line->SetLineWidth(2);

    TCanvas *c1 = new TCanvas("c1", "c1", 700, 500);
    gStyle->SetOptStat("nemr");

    h1->Draw();
    line->Draw();

    c1->SaveAs("graphs.png");
     
}
 
#ifndef __CINT__
int main(){
    random_final();
}
#endif

