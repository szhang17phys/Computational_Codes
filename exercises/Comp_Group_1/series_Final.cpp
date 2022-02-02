#include<iostream>
#include<iomanip>//Standard output
#include<cmath>
using namespace std;

double Sdownf(int N){
    double s=0.0;
    for (int n=N; n>=1; n--)
	s = s + (1.0/n);

    return s; 
}

double Supf(int N){
    double s=0.0;
    for (int n=1; n<=N; ++n)
	s = s + (1.0/n);

    return s;
}

double Subtract(double up, double down){
    double sub=0.0;
    sub=log10(1.0*(up-down)/down);

    return sub;
}



int main (){
    int N=1.0E+08;
    double sumUp=0.0;
    double sumDown=0.0;

    sumUp = Supf(N);
    sumDown = Sdownf(N);

//Output==========
    cout << "--------------------------------------" << endl;
    cout << "For N = " << N << ":" << endl;
    cout << "Sup   = " << fixed << setprecision(20) << sumUp << endl;
    cout << "Sdown = " << fixed << setprecision(20) << sumDown << endl;
    cout << "Subtraction Result: " << endl;
    cout << Subtract(sumUp, sumDown)<<endl;	
    cout << "--------------------------------------" <<endl;

    return 0;
}
