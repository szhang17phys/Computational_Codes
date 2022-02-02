#ifndef COMP_SET1_H
#define COMP_SET1_H

#include <cmath>
#define PI acos(-1)

//For problem 1------
double Energy(double s, double dsdt ){
    double ene = 0.0;
    double T = 0.0;
    double V = 0.0;
    
    T = 1.0/(32*pow(PI, 4))*pow(dsdt, 2); 
    V = 1.0/(4*pow(PI, 2))*(1-cos(s));
    ene = T+V;
    ene = ene*10000;
    return ene;
}

//For problem 1, about pendulum------
class pend{
    public:
        double s; //theta---
        double dsdt; //d theta / dt---
        double ddsddt; //d^2 theta / dt^2---
        double energy(void);//energy---

};

//For problem 1, energy, element function---
double pend::energy(void){
    
    return Energy(s, dsdt);
}

//For problem 1, motion------
double f(double theta){
    double result = 0.0;
    result = -4*PI*PI*sin(theta);

    return result;
}

//For proble 1, small angle approximation------
double theta(double Amp, double t){
    double angle = 0.0;
    double omega = 2*PI;
    
    angle = Amp*cos(omega*t);

    return angle;
}

//For preblem 1, calculating devation of the array------
double Dev(const double arr[], int size){
    double mean = 0.0;
    double dev = 0.0;

    for(int i=0; i<size; ++i){
        mean += arr[i];
    }
    mean = 1.0*mean/size;

    for(int i=0; i<size; ++i){
        dev += pow((arr[i] - mean), 2);
    }
    dev = 1.0*dev/size;    

    return dev;
}

//For problem 2------
const double g = 9.8;
const double C2_m = 0.014;

//For problem 2, calculate dvx/dt------
double func_x(double vx, double vy){
    double v2 = 0.0;
    double ax = 0.0;

    v2 = sqrt(vx*vx + vy*vy);
    ax = - C2_m*v2*vx;

    return ax;
}

//For problem 2, calculate dvy/dt------
double func_y(double vx, double vy){
    double v2 = 0.0;
    double ay = 0.0;

    v2 = sqrt(vx*vx + vy*vy);
    ay = - (g + C2_m*v2*vy);

    return ay;
}

//For problem 2, calclate the range------
double Range(double vx, double vy){
    double range = 0.0;

    range = 2.0*vx*vy/g;

    return range;
}

//For problem 2, about ballistic motion------
class ball{
    public:
        double x; //x---
        double y;//y---
	double dxdt;
	double dydt;
	double ddxddt;
	double ddyddt;
        double range(void);//range---

};


#endif
