#ifndef COMP_SET2_H
#define COMP_SET2_H

#include <cmath>
#define PI acos(-1)

//===For prob 1; transformed method=========
double f1(double x){
    return exp(-x*x)*2.0/sqrt(PI);
}

double f2(double x){
    return exp(-x);
}

double f3(double x){
    return exp(-x)*cos(x*1.0);
}

double f1_trans(double x){
    double part1 = 0.0;
    double part2 = 0.0;

    part1 = 2.0/sqrt(PI)*exp(-1.0*x*x);
    part2 = 2.0/sqrt(PI)*exp(-1.0/(1.0*x*x))/(1.0*x*x);

    return (part1+part2);
}

double f2_trans(double x){
    double part1 = 0.0;
    double part2 = 0.0;

    part1 = exp(-1.0*x);
    part2 = 1.0*exp(-1.0/x)/(1.0*x*x);

    return (part1+part2);
}

double f3_trans(double x){
    double part1 = 0.0;
    double part2 = 0.0;

    part1 = exp(-1.0*x)*cos(x*1.0);
    part2 = 1.0*exp(-1.0/x)*cos(1.0/x)/(1.0*x*x);

    return (part1+part2);
}

double triple(double (*func)(double), int N ){
    double point_x = 0.0;//the value of point---
    double sum = 0.0;//to collect the sum---
    int pointNum = 0;//num of points---

    pointNum = pow(3, N);    

    for(int i=0; i<pointNum; ++i){
	point_x = (i+0.5)/pointNum;
	sum += func(point_x);    
    }
    sum = 1.0*sum/pointNum;

    return sum;
}

//===For prob 1; Simpson method=========
double simpson(double (*func)(double), double a, double b, int N){
//a and b are peripheral points of the function---
    double point = 0.0;
    double h = 0.0;//h is the length of step---
    double sum_2 = 0.0;
    double sum_4 = 0.0;
    double sum_total = 0.0;
    int pointNum = 0;

    pointNum = pow(2, N);
    h = (b-a)*1.0/pointNum;

    for(int i=1; i<pointNum; ++i){
	point = a+i*h;
        sum_2 += 2*func(point);
    }
    for (int i=1; i<pointNum; i=i+2){
        point = a+i*h;
        sum_4 += 2*func(point);   
    }
    sum_total = (sum_2 + sum_4 + func(a) + func(b))*h/3.0;
    
    return sum_total;
}


//===For prob 2===========================
double delta(double x, double ep){//original approximate---
    return 2.0*ep/(x*x+ep*ep)/PI;
}

double simpson2(double (*func)(double, double), double epsilon, double range, int N){
//epsilon, region: [0, range], N is number of integrals---
    double point = 0.0;
    double h = 0.0;//h is the length of step---
    double sum_2 = 0.0;
    double sum_4 = 0.0;
    double sum_total = 0.0;
    int pointNum = 0;

    pointNum = pow(2, N);//number of integrals---
    h = 1.0*range/pointNum;//width of the bin---

    for(int i=1; i<pointNum; ++i){
	point = i*h;
        sum_2 += 2*func(point, epsilon);
    }
    for (int i=1; i<pointNum; i=i+2){
        point = i*h;
        sum_4 += 2*func(point, epsilon);   
    }
    sum_total = (sum_2 + sum_4 + func(0.0, epsilon) + func(range, epsilon))*h/3.0;
    
    return sum_total;

}

/*
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
*/

#endif
