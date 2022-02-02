#ifndef COMP_SET3_H
#define COMP_SET3_H

#include <cmath>
#define PI acos(-1)

using namespace std;

//===For problem 1====================
//---analytical volume------
double V_ana(int d){
    double vol = 0.0;
    double nomi = 0.0;//nominator---
    double deno = 0.0;//denominator---
    
    nomi = pow(PI, (0.5*d));
    switch(d){
        case 2:
	    deno = 1.0;
	    break;
	case 3:
	    deno = 0.75*sqrt(PI);
	    break;
	case 4:
	    deno = 2.0;
	    break;
	case 5:
	    deno = 1.875*sqrt(PI);
	    break;
	case 6:
	    deno = 6.0;
    }
    vol = 1.0*nomi/deno;

    return vol;
}


//===For problem 2====================
double f(double x){
    return sqrt(1 - x*x);
}

double p(double x){
    return 2*(1 - x);
}

double xDis(double m){
    return (1.0 - sqrt(1 - m));
}


#endif

