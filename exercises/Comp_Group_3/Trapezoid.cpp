// Group Problem 3
// Brooks, Hagen, Zhang
// 10/11/2021

#include <cmath>
#include <iostream>
#include <iomanip>
#include "cpFunctions.hpp"
#include "cpQuadrature.hpp"

using namespace std;

double Func(double x)
{	
	double PI = 3.14;
	return x*((x*x)-1); // First Integral
//	return 2*PI*cos(2*PI*x); //Second Integral
}
double Trapezoid( double limA, double limB, double nSteps ) {
	double sum = 0; // start sum at 0
	double trap = 0;
	double h = abs(limB-limA)/nSteps; // calculate step size

	//Evaluating function at the two endpoints
	sum =h*0.5*( Func(limA) + Func(limB));

	//Now we need to calculate the middle points
	//Document 05_Integration
	//Equation 5.3
	for (double i = 1; i< nSteps-1; i++) { // Number of interations
		h = h/2
		// Calculate the middle step points
		double x0 = limA + h;
		double new_point = 0;
		for(double j=0;;j++){ 
			trap = h*(Func(limA+(i*h)));

		// Update sum to include the middle points for total trap value
		sum  += trap;
		std::cout<< "Current Trap Value: "<<trap<<std::endl;
		std::cout<< "Current Sum Value: "<<h*sum<<std::endl;
		return sum;
	};

	//return sum;
}

int main()
{
	
	double TrapValue = Trapezoid(0,1,6.2);
	std::cout<< "Integral Value: "<<TrapValue<<std::endl;

}


