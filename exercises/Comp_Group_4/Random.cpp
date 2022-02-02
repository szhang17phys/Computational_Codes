#include "cpRandom.hpp"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <chrono>


using namespace std::chrono;
int main() 
{
	//Setting Variable 
	double lo = -1.0;
	double hi = 1.0;
        double number = 0.0;

	//Generating a Random Seed value with respect to time
	Ullong seed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	std::cout << seed << std::endl;

	//Messener PRNG Function
	cpRandMT Check(seed);

	for(int i=0; i<= 10000; i++)
	{

		//Messener PRNG with Range -1, 1
		
	       number = Check.randDouble()*2 -1;
    		std::cout<<number<<std::endl;
		
		//Plotting PRNG Output
	}	

	return 0;



}


