#include<iostream>
using namespace std;

int main() {

	int number = 1;
	double summed_number = 0.0;
	double calculation = 0.0;
	
	while( number < 10){
		cout<<"Number \n";
		cout << number<< "\n";
		calculation = 1.0/number;
		summed_number += calculation;
		cout<<"Summed Number \n";
		cout << summed_number<< "\n";
		number =number+1;
		}	
return 0;

}

#include <iostream>
#include <iomanip>
using namespace std;

float Sdownf (int N)
{
	float s=0;
	for (int n=N; n>=1; n--)
		s = s + (1.0/n);
	return s; 
}

double Sdownd (int N)
{
	double s=0;
	for (int n=N; n>=1; n--)
		s = s + (1.0/n);
	return s;
}

int main ()
{
	float sumf;
	double sumd;
	int N=1000;
	sumf = Sdownf (N);
	sumd = Sdownd (N);
	cout << "Sdown float = " << sumf << '\n';
	cout << "Sdown double = " << sumd << '\n';
	return 0;
}
