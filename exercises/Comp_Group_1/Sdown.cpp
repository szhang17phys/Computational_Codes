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
