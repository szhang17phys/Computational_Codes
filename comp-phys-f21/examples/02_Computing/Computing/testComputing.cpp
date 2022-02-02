#include <iostream>
// #include <fstream>
#include <sstream>
#include <string>

using namespace std;

int main( int argc, char* argv[] )
{
  cout << "Hello World" << endl;     // write to the terminal
  cout << "Enter a number: ";

  double x;                          // user input from the terminal
  cin >> x;
  cout << "x = " << x << endl;

  // double xx;                         // read from a file
  // ifstream ifl( "myfile.dat" );
  // ifl >> xx;

  string s( "Hello World " );        // writing to a string
  int i = 2;
  ostringstream oss;
  oss << s << i;
  cout << oss.str() << endl;         // output: Hello World 2 

  return 0;
}
