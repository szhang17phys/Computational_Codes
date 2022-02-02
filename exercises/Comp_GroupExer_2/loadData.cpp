#include <iostream>
#include <cmath>
#include <string>
#include <fstream> //read data from file && output data to file
using namespace std;

int main(){
    //Defining Variables
    string line;
    int count = 0;
    
    ifstream myfile ("data.txt"); //ifstream, input part of fstream
    //This while loop gets the number of lines in the text file
    while(!myfile.eof()){
	(getline(myfile,line));
	count++;
          }

    //Data Arrays
    int x[] = {};	
    int y1[] = {};	
    int y2[] = {};	
    string data;
    if(myfile.is_open()){
        while (getline(myfile, data)){ 
            cout<< data <<'\n';
        }
        myfile.close();
    }


    return 0;
}
