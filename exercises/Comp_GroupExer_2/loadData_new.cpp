#include <iostream>
#include <cmath>
#include <string>
#include <fstream> //read data from file && output data to file
using namespace std;

//szhang, Sep 15, 2021
//read file, put different columns into different arrays---

int main(){
    int column = 0;
    int row = 0;
    int value;
    int total[100][3];
    int x[100];
    int y1[100];
    int y2[100];
    ifstream myfile; //ifstream, input part of fstream
    myfile.open("data.txt");

    while(myfile >> value){//read word by word; >>: transfer
        total[row][column] = value;
        if(column == 2){
            row++;
            column=0;
        }
        else
            column++;
    }

    for(int i=0; i<100; i++){
        x[i] = total[i][0];
        y1[i] = total[i][1];
        y2[i] = total[i][2];
    }

    for(int i=0; i<100; i++){
        cout<<y2[i]<<endl;
    }

/*    for(int i=0; i<100; i++){
        for(int j=0; j<3; j++){
            cout<<total[i][j];
            cout<<"  ";
        }
        cout<<"\n";
    }
*/
    myfile.close();
    return 0;
}
