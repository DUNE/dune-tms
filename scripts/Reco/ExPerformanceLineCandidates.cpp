//These will be important
#include <iostream>
#include <fstream>
#include<cmath>
#include <vector>
#include<string>
using namespace std;

// this function converts strings of numbers in csv format into c++ vectors. A string is inputted, and the vector form outputted
vector<double> convert_string_to_vector(string string){
    // create our vector
    vector<double> vec;
    //create an integer which will store the value of the numbers in string
    double num=0;
    //the following are useful for dealing with decimals and negatives
    bool decimal = false;
    int decimal_place=1;
    bool negative = false;
    // sequentially move along the string
    for (int i = 0; i<string.length(); i++) {
        // if we see a comma or first square bracket just ignore it
        if (string[i] == ',')
            continue;
        if (string[i] == '[')
            continue;
        // if we see a negative number, set our boolean negative to be true
        if (string[i] == '-'){
            negative = true;
            continue;
        }
        // if we see a full stop, indicate that we are dealing with decimals
        if (string[i] == '.'){
            decimal = true;
            continue;
        }
        // if we come across a space or last square bracket, concatenate our number the vector and go to next array index
        if (string[i] == ' '){
            if (negative==true){
                num=num*(-1);
            }
            vec.push_back(num);
            //cout<<num<<endl;
            // we reset all relevant quantities too
            num=0;
            decimal = false;
            decimal_place=1;
            negative = false;
            continue;
        }
        // break when we come to the end (right square bracket)
        if (string[i] == ']'){
            if (negative==true){
                num=num*(-1);
            }
            vec.push_back(num);
            //cout<<num<<endl;
            break;
        }
        // when it sees a number:
        else {
            // subtract str[i] by 48 to convert it to int
            // start with ones, if there is another digit, we x10, add next digit, and repeat
            // if we get to decimals, we multiply by 100, add number, then divide by 100 for first decimal, 1000 for next, and so on
            num=num*pow(10,decimal_place)+(string[i]-48);
            if (decimal==true){
                num=num/pow(10,decimal_place);
                // decimal;_place indicates which decimal place we have gotten to
                decimal_place++;
            }
            continue;
        }
    }
    //finally we output the vector we make from our line of string
    return vec;
}

// a function extracting the performance metric from an Output txt file. Input the name of a valid .txt file and if we want to evaluate the Start or End and X,Y,or Z co-ordinate
vector<double> Performance(string filename,string U_or_V_or_X,string Start_or_End,string direction) {
    ifstream file(filename);
    //check to see if file is valid
    if (!file) {
        cout<<"I'm sorry I can't open this file";
        exit(1);
    }
    //interpret the direction and Start/End input as an index in the read-in file to start from
    int U_V_or_X;
    if(U_or_V_or_X=="U"){
        U_V_or_X=0;
    } else if (U_or_V_or_X=="V"){
        U_V_or_X=1;
    } else if (U_or_V_or_X=="X"){
        U_V_or_X=2;
    } else{
        cout<<"I'm sorry I didn't get direction";
        exit(1);
    }
    int S_or_E;
    if (Start_or_End=="Start"){
        S_or_E=0;
    } else if (Start_or_End=="End"){
        S_or_E=1;
    } else{
        cout<<"I'm sorry I didn't get Start/End";
        exit(1);
    }
    int X_or_Y_or_Z;
    int x_Y_Z_fork;
    if(direction=="X"){
        X_or_Y_or_Z=0;
        x_Y_Z_fork=0;
    } else if (direction=="Y"){
        X_or_Y_or_Z=0;
        x_Y_Z_fork=1;
    }else if (direction=="Z"){
        X_or_Y_or_Z=1;
        x_Y_Z_fork=2;
    } else{
        cout<<"I'm sorry I didn't get direction";
        exit(1);
    }
    //this gives us our final index, for which data to use, i for the reco, k for the truth
    int i=14+4*U_V_or_X+2*S_or_E+X_or_Y_or_Z;
    int k = 3+x_Y_Z_fork+7*S_or_E;
    //prepare a string, vector, and matrix
    string MyData;
    vector<double> MyVector;
    vector<vector<double>> Performance_Data;
    vector<double> p;
    //the string is filled by a line from the txt file, one at a time
    while (getline (file,MyData)) {
        //convert the string to a vector
        MyVector= convert_string_to_vector(MyData);
        //add the vector created by each line to a new row in our matrix
        Performance_Data.push_back(MyVector);
    }
    double minPosition=0;
    int nTracks = MyVector.size();
    //loop through each column in our matrix (each track)
    for (int j=0;j<nTracks;j++) {
        //We evaluate reco-truth
        double End_Separation=Performance_Data[i][j]-Performance_Data[k][j];
        if (Performance_Data[i][j]==-9999||Performance_Data[k][j]==-9999){   
            continue;
        }else if (abs(Performance_Data[i][j]+Performance_Data[k][j])<500 && abs(End_Separation>1500)){
            //cout<<End_Separation<<"    "<<Performance_Data[i][j]<<"    "<<Performance_Data[k][j]<<endl;
            continue;
        }else if (Performance_Data[2][j]-Performance_Data[5][j]<-500||Performance_Data[15][j]-Performance_Data[5][j]<-500||Performance_Data[19][j]-Performance_Data[5][j]<-500||Performance_Data[23][j]-Performance_Data[5][j]<-500){
            //cout<<End_Separation<<"    "<<Performance_Data[i][j]<<"    "<<Performance_Data[k][j]<<endl;
            continue;
        } else {
            // weighted sum of the start and end separations (we care more about end separation, so use EndW>StartW usually)
            p.push_back(End_Separation);
            if(direction=="Z"){
                p.push_back(Performance_Data[i][j]);
            }
            minPosition = minPosition+End_Separation;
            continue;
        }
    }
    p.push_back(minPosition/nTracks);
    //average it
    //our returned value is the final performance metric for this dataset. Better performance has smaller p.
    return p;
}   

void writer (string U_or_V_or_X, string Start_or_End, string direction){
    //loop over parameter choice
    for (int i=0;i<9;i++){
        for (int j=0;j<9;j++){
            //call in the data open it as a file before feeding it into Performance function
            string filename = "Histogram"+U_or_V_or_X+Start_or_End+direction+"Data"+to_string(i+1)+to_string(j+15)+".txt";
            ofstream MyFile(filename);
            if (MyFile.is_open()){
                for (const auto& element: Performance("NewDataRun"+to_string(i+1)+to_string(j+15)+".txt",U_or_V_or_X,Start_or_End,direction)){
                MyFile<<element<<endl;
                }
            }
            //close the file and save it
            MyFile.close();
            cout<<"File created as: "<<filename<<endl;
        }
    }
    return;
}

int main(){
    //here we create a .txt file for the data in the Perfromance X,Y and Z metrics.
    writer("U","Start","X");
    writer("U","Start","Z");
    writer("U","End","X");
    writer("U","End","Z");
    writer("V","Start","X");
    writer("V","Start","Z");
    writer("V","End","X");
    writer("V","End","Z");
    writer("X","Start","Y");
    writer("X","Start","Z");
    writer("X","End","Y");
    writer("X","End","Z");
    return 0;
}
