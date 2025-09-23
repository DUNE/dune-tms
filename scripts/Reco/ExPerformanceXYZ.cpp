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
vector<double> Performance(string filename,string Start_or_End,string direction) {
    ifstream file(filename);
    //check to see if file is valid
    if (!file) {
        cout<<"I'm sorry I can't open this file";
        exit(1);
    }
    //interpret the direction and Start/End input as an index in the read-in file to start from
    int S_or_E;
    if (Start_or_End=="Start"){
        S_or_E=0;
    } else if (Start_or_End=="End"){
        S_or_E=7;
    } else{
        cout<<"I'm sorry I didn't get Start/End";
        exit(1);
    }
    int X_or_Y_or_Z;
    if(direction=="X"){
        X_or_Y_or_Z=0;
    } else if (direction=="Y"){
        X_or_Y_or_Z=1;
    } else if (direction=="Z"){
        X_or_Y_or_Z=2;
    } else{
        cout<<"I'm sorry I didn't get direction";
        exit(1);
    }
    //this gives us our final index, for which data to use
    int i=S_or_E+X_or_Y_or_Z;
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
    //we want to get the number of columns in our vector, too
    int nTracks = MyVector.size();
    //we evaluate a performance metric for position from first file 
    double minPosition=0;
    //we loop through each column in our matrix (each track)
    for (int j=0;j<nTracks;j++) {
        //We evaluate separation of end position from the truth (loop over x,y,z)
        double End_Separation=Performance_Data[i][j]-Performance_Data[i+3][j];
        // check if we should add this point to the data or not - we shouldn't if we get -9999.0 in our entries
        if (Performance_Data[i][j]==-9999||Performance_Data[i+3][j]==-9999){
            continue;
        }else if (Performance_Data[2][j]-Performance_Data[5][j]<-500){
            continue;
        }else {
            p.push_back(End_Separation);
            if(direction=="Z"){
                p.push_back(Performance_Data[i][j]);
            }
            // weighted sum of the start and end separations (we care more about end separation, so use EndW>StartW usually)
            minPosition = minPosition+End_Separation;
            continue;
        }  
    }
    p.push_back(minPosition/nTracks);
    //average it
    //our returned value is the final performance metric for this dataset. Better performance has smaller p.
    return p;
}   

//use the writer on our set of datasets from different parameters to produce all files at once
void writer (string Start_or_End, string direction){
    //loop over parameter choice
    for (int i=0;i<9;i++){
        for (int j=0;j<9;j++){
            //call in the data open it as a file before feeding it into Performance function
            string filename = "Histogram"+Start_or_End+direction+"Data"+to_string(i+1)+to_string(j+15)+".txt";
            ofstream MyFile(filename);
            if (MyFile.is_open()){
                for (const auto& element: Performance("NewDataRun"+to_string(i+1)+to_string(j+15)+".txt",Start_or_End,direction)){
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
    writer("Start","X");
    writer("Start","Y");
    writer("Start","Z");
    writer("End","X");
    writer("End","Y");
    writer("End","Z");
    return 0;
}
