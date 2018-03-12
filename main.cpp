//
//  main.cpp
//
//
//  Created by Xuan Huang on 6/23/16.
//
//

#include "functions.h"
#include <string.h>
#include <math.h>

#define PI 3.1415926
#define ARLARGE 8
#define ARPRT 5
#define ANGPRT 9

using namespace std;

void printStats(vector<vertex> &v, vector<face> &f){
    double arStats[ARPRT];
    double angMaxStats[ANGPRT];
    double angMinStats[ANGPRT];
    double averAR = 0, averMin =0, averMax = 0;

    double angleInter = 90.0/ANGPRT;

    for (int i = 0; i < ARPRT; ++i)
        arStats[i] = 0;

    for (int i = 0; i < ANGPRT; ++i){
        angMaxStats[i] = 0;
        angMinStats[i] = 0;
    }

    for (int i = 0; i < f.size(); ++i){
        double tempAR = f[i].aspectR;
        double tempMin = f[i].minAng*180/PI;
        double tempMax = f[i].maxAng*180/PI;
        averAR += tempAR/f.size();
        averMin += tempMin/f.size();
        averMax += tempMax/f.size();

        //cout << "faces: " << tempAR <<" "<<tempMin <<" "<<tempMax <<endl;

        for (double j = ARPRT-1; j > -1; j--){
            if(tempAR > ARLARGE*j/(ARPRT-1)){
                arStats[(int)j]++;
                break;
            }
        }

        for (int j = 0; j < ANGPRT; ++j){
            if((tempMin > j*angleInter) && (!(tempMin > (j+1)*angleInter)) )
                angMinStats[j]++;
            if((tempMax > j*angleInter+90) && (!(tempMax > (j+1)*angleInter+90)))
                angMaxStats[j]++;
        }
        
    }

    cout <<endl<< "aver: AR " << averAR <<" min "<<averMin <<" max "<<averMax<<endl<<endl;

    for (int i = 1; i < (ARPRT); ++i){
        cout << "ar below " << i*ARLARGE*0.25 << ": " << arStats[i-1] 
        <<"   \t" <<int(arStats[i-1]*100/f.size())<<"%\t";
        for (int j = 0; j < 20; ++j)
            if(arStats[i-1]*100/f.size()> j*5) cout <<"*";
        cout<<endl;
    }
    int i = ARPRT;
    cout << "ar above " << (i-1)*ARLARGE*0.25 << ": " << arStats[i-1] 
        <<"   \t" <<int(arStats[i-1]*100/f.size())<<"%\t";
    for (int j = 0; j < 20; ++j)
        if(arStats[i-1]*100/f.size()> j*5) cout <<"*";
    cout<<endl;

    cout<<endl;
    for (int i = 1; i < (ANGPRT+1); ++i){
        cout << "minang " << (i-1)*angleInter <<"-" <<(i)*angleInter << ": " << angMinStats[i-1] 
        <<"  \t" <<int(angMinStats[i-1]*100/f.size())<<"%\t";
        for (int j = 0; j < 20; ++j)
            if(angMinStats[i-1]*100/f.size()> j*5) cout <<"*";
        cout<<endl;
    }
    cout<<endl;
    for (int i = 1; i < (ANGPRT+1); ++i){
        cout << "maxang " << (i-1)*angleInter+90 <<"-" <<(i)*angleInter+90 << ": " << angMaxStats[i-1] 
        <<"  \t" <<int(angMaxStats[i-1]*100/f.size())<<"%\t";
        for (int j = 0; j < 20; ++j)
            if(angMaxStats[i-1]*100/f.size()> j*5) cout <<"*";
        cout<<endl;
    }
}



int main(int argc, char* argv[]){
    if(argc != 5){
        cout << "Usage: ./test [flag: -Lap -s -q -qStar -qStar2] inputOff outputOff %%itr" <<endl;
        return 1; 
    }

    string runFlag = argv[1];
    string infi = argv[2];
    string outfi = argv[3];
    int itr = stoi(argv[4]);

    ifstream inputFile(infi);
    
    if(!inputFile){
        cout << "Cannot open file" << endl;
        return 1;
    }

    vector<vertex> v;
    vector<edge> e;
    vector<face> f;
    
    readIn(v, e, f, infi);

    double maxAngle, minAngle, med, aspectratio;
    maxminAng(v, f, maxAngle, minAngle);
    aspectratio = aspectR(v, f, med);

    cout << "\noriginal :" <<endl;
    cout << "max: " << maxAngle*180/PI << " min: " << minAngle*180/PI<<endl;
    cout << "aspectR: " << aspectratio<<endl;

    printStats(v,f);
    
    int a =1;
    float thresh = 3;
    while(a >0 ){
        if(!runFlag.compare("-Lap")){
            smoothLapAng(v, f);
        }else if(!runFlag.compare("-s")){
            for (int i = 0; i < 8; ++i){
                if(aspectratio < 4){
                    thresh = 1.2;
                }else if(aspectratio < 8 ){
                    thresh = aspectratio /4;
                }
                smooth2Star(v, f, thresh);
            }
        }else if(!runFlag.compare("-q")){
            printf("\n------The Smoothing-----\n\n" );
            for (int i = 0; i < itr; ++i){
                /*if(aspectratio < 4){
                    thresh = 1.2;
                }else if(aspectratio < 8 ){
                    thresh = aspectratio /4;
                }*/
                smooth2Q(v, f, thresh);
                cout <<"*";
            }
            cout <<"\n";
            printf("\n------Ends-----\n" );

        }else if(!runFlag.compare("-qStar")){
            printf("\n------The Smoothing-----\n\n" );
            for (int i = 0; i < itr; ++i){
                /*if(aspectratio < 4){
                    thresh = 1.2;
                }else if(aspectratio < 8 ){
                    thresh = aspectratio /4;
                }*/
                smooth2QStar(v, f, thresh);
                cout <<"*";
            }
            cout <<"\n";
            printf("\n------Ends-----\n" );

        }else if(!runFlag.compare("-qStar2")){
            printf("\n------The Smoothing-----\n\n" );
            for (int i = 0; i < itr; ++i){
                /*if(aspectratio < 4){
                    thresh = 1.2;
                }else if(aspectratio < 8 ){
                    thresh = aspectratio /4;
                }*/
                smooth2QStar2(v, f, thresh);
                cout <<"*";
            }
            cout <<"\n";
            printf("\n------Ends-----\n" );

        }else{
            cout << "Wrong flag!" <<endl;
            return 1;
        }
        a--;
    }


    cout << "\nafter smooth : "<<endl;
    aspectratio = aspectR(v, f, med);
    maxminAng(v, f, maxAngle, minAngle);
    cout << "max: " << maxAngle*180/PI << " min: " << minAngle*180/PI<<endl;
    cout << "aspectR: " << aspectratio<<endl;
    /*
    for (int i=0; i<v.size(); i++) {
        cout << v[i].x << " "<< v[i].y << " "<< v[i].z << " "<< endl;
    }
    for (int i=0; i<e.size(); i++) {
        cout << e[i].node1 << " "<< e[i].node2<< endl;
    }
    cout << f.size() <<endl;*/

    printStats(v,f);

    ofstream outputfile;
    outputfile.open(outfi, ios::out | ios::trunc);
    
    outputfile << "OFF" << endl;
    outputfile << v.size()<<" "<< f.size()<<" "<<e.size()<<" "<< endl;
    
    for (int i =0; i<v.size(); i++){
        outputfile << v[i].x << " " << v[i].y << " " << v[i].z <<endl;
    }
    
    for (int i =0; i<f.size(); i++){
        outputfile << f[i].listOfV.size() << " "<< f[i].listToS();
        outputfile <<"\n";
    }
    outputfile.close();


    //double s = lineDistPoint(1, 1, 2, 2, 3, 4);
    //cout << s <<endl;
    
    return 0;
    
}