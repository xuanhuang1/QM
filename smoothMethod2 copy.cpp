#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.1415926

void smooth2StarQ(std::vector<vertex> &v, std::vector<face> &f, double ar){
    // for each face 
    //  find largest angle && divide into two triangles
    //  if angle > 120, do
    // if the another triangle largest angle > 120 do, if smallest angel < 30 do

}


void smooth2StarQ2(std::vector<vertex> &v, std::vector<face> &f, double ar){
    for (int i=0; i<v.size(); i++) {// for each vertex in mesh
        double xnew = 0;
        double ynew = 0;
        vertex nextV(0,0,0), thisV(0,0,0);

        if(v[i].onBound == 0){  
            for(int j= 0; j<v[i].numOfNeighborFace; j++){ 
                int faceSize = v[i].neighbors.size()/v[i].numOfNeighborFace + 1;
                //cout << faceSize <<endl;
                // for each neihgbor face
                // get the next 2 vertices and forms a triangle with v[i]
                thisV = v[v[i].neighbors[j*faceSize]];
                nextV = v[v[i].neighbors[(j*faceSize+1)]];


                double vecx,vecy, angle = 0.01;
                double lastEdge, nextEdge,thisEdge, tempAR, minAng, thisAng;
                int self, next, ver2;
                double centerX, centerY;


                //for triangle only!!!
                centerX = GetCircumCenterX(thisV.x, thisV.y,nextV.x, nextV.y,v[i].x, v[i].y);
                centerY = GetCircumCenterY(thisV.x, thisV.y,nextV.x, nextV.y,v[i].x, v[i].y);

                std::vector<int> tempInts;
                tempInts.push_back(i);
                tempInts.push_back(v[i].neighbors[j*faceSize]);
                tempInts.push_back(v[i].neighbors[(j*faceSize+1)]);
                cout << i <<" "<<v[i].neighbors[j] <<" "<< v[i].neighbors[(j+1)] <<" "<< endl;

                // in that triangle
                face aFaceInStar = face(tempInts);
                //cout <<"toS "<<aFaceInStar.listToS()<<endl;;
                for(int k=0; k<aFaceInStar.listOfV.size(); k++){ 
                    self = aFaceInStar.listOfV[k];
                    next = aFaceInStar.listOfV[(k+1)%aFaceInStar.listOfV.size()];
                    int last = aFaceInStar.listOfV[(aFaceInStar.listOfV.size()+k-1)%aFaceInStar.listOfV.size()];
                    //cout << self <<" "<< endl;

                    vecx = v[next].x - v[self].x;   
                    vecy = v[next].y - v[self].y; 
                    double vec1x = v[last].x - v[self].x;
                    double vec1y = v[last].y - v[self].y;

                    thisEdge = sqrt(pow(vecx,2)+ pow(vecy,2));
                    thisAng = acos((vec1x*vecx + vec1y*vecy)
                             /(sqrt(pow(vec1x,2)+pow(vec1y,2))
                               *sqrt(pow(vecx,2)+pow(vecy,2)) )  );

                    // find max and min edge
                    if(k==0){
                        minAng = thisAng;
                        lastEdge = sqrt(pow(vecx,2)+ pow(vecy,2));
                        nextEdge = sqrt(pow(vec1x,2)+ pow(vec1y,2));

                        if(lastEdge < nextEdge){
                            angle = -0.01;
                            ver2 = last;
                        }else{
                            angle = 0.01;
                            ver2 = next;
                        }

                    }else{
                        if(thisAng < minAng)
                            minAng = thisAng;
                    }
                }

                double xTobeMove = v[i].x, yTobeMove = v[i].y;
                //cout << xTobeMove <<" " <<yTobeMove <<endl;
                // 2 for fine 1.5 for good mesh

                if ( max(lastEdge, nextEdge)/min(lastEdge, nextEdge) < ar )
                {
                    //cout <<"ratoi" <<endl;
                    //cout << max(lastEdge, nextEdge)/min(lastEdge, nextEdge) <<endl;
                }else{

                //if(minAng < 30*PI/180){

                    //cout << "min " <<minAng*PI/180 <<endl;
                    //cout << minAng << " "<<endl;
                    double count = 3;

                    while(count > 0){
                        xTobeMove = movePX(xTobeMove, yTobeMove, centerX, centerY, angle, 1);
                        yTobeMove = movePY(xTobeMove, yTobeMove, centerX, centerY, angle, 1);
                        count --;
                //cout << xTobeMove <<" " <<yTobeMove <<endl;

                        //cout << "thisEdge+0.01*max: " <<thisEdge+0.01*max << " currentDist: " << currentDist<<endl;
                    }
                    
                }


                //}
                double t1, t2;
                double ordis = findShortestDistInStar(v, f, v[i].x, v[i].y, v[i].neighbors, t1);
                double dis = findShortestDistInStar(v, f, xTobeMove, yTobeMove, v[i].neighbors, t2);
                /*if(dis < ordis){
                    xTobeMove = v[i].x;
                    yTobeMove = v[i].y;
                }*/


                xnew += xTobeMove;
                ynew += yTobeMove;


            }

            v[i].x = xnew/v[i].numOfNeighborFace;
            v[i].y = ynew/v[i].numOfNeighborFace;

        }
    }


}


/* backup
void smooth2Star(std::vector<vertex> &v, std::vector<face> &f, double ar){
    for (int i=0; i<v.size(); i++) {// for each vertex in mesh
        double xnew = 0;
        double ynew = 0;
        vertex nextV(0,0,0), thisV(0,0,0);

        if(v[i].onBound == 0){  
            for(int j= 0; j<v[i].numOfNeighborFace; j++){ // for each vertex in neihgbor face
                // get the next vertex and forms a triangle with v[i]
                thisV = v[v[i].neighbors[j]];
                nextV = v[v[i].neighbors[(j+1)%v[i].numOfNeighborFace]];


                double vecx,vecy, angle = 0.01;
                double lastEdge, nextEdge,thisEdge, tempAR, minAng, thisAng;
                int self, next, ver2;
                double centerX, centerY;


                //for triangle only!!!
                centerX = GetCircumCenterX(thisV.x, thisV.y,nextV.x, nextV.y,v[i].x, v[i].y);
                centerY = GetCircumCenterY(thisV.x, thisV.y,nextV.x, nextV.y,v[i].x, v[i].y);

                std::vector<int> tempInts;
                tempInts.push_back(i);
                tempInts.push_back(v[i].neighbors[j]);
                tempInts.push_back(v[i].neighbors[(j+1)%v[i].numOfNeighborFace]);

                // in that triangle
                face aFaceInStar = face(tempInts);
                for(int k=0; k<aFaceInStar.listOfV.size(); k++){ 
                    self = aFaceInStar.listOfV[k];
                    next = aFaceInStar.listOfV[(k+1)%aFaceInStar.listOfV.size()];
                    int last = aFaceInStar.listOfV[(aFaceInStar.listOfV.size()+k-1)%aFaceInStar.listOfV.size()];
                    //cout << self <<" "<< next <<" "<< endl;

                    vecx = v[next].x - v[self].x;   
                    vecy = v[next].y - v[self].y; 
                    double vec1x = v[last].x - v[self].x;
                    double vec1y = v[last].y - v[self].y;

                    thisEdge = sqrt(pow(vecx,2)+ pow(vecy,2));
                    thisAng = acos((vec1x*vecx + vec1y*vecy)
                             /(sqrt(pow(vec1x,2)+pow(vec1y,2))
                               *sqrt(pow(vecx,2)+pow(vecy,2)) )  );

                    // find max and min edge
                    if(k==0){
                        minAng = thisAng;
                        lastEdge = sqrt(pow(vecx,2)+ pow(vecy,2));
                        nextEdge = sqrt(pow(vec1x,2)+ pow(vec1y,2));

                        if(lastEdge < nextEdge){
                            angle = -0.01;
                            ver2 = last;
                        }else{
                            angle = 0.01;
                            ver2 = next;
                        }

                    }else{
                        if(thisAng < minAng)
                            minAng = thisAng;
                    }
                }

                double xTobeMove = v[i].x, yTobeMove = v[i].y;

                // 2 for fine 1.5 for good mesh

                if ( max(lastEdge, nextEdge)/min(lastEdge, nextEdge) < ar )
                {
                    //cout <<"ratoi" <<endl;
                    //cout << max(lastEdge, nextEdge)/min(lastEdge, nextEdge) <<endl;
                }else{


                //cout << max << " "<< min <<endl;

                //if(minAng < 30*PI/180){

                    //cout << "min " <<minAng*PI/180 <<endl;
                    //cout << minAng << " "<<endl;
                    double count = 3;

                    while(count > 0){
                        xTobeMove = movePX(xTobeMove, yTobeMove, centerX, centerY, angle);
                        yTobeMove = movePY(xTobeMove, yTobeMove, centerX, centerY, angle);
                        count --;
                        //cout << "thisEdge+0.01*max: " <<thisEdge+0.01*max << " currentDist: " << currentDist<<endl;
                    }
                    
                }


                //}
                double ordis = findShortestDistInStar(v, f, v[i].x, v[i].y, v[i].neighbors);
                double dis = findShortestDistInStar(v, f, xTobeMove, yTobeMove, v[i].neighbors);
                if(dis < ordis){
                    xTobeMove = v[i].x;
                    yTobeMove = v[i].y;
                }


                xnew += xTobeMove;
                ynew += yTobeMove;


            }

            v[i].x = xnew/v[i].numOfNeighborFace;
            v[i].y = ynew/v[i].numOfNeighborFace;

        }
    }


}*/
