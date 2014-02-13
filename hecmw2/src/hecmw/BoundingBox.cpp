/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BoundingBox.cpp
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "Vertex.h"
#include "BoundingBox.h"
using namespace pmw;
CBoundingBox::CBoundingBox()
{
    mMinCoord.resize(3);
    mMaxCoord.resize(3);
    mvOBBCenter.resize(3);
    mvOBBX.resize(3);
    mvOBBY.resize(3);
    mvOBBZ.resize(3);
    mvE.resize(3);
}
CBoundingBox::~CBoundingBox()
{
}
void CBoundingBox::sizingABB(CSkinFace* pFace)
{
    uiint numOfVert= pFace->getNumOfNode();
    uiint ivert;
    CContactNode* pConNode;
    pConNode= pFace->getNode(0);
    mMinCoord[0]= pConNode->getX();
    mMinCoord[1]= pConNode->getY();
    mMinCoord[2]= pConNode->getZ();
    mMaxCoord[0]= pConNode->getX();
    mMaxCoord[1]= pConNode->getY();
    mMaxCoord[2]= pConNode->getZ();
    for(ivert=1; ivert< numOfVert; ivert++) {
        pConNode= pFace->getNode(ivert);
        if(mMinCoord[0] > pConNode->getX()) mMinCoord[0]= pConNode->getX();
        if(mMinCoord[1] > pConNode->getY()) mMinCoord[1]= pConNode->getY();
        if(mMinCoord[2] > pConNode->getZ()) mMinCoord[2]= pConNode->getZ();
        if(mMaxCoord[0] < pConNode->getX()) mMaxCoord[0]= pConNode->getX();
        if(mMaxCoord[1] < pConNode->getY()) mMaxCoord[1]= pConNode->getY();
        if(mMaxCoord[2] < pConNode->getZ()) mMaxCoord[2]= pConNode->getZ();
    };
}
bool CBoundingBox::judgeABB(CSkinFace* pFace)
{
    CContactNode* pConNode;
    uiint numOfConNode= pFace->getNumOfNode();
    uiint icnode;
    bool bCheck(true);
    for(icnode=0; icnode< numOfConNode; icnode++) {
        pConNode= pFace->getNode(icnode);
        if(!judgeABB(pConNode)) {
            bCheck=false;
            break;
        }
    };
    return bCheck;
}
bool CBoundingBox::judgeABB(CContactNode* pConNode)
{
    bool bCheck(false);
    double X,Y,Z;
    X=pConNode->getX();
    Y=pConNode->getY();
    Z=pConNode->getZ();
    if(X >= mMinCoord[0] && X <= mMaxCoord[0]) {
        if(Y >= mMinCoord[1] && Y <= mMaxCoord[1]) {
            if(Z >= mMinCoord[2] && Z <= mMaxCoord[2]) {
                bCheck=true;
            }
        }
    }
    return bCheck;
}
void CBoundingBox::sizingOBB(CSkinFace* pFace)
{
    uiint numOfVert= pFace->getNumOfNode();
    uiint ivert;
    CContactNode* pConNode;
    mvOBBCenter[0]=0.0;
    mvOBBCenter[1]=0.0;
    mvOBBCenter[2]=0.0;
    for(ivert=0; ivert< numOfVert; ivert++) {
        pConNode= pFace->getNode(ivert);
        mvOBBCenter[0]+= pConNode->getX();
        mvOBBCenter[1]+= pConNode->getY();
        mvOBBCenter[2]+= pConNode->getZ();
    };
    mvOBBCenter[0]/= (double)numOfVert;
    mvOBBCenter[1]/= (double)numOfVert;
    mvOBBCenter[2]/= (double)numOfVert;
    CContactNode *pVert0,*pVert1,*pVert2;
    pVert0= pFace->getNode(0);
    pVert1= pFace->getNode(1);
    pVert2= pFace->getNode(2);
    double X1= pVert1->getX() - pVert0->getX();
    double Y1= pVert1->getY() - pVert0->getY();
    double Z1= pVert1->getZ() - pVert0->getZ();
    double X2= pVert2->getX() - pVert0->getX();
    double Y2= pVert2->getY() - pVert0->getY();
    double Z2= pVert2->getZ() - pVert0->getZ();
    mvOBBZ[0]= Y1*Z2 - Z1*Y2;
    mvOBBZ[1]= Z1*X2 - X1*Z2;
    mvOBBZ[2]= X1*Y2 - Y1*X2;
    double absValue= sqrt(mvOBBZ[0]*mvOBBZ[0] + mvOBBZ[1]*mvOBBZ[1] + mvOBBZ[2]*mvOBBZ[2]);
    mvOBBZ[0] /= absValue;
    mvOBBZ[1] /= absValue;
    mvOBBZ[2] /= absValue;
    mvOBBX[0]= (pVert0->getX() + pVert1->getX())*0.5;
    mvOBBX[1]= (pVert0->getY() + pVert1->getY())*0.5;
    mvOBBX[2]= (pVert0->getZ() + pVert1->getZ())*0.5;
    mvOBBX[0] -= mvOBBCenter[0];
    mvOBBX[1] -= mvOBBCenter[1];
    mvOBBX[2] -= mvOBBCenter[2];
    absValue= sqrt(mvOBBX[0]*mvOBBX[0] + mvOBBX[1]*mvOBBX[1] + mvOBBX[2]*mvOBBX[2]);
    mvOBBX[0] /= absValue;
    mvOBBX[1] /= absValue;
    mvOBBX[2] /= absValue;
    X1= mvOBBZ[0];
    Y1= mvOBBZ[1];
    Z1= mvOBBZ[2];
    X2= mvOBBX[0];
    Y2= mvOBBX[1];
    Z2= mvOBBX[2];
    mvOBBY[0]= Y1*Z2 - Z1*Y2;
    mvOBBY[1]= Z1*X2 - X1*Z2;
    mvOBBY[2]= X1*Y2 - Y1*X2;
    absValue= sqrt(mvOBBY[0]*mvOBBY[0] + mvOBBY[1]*mvOBBY[1] + mvOBBY[2]*mvOBBY[2]);
    mvOBBY[0] /= absValue;
    mvOBBY[1] /= absValue;
    mvOBBY[2] /= absValue;
    mvE[0]=0.0;
    mvE[1]=0.0;
    mvE[2]=0.0;
    double obbX,obbY;
    for(ivert=0; ivert< numOfVert; ivert++) {
        pConNode= pFace->getNode(ivert);
        X1= pConNode->getX() - mvOBBCenter[0];
        Y1= pConNode->getY() - mvOBBCenter[1];
        Z1= pConNode->getZ() - mvOBBCenter[2];
        obbX= abs(X1*mvOBBX[0] + Y1*mvOBBX[1] + Z1*mvOBBX[2]);
        if(ivert==0) mvE[0]=obbX;
        if(obbX > mvE[0]) mvE[0]=obbX;
        obbY= abs(X1*mvOBBY[0] + Y1*mvOBBY[1] + Z1*mvOBBY[2]);
        if(ivert==0) mvE[1]=obbY;
        if(obbY > mvE[1]) mvE[1]=obbY;
        mvE[2]=(mvE[0]+mvE[1])*0.01;
    };
    mvE[0] *= 1.10;
    mvE[1] *= 1.10;
    mvE[2] *= 1.10;
}
bool CBoundingBox::judgeOBB(CContactNode* pConNode)
{
    bool bCheck(false);
    double X1= pConNode->getX() - mvOBBCenter[0];
    double Y1= pConNode->getY() - mvOBBCenter[1];
    double Z1= pConNode->getZ() - mvOBBCenter[2];
    double projX,projY,projZ;
    projX= abs(X1*mvOBBX[0] + Y1*mvOBBX[1] + Z1*mvOBBX[2]);
    projY= abs(X1*mvOBBY[0] + Y1*mvOBBY[1] + Z1*mvOBBY[2]);
    projZ= abs(X1*mvOBBZ[0] + Y1*mvOBBZ[1] + Z1*mvOBBZ[2]);
    if(projX <= mvE[0]) {
        if(projY <= mvE[1]) {
            if(projZ <= mvE[2]) bCheck=true;
        }
    }
    return bCheck;
}
