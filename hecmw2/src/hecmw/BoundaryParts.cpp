/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BoundaryParts.cpp
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
#include "BoundaryParts.h"
using namespace pmw;
CBoundaryParts::CBoundaryParts()
{
    ;
}
CBoundaryParts::~CBoundaryParts()
{
    ;
}
void CBoundaryParts::resizeBNode(const uiint& res_size)
{
    mvBNode.resize(res_size);
}
void CBoundaryParts::setBNode(const uiint& ivert, CBoundaryNode* pBNode)
{
    mvBNode[ivert]= pBNode;
    mmBNodeID2Index[pBNode->getID()] = ivert;
}
void CBoundaryParts::setupVertexElemID()
{
    CBoundaryNode* pBNode;
    uiint numOfBNode= mvBNode.size();
    uiint ibnode;
    for(ibnode=0; ibnode < numOfBNode; ibnode++) {
        pBNode= mvBNode[ibnode];
        pBNode->setAggElemID(mnID);
    };
}


double CBoundaryParts::getCenterX()
{
    double cX(0.0);
    for(uiint i=0; i < mvBNode.size(); i++) {
        CNode *pNode = mvBNode[i]->getNode();
        cX += pNode->getX();
    }
    cX /= mvBNode.size();

    return cX;
}
double CBoundaryParts::getCenterY()
{
    double cY(0.0);
    for(uiint i=0; i < mvBNode.size(); i++) {
        CNode *pNode = mvBNode[i]->getNode();
        cY += pNode->getY();
    }
    cY /= mvBNode.size();

    return cY;
}
double CBoundaryParts::getCenterZ()
{
    double cZ(0.0);
    for(uiint i=0; i < mvBNode.size(); i++) {
        CNode *pNode = mvBNode[i]->getNode();
        cZ += pNode->getZ();
    }
    cZ /= mvBNode.size();

    return cZ;
}
vdouble CBoundaryParts::getCenterCoord()
{
    vdouble vCoord;
    vCoord.resize(3);

    vCoord[0]= getCenterX();
    vCoord[1]= getCenterY();
    vCoord[2]= getCenterZ();

    return vCoord;
}
//--
// ディレクレ境界に数式処理をする場合の関数
//--
double CBoundaryParts::getCalcValue(CBoundaryNode* pBNode, const double& ent_val, const uiint& dof, vuint& vPolandDOF, map<uiint,CPoland*>& mPoland)
{
    bool bExistFlg(false);

    if(!mPoland.empty()) {
        uiint low = 0, high= vPolandDOF.size()-1, ix;
        uiint limit= vPolandDOF.size()-1;
        ////cout << "BoundaryParts::getCalcValue --- A" << endl;

        while(low <= high) {
            ix= (low+high)/2;
            if(dof==vPolandDOF[ix]) {
                bExistFlg=true;
                break;
            } else if(dof < vPolandDOF[ix]) {
                if(ix==0) break;
                high= ix-1;
            } else {
                if(ix >= limit) break;
                low = ix+1;
            }
        };
        ////cout << "BoundaryParts::getCalcValue --- B    ix:" << ix << endl;
    }

    if(bExistFlg) {
        double x=pBNode->getX(),y=pBNode->getY(),z=pBNode->getZ();
        CCalc* pCalc=CCalc::Instance();

        pCalc->setElementParam( ent_val, x,y,z );
        double calcVal= pCalc->Exec(mPoland[dof]);

        return calcVal;
    } else {
        return ent_val;
    }
}



