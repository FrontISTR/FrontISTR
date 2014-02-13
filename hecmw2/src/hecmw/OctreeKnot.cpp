/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/OctreeKnot.cpp
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
#include "ContactNode.h"
#include "OctreeKnot.h"
using namespace pmw;
COctreeKnot::COctreeKnot()
{
    mvMasterNode.clear();
    mvSlaveNode.clear();
    mvMasterNode.resize(0);
    mvSlaveNode.resize(0);
}
COctreeKnot::~COctreeKnot()
{
}
void COctreeKnot::createChildKnot()
{
    COctreeKnot *pKnot;
    mvChildKnot.reserve(8);
    uiint i;
    for(i=0; i< 8; i++) {
        pKnot= new COctreeKnot;
        pKnot->setParentKnot(this);
        mvChildKnot.push_back(pKnot);
    };
    double inX,inY,inZ;
    inX=(maxX + minX)*0.5;
    inY=(maxY + minY)*0.5;
    inZ=(maxZ + minZ)*0.5;
    mvChildKnot[0]->setX(minX,inX);
    mvChildKnot[0]->setY(minY,inY);
    mvChildKnot[0]->setZ(minZ,inZ);
    mvChildKnot[1]->setX(inX,maxX);
    mvChildKnot[1]->setY(minY,inY);
    mvChildKnot[1]->setZ(minZ,inZ);
    mvChildKnot[2]->setX(inX,maxX);
    mvChildKnot[2]->setY(inY,maxY);
    mvChildKnot[2]->setZ(minZ,inZ);
    mvChildKnot[3]->setX(minX,inX);
    mvChildKnot[3]->setY(inY,maxY);
    mvChildKnot[3]->setZ(minZ,inZ);
    mvChildKnot[4]->setX(minX,inX);
    mvChildKnot[4]->setY(minY,inY);
    mvChildKnot[4]->setZ(inZ,maxZ);
    mvChildKnot[5]->setX(inX,maxX);
    mvChildKnot[5]->setY(minY,inY);
    mvChildKnot[5]->setZ(inZ,maxZ);
    mvChildKnot[6]->setX(inX,maxX);
    mvChildKnot[6]->setY(inY,maxY);
    mvChildKnot[6]->setZ(inZ,maxZ);
    mvChildKnot[7]->setX(minX,inX);
    mvChildKnot[7]->setY(inY,maxY);
    mvChildKnot[7]->setZ(inZ,maxZ);
}
void COctreeKnot::distItem()
{
    double inX,inY,inZ;
    inX=(maxX + minX)*0.5;
    inY=(maxY + minY)*0.5;
    inZ=(maxZ + minZ)*0.5;
    bool bX,bY,bZ;
    bool beX,beY,beZ;
    CContactNode* pConNode;
    uiint numOfItem= mvSlaveNode.size();
    uiint i;
    for(i=0; i< numOfItem; i++) {
        pConNode= mvSlaveNode[i];
        bX=false;
        bY=false;
        bZ=false;
        beX=false;
        beY=false;
        beZ=false;
        uiint distID;
        if(pConNode->getX() < inX) bX=true;
        if(pConNode->getY() < inY) bY=true;
        if(pConNode->getZ() < inZ) bZ=true;
        if( bX &&  bY &&  bZ) {
            mvChildKnot[0]->addSlaveNode(pConNode);
            distID=0;
        }
        if(!bX &&  bY &&  bZ) {
            mvChildKnot[1]->addSlaveNode(pConNode);
            distID=1;
        }
        if(!bX && !bY &&  bZ) {
            mvChildKnot[2]->addSlaveNode(pConNode);
            distID=2;
        }
        if( bX && !bY &&  bZ) {
            mvChildKnot[3]->addSlaveNode(pConNode);
            distID=3;
        }
        if( bX &&  bY && !bZ) {
            mvChildKnot[4]->addSlaveNode(pConNode);
            distID=4;
        }
        if(!bX &&  bY && !bZ) {
            mvChildKnot[5]->addSlaveNode(pConNode);
            distID=5;
        }
        if(!bX && !bY && !bZ) {
            mvChildKnot[6]->addSlaveNode(pConNode);
            distID=6;
        }
        if( bX && !bY && !bZ) {
            mvChildKnot[7]->addSlaveNode(pConNode);
            distID=7;
        }
    };
    numOfItem= mvMasterNode.size();
    for(i=0; i< numOfItem; i++) {
        pConNode= mvMasterNode[i];
        bX=false;
        bY=false;
        bZ=false;
        uiint distID;
        if(pConNode->getX() < inX) bX=true;
        if(pConNode->getY() < inY) bY=true;
        if(pConNode->getZ() < inZ) bZ=true;
        if( bX &&  bY &&  bZ) {
            mvChildKnot[0]->addMasterNode(pConNode);
            distID=0;
        }
        if(!bX &&  bY &&  bZ) {
            mvChildKnot[1]->addMasterNode(pConNode);
            distID=1;
        }
        if(!bX && !bY &&  bZ) {
            mvChildKnot[2]->addMasterNode(pConNode);
            distID=2;
        }
        if( bX && !bY &&  bZ) {
            mvChildKnot[3]->addMasterNode(pConNode);
            distID=3;
        }
        if( bX &&  bY && !bZ) {
            mvChildKnot[4]->addMasterNode(pConNode);
            distID=4;
        }
        if(!bX &&  bY && !bZ) {
            mvChildKnot[5]->addMasterNode(pConNode);
            distID=5;
        }
        if(!bX && !bY && !bZ) {
            mvChildKnot[6]->addMasterNode(pConNode);
            distID=6;
        }
        if( bX && !bY && !bZ) {
            mvChildKnot[7]->addMasterNode(pConNode);
            distID=7;
        }
    };
}
void COctreeKnot::setItemProp()
{
    CContactNode *pConNode;
    uiint numOfItem;
    uiint inode;
    numOfItem= mvMasterNode.size();
    for(inode=0; inode < numOfItem; inode++) {
        pConNode= mvMasterNode[inode];
        pConNode->setOctreeID(mLayer, mID);
    };
    numOfItem= mvSlaveNode.size();
    for(inode=0; inode < numOfItem; inode++) {
        pConNode= mvSlaveNode[inode];
        pConNode->setOctreeID(mLayer, mID);
    };
}
