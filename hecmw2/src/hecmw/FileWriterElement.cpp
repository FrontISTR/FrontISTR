/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileWriterElement.cpp
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
#include "FileWriterElement.h"
using namespace FileIO;
CFileWriterElement::CFileWriterElement()
{
    ;
}
CFileWriterElement::~CFileWriterElement()
{
    ;
}
string CFileWriterElement::StrType(const uiint& elem_type)
{
    string sType;
    switch(elem_type) {
    case(pmw::ElementType::Hexa):
        sType= "Hexa";
        break;
    case(pmw::ElementType::Hexa2):
        sType= "Hexa2";
        break;
    case(pmw::ElementType::Tetra):
        sType= "Tetra";
        break;
    case(pmw::ElementType::Tetra2):
        sType= "Tetra2";
        break;
    case(pmw::ElementType::Prism):
        sType= "Prism";
        break;
    case(pmw::ElementType::Prism2):
        sType= "Prism2";
        break;
    case(pmw::ElementType::Quad):
        sType= "Quad";
        break;
    case(pmw::ElementType::Quad2):
        sType= "Quad2";
        break;
    case(pmw::ElementType::Triangle):
        sType= "Triangle";
        break;
    case(pmw::ElementType::Triangle2):
        sType= "Triangle2";
        break;
    case(pmw::ElementType::Beam):
        sType= "Beam";
        break;
    case(pmw::ElementType::Beam2):
        sType= "Beam2";
        break;
    default:
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "CFileWriterElement ", elem_type);
        break;
    }
    return sType;
}
void CFileWriterElement::WriteDebug(ofstream& ofs, const uiint& mgLevel)
{
    pmw::CAssyModel *pAssyModel;
    pmw::CMesh *pMesh;
    pmw::CElement *pElem;
    pmw::CNode   *pNode;
    string white(" ");
    pAssyModel= mpGMGModel->getAssyModel(mgLevel);
    uiint numOfPart= pAssyModel->getNumOfMesh();
    uiint numOfElem, numOfNode;
    uiint i,ii,iii;
    for(i=0; i< numOfPart; i++) {
        pMesh= pAssyModel->getMesh(i);
        ofs << " -- Element Block Start -- " << ", mgLevel == " << mgLevel<< ", Mesh ID==" << pMesh->getMeshID() << endl;
        numOfElem= pMesh->getNumOfElement();
        for(ii=0; ii< numOfElem; ii++) {
            pElem= pMesh->getElementIX(ii);
            ofs << white
                << pElem->getID() << white
                << StrType(pElem->getType()) << white;
            numOfNode= pElem->getNumOfNode();
            for(iii=0; iii< numOfNode; iii++) {
                pNode= pElem->getNode(iii);
                ofs << pNode->getID() << white;
            };
            {
                ofs << ",Edge Node ";
                uiint numOfEdge,numOfFace, iedge, iface;
                numOfEdge= pElem->getEdgeInterNodeSize();
                for(iedge=0; iedge< numOfEdge; iedge++) {
                    if(pElem->getEdgeInterNode(iedge)) ofs << pElem->getEdgeInterNode(iedge)->getID() << white;
                };
                ofs << ",Face Node ";
                numOfFace= pElem->getFaceNodeSize();
                for(iface=0; iface< numOfFace; iface++) {
                    if(pElem->getFaceNode(iface)) ofs << pElem->getFaceNode(iface)->getID() << white;
                };
            }
            ofs << endl;
        };
        ofs << " -- Element Block End -- " << endl;
        pmw::CCommMesh *pCommMesh;
        pmw::CCommElement *pCommElem,*pDCommElem;
        uiint numOfCommMesh= pMesh->getNumOfCommMesh();
        uiint icom;
        for(icom=0; icom< numOfCommMesh; icom++) {
            pCommMesh= pMesh->getCommMesh(icom);
            ofs << " -- CommElement -- " << endl;
            uiint numOfCommElem= pCommMesh->getNumOfCommElement();
            uiint icomelem;
            for(icomelem=0; icomelem< numOfCommElem; icomelem++) {
                pCommElem= pCommMesh->getCommElement(icomelem);
                ofs << "CommElem in ElemID => " << pCommElem->getElement()->getID() << endl;
            }
            ofs << " -- DCommElement -- " << endl;
            uiint numOfDCommElem= pCommMesh->getNumOfDCommElement();
            for(icomelem=0; icomelem< numOfDCommElem; icomelem++) {
                pDCommElem= pCommMesh->getDCommElement(icomelem);
                ofs << "DCommElem in ElemID => " << pDCommElem->getElement()->getID() << endl;
            }
        }
    };
}
