//
//  FileWriterElement.cpp
//
//
//
//                  2009.07.23
//                  2009.07.23
//                  k.Takeda
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

// method
// --
string CFileWriterElement::StrType(const uint& elem_type)
{
    string sType;

    switch(elem_type){
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
void CFileWriterElement::Write(ofstream& ofs, const uint& mgLevel)
{
   pmw::CAssyModel *pAssyModel;
   pmw::CMesh *pMesh;
   pmw::CElement *pElem;
   pmw::CNode   *pNode;
   string white(" ");

   pAssyModel= mpGMGModel->getAssyModel(mgLevel);

   uint numOfPart= pAssyModel->getNumOfMesh();
   uint numOfElem, numOfNode;
   uint i,ii,iii;
   for(i=0; i< numOfPart; i++){
       pMesh= pAssyModel->getMesh(i);

       ofs << " -- Element Block Start -- " << ", mgLevel == " << mgLevel<< ", Mesh ID==" << pMesh->getMeshID() << endl;

       ////debug
       //cout << "FileWriterElement::Write, numOfElem => " << pMesh->getNumOfElement() << endl;
       
       numOfElem= pMesh->getNumOfElement();
       for(ii=0; ii< numOfElem; ii++){
           
           pElem= pMesh->getElementIX(ii);

           //ofs << pMesh->getMeshID() << white
           ofs << white  //Visualでの確認の為にMeshID出力はとりあえず中止
               << pElem->getID() << white
               << StrType(pElem->getType()) << white;

           numOfNode= pElem->getNumOfNode();
           for(iii=0; iii< numOfNode; iii++){
               pNode= pElem->getNode(iii);

               ofs << pNode->getID() << white;
           };
           //debug
           {
               //ofs << ",Edge Node ";
               //uint numOfEdge,numOfFace, iedge, iface;
               //numOfEdge= pElem->getNumOfEdge();
               //for(iedge=0; iedge< numOfEdge; iedge++){
               //    if(pElem->getEdgeInterNode(iedge)) ofs << pElem->getEdgeInterNode(iedge)->getID() << white;
               //};
               //ofs << ",Face Node ";
               //numOfFace= pElem->getNumOfFace();
               //for(iface=0; iface< numOfFace; iface++){
               //    if(pElem->getFaceNode(iface)) ofs << pElem->getFaceNode(iface)->getID() << white;
               //};
           }
           //debug end

           ofs << endl;

       };//Elementループエンド
       ofs << " -- Element Block End -- " << endl;

       pmw::CCommMesh *pCommMesh;
       pmw::CCommElement *pCommElem,*pDCommElem;
       uint numOfCommMesh= pMesh->getNumOfCommMesh();
       uint icom;
       for(icom=0; icom< numOfCommMesh; icom++){
           pCommMesh= pMesh->getCommMesh(icom);/////////////// <- 本来はCommID,修正の必要あり.

           ofs << " -- CommElement -- " << endl;
           uint numOfCommElem= pCommMesh->getNumOfCommElement();
           uint icomelem;
           for(icomelem=0; icomelem< numOfCommElem; icomelem++){
               pCommElem= pCommMesh->getCommElement(icomelem);
               ofs << "CommElem in ElemID => " << pCommElem->getElement()->getID() << endl;
           }
           ofs << " -- DCommElement -- " << endl;
           uint numOfDCommElem= pCommMesh->getNumOfDCommElement();
           for(icomelem=0; icomelem< numOfDCommElem; icomelem++){
               pDCommElem= pCommMesh->getDCommElement(icomelem);
               ofs << "DCommElem in ElemID => " << pDCommElem->getElement()->getID() << endl;
           }
       }//CommMeshループエンド

   };//Meshループエンド
}














