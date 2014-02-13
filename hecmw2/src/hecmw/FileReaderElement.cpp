/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderElement.cpp
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
#include "ElementType.h"
#include "FileReaderElement.h"
using namespace FileIO;
CFileReaderElement::CFileReaderElement()
{
}
CFileReaderElement::~CFileReaderElement()
{
}
string CFileReaderElement::Name()
{
    return  "FileReaderElement";
}

bool CFileReaderElement::Read(ifstream& ifs, string& sLine)
{
    string sTypeName;
    uiint   nElementID, numOfElem, nMeshID, maxID, minID;
    vint   vLocalNode;
    uiint   mgLevel(0);
    uiint  i;
    if(TagCheck(sLine, FileBlockName::StartElement()) ) {
        sLine = getLineSt(ifs);
        istringstream iss(sLine.c_str());
        iss >> numOfElem >> nMeshID >> maxID >> minID;
        mpFactory->reserveElement(mgLevel, nMeshID, numOfElem);
        mpFactory->initBucketElement(mgLevel, nMeshID, maxID, minID);
        uiint nCount(0);
        while(!ifs.eof()) {
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndElement()) ) break;
            istringstream iss(sLine.c_str());
            iss >> sTypeName;
            vLocalNode.clear();
            if(sTypeName=="Tetra") {
                vLocalNode.resize(4);

                iss >> nElementID;
                for(i=0; i< vLocalNode.size(); i++) {
                    iss >> vLocalNode[i];
                }

                mpFactory->GeneElement(mgLevel, nMeshID, nElementID, pmw::ElementType::Tetra, vLocalNode);
                mpFactory->setIDBucketElement(mgLevel, nMeshID, nElementID, nCount);
                nCount++;
            }
            if(sTypeName=="Hexa") {
                vLocalNode.resize(8);
                iss >> nElementID;
                for(i=0; i< vLocalNode.size(); i++) {
                    iss >> vLocalNode[i];
                }
                mpFactory->GeneElement(mgLevel, nMeshID, nElementID, pmw::ElementType::Hexa,vLocalNode);
                mpFactory->setIDBucketElement(mgLevel, nMeshID, nElementID, nCount);
                nCount++;
            }
            if(sTypeName=="Prism") {
                vLocalNode.resize(6);
                iss >> nElementID;
                for(i=0; i< vLocalNode.size(); i++) {
                    iss >> vLocalNode[i];
                }
                mpFactory->GeneElement(mgLevel, nMeshID, nElementID, pmw::ElementType::Prism,vLocalNode);
                mpFactory->setIDBucketElement(mgLevel, nMeshID, nElementID, nCount);
                nCount++;
            }
            if(sTypeName=="Triangle") {
                vLocalNode.resize(3);
                iss >> nElementID;
                for(i=0; i< vLocalNode.size(); i++) {
                    iss >> vLocalNode[i];
                }
                mpFactory->GeneElement(mgLevel, nMeshID, nElementID, pmw::ElementType::Triangle,vLocalNode);
                mpFactory->setIDBucketElement(mgLevel, nMeshID, nElementID, nCount);
                nCount++;
            }
            if(sTypeName=="Quad") {
                vLocalNode.resize(4);
                iss >> nElementID;
                for(i=0; i< vLocalNode.size(); i++) {
                    iss >> vLocalNode[i];
                }
                mpFactory->GeneElement(mgLevel, nMeshID, nElementID, pmw::ElementType::Quad,vLocalNode);
                mpFactory->setIDBucketElement(mgLevel, nMeshID, nElementID, nCount);
                nCount++;
            }
            if(sTypeName=="Beam") {
                vLocalNode.resize(2);
                iss >> nElementID;
                for(i=0; i< vLocalNode.size(); i++) {
                    iss >> vLocalNode[i];
                }
                mpFactory->GeneElement(mgLevel, nMeshID, nElementID, pmw::ElementType::Beam,vLocalNode);
                mpFactory->setIDBucketElement(mgLevel, nMeshID, nElementID, nCount);
                nCount++;
            }
            if(sTypeName=="Tetra2") {
                vLocalNode.resize(pmw::NumberOfNode::Tetra2());
                iss >> nElementID;
                for(i=0; i< vLocalNode.size(); i++) {
                    iss >> vLocalNode[i];
                }
                mpFactory->GeneElement(mgLevel, nMeshID, nElementID, pmw::ElementType::Tetra2, vLocalNode);
                mpFactory->setIDBucketElement(mgLevel, nMeshID, nElementID, nCount);
                nCount++;
            }
            if(sTypeName=="Hexa2") {
                vLocalNode.resize(pmw::NumberOfNode::Hexa2());
                iss >> nElementID;
                for(i=0; i< vLocalNode.size(); i++) {
                    iss >> vLocalNode[i];
                }
                mpFactory->GeneElement(mgLevel, nMeshID, nElementID, pmw::ElementType::Hexa2, vLocalNode);
                mpFactory->setIDBucketElement(mgLevel, nMeshID, nElementID, nCount);
                nCount++;
            }
            if(sTypeName=="Prism2") {
                vLocalNode.resize(pmw::NumberOfNode::Prism2());
                iss >> nElementID;
                for(i=0; i< vLocalNode.size(); i++) {
                    iss >> vLocalNode[i];
                }
                mpFactory->GeneElement(mgLevel, nMeshID, nElementID, pmw::ElementType::Prism2, vLocalNode);
                mpFactory->setIDBucketElement(mgLevel, nMeshID, nElementID, nCount);
                nCount++;
            }
            if(sTypeName=="Triangle2") {
                vLocalNode.resize(pmw::NumberOfNode::Triangle2());
                iss >> nElementID;
                for(i=0; i< vLocalNode.size(); i++) {
                    iss >> vLocalNode[i];
                }
                mpFactory->GeneElement(mgLevel, nMeshID, nElementID, pmw::ElementType::Triangle2, vLocalNode);
                mpFactory->setIDBucketElement(mgLevel, nMeshID, nElementID, nCount);
                nCount++;
            }
            if(sTypeName=="Quad2") {
                vLocalNode.resize(pmw::NumberOfNode::Quad2());
                iss >> nElementID;
                for(i=0; i< vLocalNode.size(); i++) {
                    iss >> vLocalNode[i];
                }
                mpFactory->GeneElement(mgLevel, nMeshID, nElementID, pmw::ElementType::Quad2, vLocalNode);
                mpFactory->setIDBucketElement(mgLevel, nMeshID, nElementID, nCount);
                nCount++;
            }
            if(sTypeName=="Beam2") {
                vLocalNode.resize(pmw::NumberOfNode::Beam2());
                iss >> nElementID;
                for(i=0; i< vLocalNode.size(); i++) {
                    iss >> vLocalNode[i];
                }
                mpFactory->GeneElement(mgLevel, nMeshID, nElementID, pmw::ElementType::Beam2, vLocalNode);
                mpFactory->setIDBucketElement(mgLevel, nMeshID, nElementID, nCount);
                nCount++;
            }
        };
        mpFactory->setupElement(mgLevel, nMeshID);
        return true;
    } else {
        return false;
    }
}
bool CFileReaderElement::Read_bin(ifstream& ifs)
{
    return true;
}
