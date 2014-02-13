/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderMaterial.cpp
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
#include "FileReaderMaterial.h"
using namespace FileIO;
CFileReaderMaterial::CFileReaderMaterial()
{
    ;
}
CFileReaderMaterial::~CFileReaderMaterial()
{
    ;
}
string CFileReaderMaterial::Name()
{
    return  "FileReaderMaterial";
}

bool CFileReaderMaterial::Read(ifstream& ifs, string& sLine)
{
    string sMaterialName;
    uiint   materialID, numOfMaterial, numOfItem, meshID;
    vstring vPropName;
    vuint   vPropType;
    vdouble vPropValue;
    if(TagCheck(sLine, FileBlockName::StartMaterial()) ) {
        sLine = getLineSt(ifs);
        istringstream iss(sLine.c_str());
        iss >> numOfMaterial;
        mpFactory->reserveMaterial(numOfMaterial);
        while(true) {
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndMaterial())) break;
            istringstream iss(sLine.c_str());
            iss >> meshID >> materialID >> sMaterialName >> numOfItem;
            vPropName.resize(numOfItem);
            vPropType.resize(numOfItem);
            vPropValue.resize(numOfItem);
            for(uiint i=0; i< numOfItem; i++) {
                iss >> vPropName[i] >> vPropValue[i];
                if(vPropName[i]=="Density") vPropType[i]= pmw::MaterialPropType::Density;
                if(vPropName[i]=="Poisson") vPropType[i]= pmw::MaterialPropType::Poisson;
                if(vPropName[i]=="YoungModule") vPropType[i]= pmw::MaterialPropType::YoungModule;
                if(vPropName[i]=="Temp_Depend_YoungModule") vPropType[i]= pmw::MaterialPropType::Temp_Depend_YoungModule;
                if(vPropName[i]=="Temp_Depend_Poisson") vPropType[i]= pmw::MaterialPropType::Temp_Depend_Poisson;
                if(vPropName[i]=="Linear_Expansion") vPropType[i]= pmw::MaterialPropType::Linear_Expansion;
                if(vPropName[i]=="Thermal_Conductivity") vPropType[i]= pmw::MaterialPropType::Thermal_Conductivity;
                if(vPropName[i]=="Heat_Transfer_Rate") vPropType[i]= pmw::MaterialPropType::Heat_Transfer_Rate;
                if(vPropName[i]=="Cp") vPropType[i]= pmw::MaterialPropType::Cp;
                if(vPropName[i]=="Cv") vPropType[i]= pmw::MaterialPropType::Cv;
                mpLogger->Monitor(Utility::LoggerMode::MWDebug,vPropType[i],vPropValue[i],vPropName[i]);
            };
            mpFactory->GeneMaterial(meshID, materialID, sMaterialName, vPropType, vPropValue);
            vPropName.clear();
            vPropType.clear();
            vPropValue.clear();
        };
        return true;
    } else {
        return false;
    }
}
bool CFileReaderMaterial::Read_bin(ifstream& ifs)
{
    return true;
}
