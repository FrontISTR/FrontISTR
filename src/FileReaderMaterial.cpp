//
//  FileReaderMaterial.cpp
//
//
//
//                  2009.07.27
//                  2009.07.27
//                  k.Takeda
#include "FileReaderMaterial.h"
using namespace FileIO;

// construct & destruct
//
CFileReaderMaterial::CFileReaderMaterial()
{
    ;
}
CFileReaderMaterial::~CFileReaderMaterial()
{
    ;
}

// method
// --
bool CFileReaderMaterial::Read(ifstream& ifs, string& sLine)
{
    string sMaterialName;
    uint   materialID, numOfMaterial, numOfItem;
    vstring vPropName;
    vuint   vPropType;
    vdouble vPropValue;
    

    if(TagCheck(sLine, FileBlockName::StartMaterial()) ){
        
        sLine = getLineSt(ifs);
        istringstream iss(sLine.c_str());
        
        // 材料数,  MeshID <= MeshIDについては検討段階
        iss >> numOfMaterial;

        mpFactory->reserveMaterial(numOfMaterial);
        
        while(true){
            sLine = getLineSt(ifs);

            if(TagCheck(sLine, FileBlockName::EndMaterial())) break;

            // MaterialID, 材料名, アイテム数, Material Value (アイテム数だけ存在)
            //
            istringstream iss(sLine.c_str());
            iss >> materialID >> sMaterialName >> numOfItem;
            vPropName.resize(numOfItem);
            vPropType.resize(numOfItem);
            vPropValue.resize(numOfItem);
            

            for(uint i=0; i< numOfItem; i++){
                iss >> vPropName[i] >> vPropValue[i];
                
                //プロパティ名による,MateriapPropType選別
                //
                if(vPropName[i]=="Density") vPropType[i]= pmw::MaterialPropType::Density;        //ρ(密度)
                if(vPropName[i]=="Poisson") vPropType[i]= pmw::MaterialPropType::Poisson;        //ν(ポアソン比)
                if(vPropName[i]=="YoungModule") vPropType[i]= pmw::MaterialPropType::YoungModule;//E(縦弾性係数)
                if(vPropName[i]=="Temp_Depend_YoungModule") vPropType[i]= pmw::MaterialPropType::Temp_Depend_YoungModule;//温度依存 E(縦弾性係数)
                if(vPropName[i]=="Temp_Depend_Poisson") vPropType[i]= pmw::MaterialPropType::Temp_Depend_Poisson;        //温度依存 ν(ポアソン比)
                if(vPropName[i]=="Linear_Expansion") vPropType[i]= pmw::MaterialPropType::Linear_Expansion;        //線膨張
                if(vPropName[i]=="Thermal_Conductivity") vPropType[i]= pmw::MaterialPropType::Thermal_Conductivity;//熱伝導   
                if(vPropName[i]=="Heat_Transfer_Rate") vPropType[i]= pmw::MaterialPropType::Heat_Transfer_Rate;//熱伝達
                if(vPropName[i]=="Cp") vPropType[i]= pmw::MaterialPropType::Cp;                                //定圧比熱
                if(vPropName[i]=="Cv") vPropType[i]= pmw::MaterialPropType::Cv;                                //定積比熱

                mpLogger->Monitor(Utility::LoggerMode::MWDebug,vPropType[i],vPropValue[i],vPropName[i]);
            };
            
            mpFactory->GeneMaterial(materialID, sMaterialName, vPropType, vPropValue);

            vPropName.clear();
            vPropType.clear();
            vPropValue.clear();
        };
        return true;
    }else{
        return false;
    }
}










