//
//	GMGModel.h
//
//				2009.05.01
//				2008.11.10
//				k.Takeda
#ifndef GMG_Model_hh_C77868BF_2116_41b3_B667_59BA392FCC5C
#define GMG_Model_hh_C77868BF_2116_41b3_B667_59BA392FCC5C

#include "CommonStd.h"
#include "AssyModel.h"

#include "Material.h"

#include "Logger.h"

namespace pmw{
class CGMGModel
{
////public:
////	static CGMGModel* Instance(){
////		if(!pInstance)
////		    pInstance = new CGMGModel();
////		return pInstance;
////	}
////private:
////	static CGMGModel* pInstance;
////	CGMGModel(void);
public:
static CGMGModel* Instance(){
	static CGMGModel gmgModel;
	return &gmgModel;
}
private:
	CGMGModel(void);
public:
    //CGMGModel(void);
    virtual ~CGMGModel(void);


protected:
    std::vector<CAssyModel*> mvAssyModel;
    uint  mNumOfLevel;//MultiGrid Level数(階層数)

    vector<CMaterial*> mvMaterial;//材質数だけのMaterialオブジェクト

    Utility::CLogger *mpLogger;

public:
    // AssemblyModel ( MultiGrid Level )
    CAssyModel* getAssyModel(const uint& mgLevel){ return mvAssyModel[mgLevel];}

    // MultiGrid Level
    //void  setNumOfLevel(const uint& num_of_level){ mNumOfLevel= num_of_level;}
    uint& getNumOfLevel(){return mNumOfLevel;}


    // GMG
    void  addModel(CAssyModel *pAssyModel, const uint& i);
    void  initAssyModel(const uint& num_of_mglevel);
    void  reserveAssyModel(const uint& num_of_mglevel);
    void  setModel(CAssyModel *pAssyModel, const uint& i);

    // test method
    void resizeAssyModel(const uint& num_of_part);



    // Material
    void reserveMaterial(const uint& res_size){ mvMaterial.reserve(res_size);}
    void setMaterial(CMaterial* pMaterial){ mvMaterial.push_back(pMaterial);}
    CMaterial* getMaterial(const uint& i){ return mvMaterial[i];}
    uint getNumOfMaterial(){ return mvMaterial.size();}

};
}
#endif
