//
// GMGModel.cpp
//				2009.05.01
//				2008.11.10
//				k.Takeda
#include "GMGModel.h"
using namespace pmw;

CGMGModel::CGMGModel(void)
{
    // mvAssyModel.resize(0);
    mpLogger = Utility::CLogger::Instance();
}

CGMGModel::~CGMGModel(void)
{
    for_each(mvAssyModel.begin(), mvAssyModel.end(), DeleteObject());

    ////debug
    //cout << "~CGMGModel" << endl;
}

//
// FileReaderRefine -> Factory::GeneAssyModelから使用
//
void CGMGModel::initAssyModel(const uiint& num_of_mglevel)
{ 
    //cout << " CGMGModel::mvAssyModel.empty()    == " << mvAssyModel.empty() << endl;
    //cout << " CGMGModel::mvAssyModel.size()     == " << mvAssyModel.size() << endl;

    mvAssyModel.resize(num_of_mglevel);
    mNumOfLevel = num_of_mglevel;

    ////debug
    //uiint num = mvAssyModel.size();
    //mpLogger->Info(Utility::LoggerMode::MWDebug, "mvAssyModel.size() == ", num);
}

//
// FileReaderRefine -> Factory::GeneAssyModelから使用
//
void CGMGModel::reserveAssyModel(const uiint& num_of_mglevel)
{
    if(num_of_mglevel > mvAssyModel.max_size()) mvAssyModel.reserve(num_of_mglevel);
}

//
//
void CGMGModel::setModel(CAssyModel* pAssyModel, const uiint& i)
{
    mvAssyModel[i] = pAssyModel;
}

//
//
void CGMGModel::addModel(CAssyModel* pAssyModel, const uiint& i)
{
    if( i < mvAssyModel.size()){
        mvAssyModel[i]=pAssyModel;
    }else{
        mvAssyModel.push_back(pAssyModel);
    }
}

//
// CMW から初期化する場合の関数
//
void CGMGModel::initAssyModel()
{
    mvAssyModel.resize(1);

    mvAssyModel[0]= new CAssyModel;
    mvAssyModel[0]->setMGLevel(0);

    mNumOfLevel = 1;
}

//
// CMW::Refine -> Factory::GeneAssyModelから使用
//
void CGMGModel::resizeAssyModel(const uiint &nNumOfMGLevel)
{
    mvAssyModel.resize(nNumOfMGLevel);
    mNumOfLevel = nNumOfMGLevel;
}




