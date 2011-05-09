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

    //debug
    cout << "~CGMGModel" << endl;
}

//
//
void CGMGModel::initAssyModel(const uint& num_of_mglevel)
{ 
    //cout << " CGMGModel::mvAssyModel.empty()    == " << mvAssyModel.empty() << endl;
    //cout << " CGMGModel::mvAssyModel.size()     == " << mvAssyModel.size() << endl;

    mvAssyModel.resize(num_of_mglevel);
    mNumOfLevel = num_of_mglevel;

    //debug
    uint num = mvAssyModel.size();
    mpLogger->Info(Utility::LoggerMode::MWDebug, "mvAssyModel.size() == ", num);
}

//
//
void CGMGModel::reserveAssyModel(const uint& num_of_mglevel)
{
    if(num_of_mglevel > mvAssyModel.max_size()) mvAssyModel.reserve(num_of_mglevel);
}

//
//
void CGMGModel::setModel(CAssyModel* pAssyModel, const uint& i)
{
    mvAssyModel[i] = pAssyModel;
}

//
//
void CGMGModel::addModel(CAssyModel* pAssyModel, const uint& i)
{
    if( i < mvAssyModel.size()){
        mvAssyModel[i]=pAssyModel;
    }else{
        mvAssyModel.push_back(pAssyModel);
    }
}

//
//
void CGMGModel::resizeAssyModel(const uint &num_of_part)
{
    mvAssyModel.resize(num_of_part);
}




