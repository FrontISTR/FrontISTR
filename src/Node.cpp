//
// Node3.cpp
//
//				2008.12.18
//				2008.12.18
//				k.Takeda
#include "Node.h"
using namespace pmw;

// 
//
CNode::CNode(void)
{
//	mb_Boundary = false;
//	mb_DOF3 = false;
//	mb_DOF6 = false;
}

CNode::~CNode(void)
{
//    //debug
//    cout << "~CNode" << endl;
}

//// Node3
////
//void CNode::InitializeNode3()
//{
//    //    mvDisp.resize(3);  // tx,ty,tz
//    //    mvStress.resize(6);// x, y, z, xy, yz, zx
//    //
//    //    mvVonMieses.resize(1);//
//    //    mvTemperature.resize(1);// Temperature
//
//    mb_DOF3 = true;
//    mb_DOF6 = false;
//    mb_ADOF = false;
//}
//
//// Node6
////
//void CNode::InitializeNode6()
//{
//    //    mvDisp.resize(6); // tx,ty,tz, rx,ry,rz
//    //    mvStress.resize(5);// Mx,My,Mxy,Qx,Qy,Qz
//    //
//    //    mvVonMieses.resize(2);
//    //    mvTemperature.resize(1);// Temperature
//
//    mb_DOF3 = false;
//    mb_DOF6 = true;
//    mb_ADOF = false;
//}
//
//// Arbitrary DOF Node Initialize(任意自由度Nodeの変数確保)
////   vParam.size() == DOF
////   vParam[0]     == 0 paramの変数の個数
////   vParam[1]     == 1 paramの変数の個数
////     ....        ....      .....
////   vParam[#]     == # paramの変数の個数
////
//void CNode::InitializeNodeADOF(const vint& vParam, const uint& num_of_param)
//{
//    Utility::CLogger *pLogger = Utility::CLogger::Instance();
//    if(num_of_param != vParam.size()) pLogger->Info(Utility::LoggerMode::MWDebug,"Parameter Mismatch! @CNode::InitializeNodeADOF");
//
//    mvvArbitVar.resize(num_of_param);
//    uint i;
//    for(i=0; i < num_of_param; i++){
//        mvvArbitVar[i].resize(vParam[i]);
//    };
//    mb_DOF3 = false;
//    mb_DOF6 = false;
//    mb_ADOF = true;
//}

////
////
//void CNode::setBool(const char* type, const bool& b_type)
//{
//    Utility::CLogger *pLogger = Utility::CLogger::Instance();
//
//    if("Boundary"==type){
//        mb_Boundary = b_type;
//    }
//    else if("DOF3"==type){
//        mb_DOF3 = b_type;
//    }
//    else if("DOF6"==type){
//        mb_DOF6 = b_type;
//    }
//    else{
//        pLogger->Info(Utility::LoggerMode::Error,"Node::setBool arg ERROR");
//    }
//}
