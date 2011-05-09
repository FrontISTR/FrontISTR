//
//  FileReaderContactMesh.cpp
//
//
//
//                      2009.10.20
//                      2009.10.20
//                      k.Takeda
#include "FileReaderContactMesh.h"
using namespace FileIO;

// construct & destruct
//
CFileReaderContactMesh::CFileReaderContactMesh()
{
    ;
}
CFileReaderContactMesh::~CFileReaderContactMesh()
{
    ;
}

// 読み込み メソッド
// --
bool CFileReaderContactMesh::Read(ifstream& ifs, string& sLine)
{
    uint numOfContactMesh, maxID, minID;
    uint my_rank, trans_rank;
    uint contactID, meshID, nodeID, rank, elemID, elemFaceID, shapeType;
    bool bmesh(false);//計算領域のContactMeshか？
    uint maslave;//マスター,スレーブ切り替え
    uint mgLevel(0);
    uint conNodeID, numOfConNode, contactFaceID, numOfContactFace;
    
    vuint vConNodeID;
    vuint quadConNodeID; quadConNodeID.resize(4);
    vuint triConNodeID;  triConNodeID.resize(3);
    vuint beamConNodeID; beamConNodeID.resize(2);
    vuint pointConNodeID;pointConNodeID.resize(1);
    
    vdouble vCoord; vCoord.resize(3);
    string shapeStr;//文字列"Quad","Triangle","Beam","Point"読み込み用途
    vstring strData; strData.resize(3);//文字列"-"の読み込み用途

    string sParamType;//ConNodeのパラメータタイプ(変位,スカラー,変位+スカラー)
    uint numOfDisp,numOfScalar;

    if(TagCheck(sLine, FileBlockName::StartContactMesh()) ){// スタート タグ

        //debug
        cout << "FileReaderContactMesh::Read" << endl;

        while(true){
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndContactMesh()) ) break; // エンド タグ
            
            istringstream iss(sLine.c_str());
            // 接合メッシュ数(MPCメッシュ数),MaxID,MinID
            iss >> numOfContactMesh >> maxID >> minID;// maxID,minIDは,使ってない…
            
            uint icont;
            for(icont=0; icont< numOfContactMesh; icont++){
                sLine= getLineSt(ifs);
                iss.clear();
                iss.str(sLine);
                // 接合MeshID::ContactMeshID, rank(接合Meshのrank)
                iss >> contactID >> my_rank >> trans_rank;
                
                mpFactory->GeneContactMesh(contactID, my_rank, trans_rank);

                sLine= getLineSt(ifs);
                iss.clear();
                iss.str(sLine);

                // ContactNode数,MaxID,MinID
                iss >> numOfConNode >> maxID >> minID;// maxID,minIDは,Factoryで使っていない…

                uint icnode;
                for(icnode=0; icnode< numOfConNode; icnode++){
                    sLine= getLineSt(ifs);
                    iss.clear();
                    iss.str(sLine);

                    // contactNodeID, X, Y, Z, meshID, nodeID, rank, マスタースレーブ, paramタイプ, 変位DOF, スカラー数
                    iss >> conNodeID >> vCoord[0] >> vCoord[1] >> vCoord[2] >> strData[0] >> strData[1] >> rank >> maslave;
                    iss >> sParamType >> numOfDisp >> numOfScalar;

//                    //debug
//                    cout << "FileReaderContactMesh::Read, conNodeID= " << conNodeID
//                            << ", x=" << vCoord[0] << ", y=" << vCoord[1] << ", z=" << vCoord[2]
//                            << ", strData[0]= " << strData[0] << ", strData[1]= " << strData[1] << endl;
                    

                    //対応するデータが存在する場合は("-"ではない場合),データを文字列から整数に変換
                    bmesh=false;//フラグ初期化
                    if(strData[0]!="-" && strData[1]!="-"){
                        bmesh= true;
                        meshID = boost::lexical_cast<unsigned int>(strData[0]);
                        nodeID = boost::lexical_cast<unsigned int>(strData[1]);
                    }
                    //const string& s_param_type, const uint& numOfVector, const uint& numOfScalar,
                    mpFactory->GeneContactNode(mgLevel, contactID, conNodeID, vCoord, 
                                                sParamType, numOfDisp, numOfScalar,
                                                bmesh, meshID, nodeID, rank, maslave);
                    
                };//icnodeループ(ContactNodeループ)
                
                // maslave::マスター,スレーブ切り替え
                for(maslave=0; maslave< 2; maslave++){
                    sLine= getLineSt(ifs);
                    iss.clear();
                    iss.str(sLine);
                    // 接合面の数,MaxID,MinID
                    iss >> numOfContactFace >> maxID >> minID;

                    uint iface;
                    for(iface=0; iface< numOfContactFace; iface++){
                        sLine= getLineSt(ifs);
                        iss.clear();
                        iss.str(sLine);

                        // 接合面ID,meshID,elemID,elemFaceID, shapeType, conNodeID.......
                        iss >> contactFaceID >> strData[0]  >> strData[1] >> strData[2] >> shapeStr;

//                        //debug
//                        cout << "FileReaderContactMesh::Read, ConFaceID= " << contactFaceID
//                                << ", strData[0]= " << strData[0]
//                                << ", strData[1]= " << strData[1] << ", strData[2]= " << strData[2] << endl;


                        //対応するデータが存在する場合は("-"ではない場合),データを文字列から整数に変換
                        bmesh=false;//フラグ初期化
                        if(strData[0]!="-" && strData[1]!="-" && strData[2]!="-"){
                            bmesh= true;
                            meshID =  boost::lexical_cast<unsigned int>(strData[0]);
                            elemID =  boost::lexical_cast<unsigned int>(strData[1]);
                            elemFaceID = boost::lexical_cast<unsigned int>(strData[2]);
                        }
                        // ContactNodeIDコネクティビティ 読み込み
                        vConNodeID.clear();
                        if(shapeStr=="Quad"){
                            shapeType= pmw::ElementType::Quad;
                            iss >> quadConNodeID[0] >> quadConNodeID[1] >> quadConNodeID[2] >> quadConNodeID[3];

                            vConNodeID= quadConNodeID;
                        }
                        if(shapeStr=="Triangle"){
                            shapeType= pmw::ElementType::Triangle;
                            iss >> triConNodeID[0] >> triConNodeID[1] >> triConNodeID[2];

                            vConNodeID= triConNodeID;
                        }
                        if(shapeStr=="Beam"){
                            shapeType= pmw::ElementType::Beam;
                            iss >> beamConNodeID[0] >> beamConNodeID[1];

                            vConNodeID= beamConNodeID;
                        }
                        if(shapeStr=="Point"){
                            shapeType= pmw::ElementType::Point;
                            iss >> pointConNodeID[0];

                            vConNodeID= pointConNodeID;
                        }
                        uint face_rank;
                        //Face rank
                        iss >> face_rank;
                        
                        switch(maslave){//マスタースレーブswitch
                            case(0)://マスター面 生成
                                mpFactory->GeneMasterFace(contactID, shapeType, contactFaceID, bmesh, 
                                        meshID,elemID,elemFaceID,vConNodeID, face_rank);
                                break;
                            case(1)://スレーブ面 生成
                                mpFactory->GeneSlaveFace(contactID,shapeType,contactFaceID, bmesh, 
                                        meshID,elemID,elemFaceID, vConNodeID, face_rank);
                                break;
                        }//switch(maslave)

                    };//ifaceループ
                };//maslave:マスター,スレーブ切り替えループ
            };//icontループ(ContactMeshループ)
        };//while ループ(読み込みループ)
        return true;
    }else{
        return false;
    }
}









