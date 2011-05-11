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
using namespace boost;

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
    uiint numOfContactMesh, maxID, minID;
    uiint my_rank, trans_rank;
    uiint nProp(0);//属性: 0:MPC, 1:接触
    uiint contactID, meshID, nodeID, rank, elemID, elemFaceID, shapeType;
    bool bmesh(false);//計算領域のContactMeshか？
    uiint maslave;//マスター,スレーブ切り替え
    uiint mgLevel(0);
    uiint conNodeID, numOfConNode, contactFaceID, numOfContactFace;
    
    vuint vConNodeID;
    vuint quadConNodeID;  quadConNodeID.resize(4);
    vuint quad2ConNodeID; quad2ConNodeID.resize(8);
    vuint triConNodeID;   triConNodeID.resize(3);
    vuint tri2ConNodeID;  tri2ConNodeID.resize(6);
    vuint beamConNodeID;  beamConNodeID.resize(2);
    vuint beam2ConNodeID; beam2ConNodeID.resize(3);
    vuint pointConNodeID; pointConNodeID.resize(1);
    
    vdouble vCoord; vCoord.resize(3);
    string shapeStr;//文字列"Quad","Triangle","Beam","Point"読み込み用途
    vstring strData; strData.resize(3);//文字列"-"の読み込み用途

    string sParamType;//ConNodeのパラメータタイプ(変位,スカラー,変位+スカラー)
    uiint numOfDisp,numOfScalar;

    if(TagCheck(sLine, FileBlockName::StartContactMesh()) ){// スタート タグ

        //debug
        //cout << "FileReaderContactMesh::Read" << endl;

        while(true){
            sLine = getLineSt(ifs);
            if(TagCheck(sLine, FileBlockName::EndContactMesh()) ) break; // エンド タグ
            
            istringstream iss(sLine.c_str());
            // 接合メッシュ数(MPCメッシュ数),MaxID,MinID
            iss >> numOfContactMesh >> maxID >> minID;// maxID,minIDは,使ってない…

            //cout << "FileReaderContactMesh::Read  --- A" << endl;
            
            uiint icont;
            for(icont=0; icont< numOfContactMesh; icont++){
                sLine= getLineSt(ifs);
                iss.clear();
                iss.str(sLine);

                // boost トークン分割
                // ----
                char_separator<char> sep(" \t\n");
                tokenizer< char_separator<char> > tokens(sLine, sep);

                uiint nCount(0);
                typedef tokenizer< char_separator<char> >::iterator Iter;
                for(Iter it=tokens.begin(); it != tokens.end(); ++it){
                    string str = *it;
                    if(nCount==0){ contactID = atoi(str.c_str());}
                    if(nCount==1){ my_rank   = atoi(str.c_str());}
                    if(nCount==2){ trans_rank= atoi(str.c_str());}
                    if(nCount==3){ nProp     = atoi(str.c_str());}//入力ファイルにnPropが無い場合は、"初期値0:MPC"のまま.
                    nCount++;
                };

                //cout << "FileReaderContactMesh::Read  --- B" << endl;
                
                
                mpFactory->GeneContactMesh(contactID, my_rank, trans_rank, nProp);//nPropが無い場合は、初期値0:MPC

                sLine= getLineSt(ifs);
                iss.clear();
                iss.str(sLine);

                // ContactNode数,MaxID,MinID
                iss >> numOfConNode >> maxID >> minID;// maxID,minIDは,Factoryで使っていない…

                uiint icnode;
                for(icnode=0; icnode< numOfConNode; icnode++){
                    sLine= getLineSt(ifs);
                    iss.clear();
                    iss.str(sLine);

                    // contactNodeID, X, Y, Z, meshID, nodeID, rank, マスタースレーブ, paramタイプ, 変位DOF, スカラー数
                    iss >> conNodeID >> vCoord[0] >> vCoord[1] >> vCoord[2] >> strData[0] >> strData[1] >> rank >> maslave;
                    iss >> sParamType >> numOfDisp >> numOfScalar;
                    

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

                //cout << "FileReaderContactMesh::Read  --- C" << endl;
                
                // maslave::マスター,スレーブ切り替え
                for(maslave=0; maslave< 2; maslave++){
                    sLine= getLineSt(ifs);
                    iss.clear();
                    iss.str(sLine);
                    // 接合面の数,MaxID,MinID
                    iss >> numOfContactFace >> maxID >> minID;

                    uiint iface;
                    for(iface=0; iface< numOfContactFace; iface++){
                        sLine= getLineSt(ifs);
                        iss.clear();
                        iss.str(sLine);

                        // 接合面ID,meshID,elemID,elemFaceID, shapeType, conNodeID.......
                        iss >> contactFaceID >> strData[0]  >> strData[1] >> strData[2] >> shapeStr;


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
                        
                        shapeType = IntElemType(shapeStr);
                        
                        switch(shapeType){
                            case(pmw::ElementType::Quad):
                                iss >> quadConNodeID[0] >> quadConNodeID[1] >> quadConNodeID[2] >> quadConNodeID[3];
                                vConNodeID= quadConNodeID;
                                break;
                                
                            case(pmw::ElementType::Quad2):
                                iss >> quad2ConNodeID[0] >> quad2ConNodeID[1] >> quad2ConNodeID[2] >> quad2ConNodeID[3]
                                    >> quad2ConNodeID[4] >> quad2ConNodeID[5] >> quad2ConNodeID[6] >> quad2ConNodeID[7];
                                vConNodeID= quad2ConNodeID;
                                break;
                                
                            case(pmw::ElementType::Triangle):
                                iss >> triConNodeID[0] >> triConNodeID[1] >> triConNodeID[2];
                                vConNodeID= triConNodeID;
                                break;
                                
                            case(pmw::ElementType::Triangle2):
                                iss >> tri2ConNodeID[0] >> tri2ConNodeID[1] >> tri2ConNodeID[2]
                                    >> tri2ConNodeID[3] >> tri2ConNodeID[4] >> tri2ConNodeID[5];
                                vConNodeID= tri2ConNodeID;
                                break;
                                
                            case(pmw::ElementType::Beam):
                                iss >> beamConNodeID[0] >> beamConNodeID[1];
                                vConNodeID= beamConNodeID;
                                break;
                                
                            case(pmw::ElementType::Beam2):
                                iss >> beam2ConNodeID[0] >> beam2ConNodeID[1] >> beam2ConNodeID[2];
                                vConNodeID= beam2ConNodeID;
                                break;
                                
                            case(pmw::ElementType::Point):
                                iss >> pointConNodeID[0];
                                vConNodeID= pointConNodeID;
                                break;
                                
                            default:
                                break;
                        }
                        

                        //Face rank
                        uiint face_rank;
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


bool CFileReaderContactMesh::Read_bin(ifstream& ifs)
{
    return true;
}






