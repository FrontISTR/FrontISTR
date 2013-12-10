/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   FileWriter.cpp
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "FileWriter.h"
using namespace boost::xpressive;
using namespace boost;
CFileWriter::CFileWriter()
{
	msDirichlet= "Dirichlet";
	msNeumann  = "Neumann";
}
CFileWriter::~CFileWriter()
{
}

//--
// Header && mesh version 0.4
//--
bool CFileWriter::WriteHeader(ofstream& ofs)
{
	string sVer = MW3_Tag::VerNumber();// バージョン番号

	ofs << "//----------------------------------------      \n"
		     "//      converted mesh data                     \n"
		     "// output format  \"MW3 mesh ver " + sVer + "\" \n"
		     "//----------------------------------------      \n" << endl;
	
	//--
	// Version 
	//--
	ofs << "Version \n"      
	       "  MW3_" + sVer + " \n"
	       "End \n" 
			<< endl;

	return true;
}
//--
// AssyModel
//--
bool CFileWriter::WriteAssyModel(CAssyModel *pAssyModel, ofstream& ofs)
{
	if(!pAssyModel->existMesh()) return false;//Meshが無い.

	CMessage *pMsg=CMessage::Instance();
	
	/*
	メッシュパーツ数  最大ID　最小ID
	メッシュパーツID  属性番号　
	…
	通信テーブル数
	メッシュパーツID  ランク・ペア1st  ランク・ペア2nd
	…
	*/

	ofs << "AssyModel" << endl;
	//メッシュ数  MaxID  MinID
	ofs << format(" %lu %lu %lu") %pAssyModel->getNumOfMesh() %pAssyModel->getMaxIDinMesh() %pAssyModel->getMinIDinMesh() << endl;
	//MeshID  属性番号
	for(size_t i=0; i < pAssyModel->getNumOfMesh(); i++){
		CMesh* pMesh=pAssyModel->getMesh(i);
		ofs << format(" %3lu %lu") %pMesh->getID()  %0 << endl;
	};
	//通信テーブル数(シングル:常に"0")
	ofs << format(" %3lu") %0 << endl;
	ofs << "End\n" << endl;


	//メッセージ表示
	pMsg->info("wrote :1 block \"AssyModel\"");
	
	return true;
}
//--
// Node
//--
bool CFileWriter::WriteNode(CAssyModel *pAssyModel, ofstream& ofs)
{
	if(!pAssyModel->existMesh()) return false;//Meshが無い

	CMessage *pMsg=CMessage::Instance();

	size_t nNumMesh= pAssyModel->getNumOfMesh();

	for(size_t i=0; i < nNumMesh; i++){
		CMesh *pMesh= pAssyModel->getMesh(i);
		/*
		節点数  メッシュパーツID　最大節点ID  最小節点ID
		種類(V,S,SV)　スカラー変数の数  ベクトル変数の数  ID  X  Y  Z
		…
		*/
		ofs << "Node" << endl;
		ofs << format(" %lu %lu %lu %lu") %pMesh->getNumOfNode() %pMesh->getID() %pMesh->getMaxNodeID() %pMesh->getMinNodeID() << endl;
		
		for(size_t ii=0; ii < pMesh->getNumOfNode(); ii++){
			CNode *pNode= pMesh->getNode(ii);
			ofs << format(" %s  %9lu  %23.15e %23.15e %23.15e") %"V 0 3"  %pNode->getID() %pNode->getX() %pNode->getY() %pNode->getZ() << endl;
		};
		ofs << "End\n" << endl;
	};


	// メッセージ表示
	string sNumMesh= lexical_cast<string>(nNumMesh);
	pMsg->info("wrote :"+sNumMesh+" block \"Node\"");

	return true;
}
//--
// Element
//--
bool CFileWriter::WriteElement(CAssyModel *pAssyModel, ofstream& ofs)
{
	if(!pAssyModel->existMesh()) return false;//Meshが無い

	CMessage *pMsg=CMessage::Instance();

	// regex
	string strRegex= FistrElemTypeS::regStandard();// MW3に存在する要素タイプ
	sregex regBasis = sregex::compile(strRegex);
	string strRegex2= FistrElemTypeS::regShellType();// シェル要素
	sregex regShell= sregex::compile(strRegex2);
	smatch res;

	size_t nNumMesh= pAssyModel->getNumOfMesh();

	for(size_t i=0; i < nNumMesh; i++){
		CMesh *pMesh= pAssyModel->getMesh(i);
		/*
		要素数  メッシュパーツID　最大要素ID  最小要素ID
		種類(Hexa,Tetra,etc.)　ID  節点ID  節点ID  節点ID  節点ID ... ...(要素節点数だけ繰り返し)
		…
		*/
		ofs << "Element" << endl;
		ofs << format(" %lu %lu %lu %lu") %pMesh->getNumOfElement() %pMesh->getID() %pMesh->getMaxElementID() %pMesh->getMinElementID() << endl;
		for(size_t ii=0; ii < pMesh->getNumOfElement(); ii++){
			
			CElement* pElem= pMesh->getElement(ii);
			string sType=pElem->getCType();
			size_t nType=pElem->getNType();

			//if( regex_search(sType, res, regBasis) ){
			if( regex_match(sType, regBasis)){
				//要素タイプ　ID  節点ID  節点ID  節点ID  節点ID ... ...
				//ofs << format(" %s %9lu ") %res.str()  %pElem->getID();
				ofs << format(" %s %9lu ") %sType  %pElem->getID();
				for(size_t iii=0; iii < pElem->getNumOfNode(); iii++){
					ofs << format(" %9lu")  %pElem->getNodeID_Fistr2MW3(iii);
				};
				ofs << endl;
			}else if( regex_search(sType, res, regShell) ){
				//シェル要素をQuad,Triangleに変換
				string sRepType;
				switch(nType){
				case(FistrElementType::TriShell):   sRepType=FistrElemTypeS::Triangle();   break;
				case(FistrElementType::TriShell2):  sRepType=FistrElemTypeS::Triangle2();  break;
				case(FistrElementType::QuadShell):  sRepType=FistrElemTypeS::Quad();  break;
				case(FistrElementType::QuadShell2): sRepType=FistrElemTypeS::Quad2(); break;
				default:
					pMsg->error(res.str()+" ? ");//regexが正しければ、この処理には入らない.
					break;
				}
				//要素タイプ(置換)　ID  節点ID  節点ID  節点ID  節点ID ... ...
				ofs << format(" %s %9lu ") %sRepType %pElem->getID();
				for(size_t iii=0; iii < pElem->getNumOfNode(); iii++){
					ofs << format(" %9lu") %pElem->getNodeID_Fistr2MW3(iii);
				};
				ofs << endl;
			}else{
				pMsg->error( res.str()+" not supported");//Error --- IFaceQuadも含む.
			}
		};
		ofs << "End\n" << endl;
	};

	// メッセージ表示
	string sNumMesh= lexical_cast<string>(nNumMesh);
	pMsg->info("wrote :"+sNumMesh+" block \"Element\"");

	return true;
}
//--
// BoundaryNodeMesh : デフォルトはDirichlet型で出力(理由:FrontISTRは境界値を使用しない)
//--
bool CFileWriter::WriteBndNodeMesh(CAssyModel *pAssyModel, ofstream& ofs)
{
	if(!pAssyModel->existMesh()) return false;//Meshが無い

	CMessage *pMsg=CMessage::Instance();

	size_t nNMesh=pAssyModel->getNumOfMesh();
	size_t nSumGrpN=0;
	string sBndType;

	for(size_t i=0; i < nNMesh; i++){
		CMesh *pMesh=pAssyModel->getMesh(i);
		/*
		タグ：BoundaryNodeMesh
		書式
		メッシュパーツID  境界条件数
		節点境界メッシュID  境界条件種類　名称 数式数(DOF数) DOF 数式 DOF 数式 ... # 数式:ver0.4以降
		*/
		size_t nNNgrp=pMesh->getNumOfNgrp();
		nSumGrpN += nNNgrp;

		if(nNNgrp > 0){
			ofs << "BoundaryNodeMesh" << endl;
			ofs << format(" %lu %lu") %pMesh->getID() %nNNgrp << endl;
		
			for(size_t ii=0; ii < nNNgrp; ii++){
				CNgroup *pNgrp= pMesh->getNgrp(ii);

				switch(pNgrp->getBndType()){ //境界種類あり？
				case( BndType::Dirichlet ): sBndType = msDirichlet; break;
				case( BndType::Neumann ): sBndType = msNeumann; break;
				case( BndType::NotUse ): sBndType = msDirichlet; break;
				default: 
					sBndType = msDirichlet;
					break;
				}
			
				ofs << format(" %lu %s %s %lu") %pNgrp->getID() %sBndType %pNgrp->getGroupName() %0 << endl;//数式数:0(元データに無い)
			};
			ofs << "End\n" << endl;
		}
	};

	
	for(size_t i=0; i < nNMesh; i++){
		CMesh *pMesh=pAssyModel->getMesh(i);
		/*
		タグ：BoundaryNode
		書式( ディレクレ&ノイマン 共通 )
		節点境界メッシュID　境界条件種類  メッシュパーツID  境界節点数
		境界節点ID  節点ID  X  Y  Z  自由度番号 境界値(FrontISTRは値を使用しない⇒常に0.0)
		*/
		size_t nNNgrp=pMesh->getNumOfNgrp();
		for(size_t ii=0; ii < nNNgrp; ii++){
			CNgroup *pNgrp= pMesh->getNgrp(ii);

			switch(pNgrp->getBndType()){ //境界種類あり？
			case( BndType::Dirichlet ): sBndType = msDirichlet; break;
			case( BndType::Neumann ): sBndType = msNeumann; break;
			case( BndType::NotUse ): sBndType = msDirichlet; break;
			default: 
				sBndType = msDirichlet;
				break;
			}

			ofs << "BoundaryNode" << endl;
			ofs << format(" %lu %s %lu %lu") %pNgrp->getID() %sBndType %pMesh->getID() %pNgrp->getNumOfNode() << endl;
			for(size_t iii=0; iii < pNgrp->getNumOfNode(); iii++){
				CNode* pNode=pNgrp->getNode(iii);
				ofs << format(" %9lu %9lu %23.15e %23.15e %23.15e %lu   %lf") %iii %pNode->getID() %pNode->getX() %pNode->getY() %pNode->getZ() %0 %0.0 << endl;
			};
			ofs << "End\n" << endl;
		};
	};


	// メッセージ表示
	string sSumGrpN= lexical_cast<string>(nSumGrpN);
	pMsg->info("wrote :"+ sSumGrpN +" block \"BoundaryNode\"");

	return true;
}
//--
// BoundaryFaceMesh : デフォルトはNeumann型で出力(理由:FrontISTRは境界値を使用しない)
//--
bool CFileWriter::WriteBndFaceMesh(CAssyModel *pAssyModel, ofstream& ofs)
{
	if(!pAssyModel->existMesh()) return false;//Meshが無い

	CMessage *pMsg=CMessage::Instance();

	size_t nNMesh=pAssyModel->getNumOfMesh();
	size_t nSumGrpN=0;
	string sBndType;
	/*
	タグ：BoundaryFaceMesh
	書式
	メッシュパーツID  境界条件数
	面境界メッシュID  境界条件種類　名称  自由度数 自由度番号 自由度番号 ...(繰り返し) 数式数(DOF数) DOF 数式 DOF 数式 ...
	*/
	for(size_t i=0; i < nNMesh; i++){
		CMesh *pMesh=pAssyModel->getMesh(i);

		size_t nNSgrp=pMesh->getNumOfSgrp();
		nSumGrpN += nNSgrp;

		if(nNSgrp > 0){
			ofs << "BoundaryFaceMesh" << endl;
			ofs << format(" %lu %lu") %pMesh->getID() %nNSgrp << endl;
			for(size_t ii=0; ii < nNSgrp; ii++){
				CSgroup* pSgrp=pMesh->getSgrp(ii);

				switch(pSgrp->getBndType()){ //境界種類あり？
				case( BndType::Dirichlet ): sBndType = msDirichlet; break;
				case( BndType::Neumann ): sBndType = msNeumann; break;
				case( BndType::NotUse ): sBndType = msNeumann; break;
				default: 
					sBndType = msNeumann;
					break;
				}

				ofs << format(" %lu %s %s %lu %lu %lu") %pSgrp->getID() %sBndType %pSgrp->getGroupName() %1 %0 %0 << endl;//自由度数:1,自由度:0,数式数:0
			};
			ofs << "End\n" << endl;
		}
	};
	//--
	// BoundaryFace
	//--
	for(size_t i=0; i < nNMesh; i++){
		CMesh *pMesh=pAssyModel->getMesh(i);
		size_t nNSgrp=pMesh->getNumOfSgrp();
		
		for(size_t ii=0; ii < nNSgrp; ii++){
			CSgroup *pSgrp= pMesh->getSgrp(ii);

			if(pSgrp->getBndType()!=BndType::Dirichlet){
				/*
				タグ：BoundaryFace
				書式( ノイマン型 )
				面境界メッシュID　境界条件種類  メッシュパーツID  境界節点数 Face数
				境界節点ID  節点ID  
				…
				面形状タイプ 境界面ID  要素ID  面番号  自由度番号 境界節点ID...(繰り返し)  境界値
				*/
				ofs << "BoundaryFace" << endl;
				ofs << format(" %lu %s %5lu %5lu %5lu") %pSgrp->getID() %msNeumann %pMesh->getID() %pSgrp->getNumOfNode() %pSgrp->getNumOfFace() << endl;
				// BNode 
				size_t nNNode=pSgrp->getNumOfNode();
				for(size_t inode=0; inode < nNNode; inode++){
					CNode* pNode = pSgrp->getNode(inode);
					ofs << format(" %5lu %5lu") %inode %pNode->getID() << endl;
				};

				// Face
				size_t nNFace=pSgrp->getNumOfFace();
				for(size_t iface=0; iface < nNFace; iface++){
				
					ofs << format(" %s %5lu %5lu %5lu %5lu") %pSgrp->getFaceType(iface) %iface %pSgrp->getElemID(iface) %pSgrp->getConvMW3FaceN(iface) %0;

					vector<size_t> vBNode = pSgrp->getBNodeN_Face(iface);
					for(size_t ibnode=0; ibnode < vBNode.size(); ibnode++) ofs << format(" %5lu") %vBNode[ibnode];
					ofs << format("   %10.3e") %0.000 << endl;//境界値
				};
				ofs << "End\n" << endl;
			}else{
				/*
				タグ：BoundaryFace
				書式( ディレクレ型 )
				面境界メッシュID　境界条件種類  メッシュパーツID  境界節点数 Face数
				境界節点ID  節点ID  DOF数　 DOF番号 境界値 DOF番号 境界値 DOF番号 境界値…
				…
				面形状タイプ 境界面ID  要素ID  面番号  　 境界節点ID...(繰り返し)  
				…
				*/
				ofs << "BoundaryFace" << endl;
				ofs << format(" %lu %s %5lu %5lu %5lu") %pSgrp->getID() %msDirichlet %pMesh->getID() %pSgrp->getNumOfNode() %pSgrp->getNumOfFace() << endl;
				// BNode 
				size_t nNNode=pSgrp->getNumOfNode();
				for(size_t inode=0; inode < nNNode; inode++){
					CNode* pNode = pSgrp->getNode(inode);
					ofs << format(" %5lu %5lu %lu %lu %lf") %inode %pNode->getID() %1 %0 %0.000 << endl;
				};
				// Face
				size_t nNFace=pSgrp->getNumOfFace();
				for(size_t iface=0; iface < nNFace; iface++){
				
					ofs << format(" %s %5lu %5lu %5lu") %pSgrp->getFaceType(iface) %iface %pSgrp->getElemID(iface) %pSgrp->getConvMW3FaceN(iface);

					vector<size_t> vBNode = pSgrp->getBNodeN_Face(iface);
					for(size_t ibnode=0; ibnode < vBNode.size(); ibnode++) ofs << format(" %5lu") %vBNode[ibnode];
					ofs << endl;
				};
				ofs << "End\n" << endl;
			}//if(BndType)
		};
	};


	// メッセージ表示
	string sSumGrpN= lexical_cast<string>(nSumGrpN);
	pMsg->info("wrote :"+ sSumGrpN +" block \"BoundaryFace\"");

	return true;
}
//--
// BoundaryEdgeMesh : デフォルトはNeumann型で出力(理由:FrontISTRは境界値を使用しない)
//--
bool CFileWriter::WriteBndEdgeMesh(CAssyModel *pAssyModel, ofstream& ofs)
{
	if(!pAssyModel->existMesh()) return false;//Meshが無い

	CMessage *pMsg=CMessage::Instance();

	size_t nNMesh=pAssyModel->getNumOfMesh();
	size_t nSumGrpN=0;
	string sBndType;
	/*
	タグ：BoundaryEdgeMesh
	書式
	メッシュパーツID  境界条件数
	辺境界メッシュID  境界条件種類　名称  自由度数 自由度番号 自由度番号 ...(繰り返し) 数式数(DOF数) DOF 数式 DOF 数式 ...
	…
	*/
	for(size_t i=0; i < nNMesh; i++){
		CMesh *pMesh=pAssyModel->getMesh(i);

		size_t nNLgrp=pMesh->getNumOfLgrp();
		nSumGrpN += nNLgrp;

		if(nNLgrp > 0){
			ofs << "BoundaryEdgeMesh" << endl;
			ofs << format(" %lu %lu") %pMesh->getID() %nNLgrp << endl;
			for(size_t ii=0; ii < nNLgrp; ii++){
				CLgroup* pLgrp=pMesh->getLgrp(ii);

				switch(pLgrp->getBndType()){ //境界種類あり？
				case( BndType::Dirichlet ): sBndType = msDirichlet; break;
				case( BndType::Neumann ): sBndType = msNeumann; break;
				case( BndType::NotUse ): sBndType = msNeumann; break;
				default: 
					sBndType = msNeumann;
					break;
				}

				ofs << format(" %lu %s %s %lu %lu %lu") %pLgrp->getID() %sBndType %pLgrp->getGroupName() %1 %0 %0 << endl;//自由度数:1,自由度:0,数式数:0
			};
			ofs << "End\n" << endl;
		}
	};


	for(size_t i=0; i < nNMesh; i++){
		CMesh *pMesh=pAssyModel->getMesh(i);

		size_t nNLgrp=pMesh->getNumOfLgrp();
		
		for(size_t ii=0; ii < nNLgrp; ii++){
			CLgroup *pLgrp= pMesh->getLgrp(ii);
			if(pLgrp->getBndType()!=BndType::Dirichlet){
				/*
				タグ：BoundaryEdge
				書式 ( ノイマン )
				辺境界メッシュID　境界条件種類  メッシュパーツID  境界節点数　Edge数
				境界節点ID  節点ID 
				…
				辺形状タイプ 境界辺ID  要素ID  辺番号  自由度番号 境界節点ID 境界節点ID 境界値
				…
				*/
				ofs << "BoundaryEdge" << endl;
				ofs << format(" %lu %s %5lu %5lu %5lu") %pLgrp->getID() %msNeumann %pMesh->getID() %pLgrp->getNumOfNode() %pLgrp->getNumOfEdge() << endl;
				// BNode 
				size_t nNNode=pLgrp->getNumOfNode();
				for(size_t inode=0; inode < nNNode; inode++){
					CNode* pNode = pLgrp->getNode(inode);
					ofs << format(" %5lu %5lu") %inode %pNode->getID() << endl;
				};
				// Edge
				size_t nNEdge=pLgrp->getNumOfEdge();
				for(size_t iedge=0; iedge < nNEdge; iedge++){
					ofs << format(" %s %5lu %5lu %5lu %5lu") %pLgrp->getEdgeType(iedge) %iedge %pLgrp->getElemID(iedge) %pLgrp->getConvMW3EdgeN(iedge) %0;

					vector<size_t> vBNode = pLgrp->getBNode_Edge(iedge);
					for(size_t ibnode=0; ibnode < vBNode.size(); ibnode++) ofs << format(" %5lu") %vBNode[ibnode];
					ofs << format("   %10.3e") %0.000 << endl;//境界値
				};
				ofs << "End\n" << endl;
			}else{
				/*
				タグ：BoundaryEdge
				書式 ( ディレクレ )
				辺境界メッシュID　境界条件種類  メッシュパーツID  境界節点数　Edge数
				境界節点ID  節点ID  DOF数　 DOF番号 境界値 DOF番号 境界値 DOF番号 境界値…
				…
				辺形状タイプ 境界辺ID  要素ID  辺番号  境界節点ID 境界節点ID 
				…
				*/
				ofs << "BoundaryEdge" << endl;
				ofs << format(" %lu %s %5lu %5lu %5lu") %pLgrp->getID() %msDirichlet %pMesh->getID() %pLgrp->getNumOfNode() %pLgrp->getNumOfEdge() << endl;
				// BNode 
				size_t nNNode=pLgrp->getNumOfNode();
				for(size_t inode=0; inode < nNNode; inode++){
					CNode* pNode = pLgrp->getNode(inode);
					ofs << format(" %5lu %5lu %lu %lu %lf") %inode %pNode->getID() %1 %0 %0.000 << endl;
				};
				// Edge
				size_t nNEdge=pLgrp->getNumOfEdge();
				for(size_t iedge=0; iedge < nNEdge; iedge++){
					ofs << format(" %s %5lu %5lu %5lu") %pLgrp->getEdgeType(iedge) %iedge %pLgrp->getElemID(iedge) %pLgrp->getConvMW3EdgeN(iedge);

					vector<size_t> vBNode = pLgrp->getBNode_Edge(iedge);
					for(size_t ibnode=0; ibnode < vBNode.size(); ibnode++) ofs << format(" %5lu") %vBNode[ibnode];
					ofs << endl;
				};
				ofs << "End\n" << endl;
			}//if(BndType)
		};
	};

	// メッセージ表示
	string sSumGrpN= lexical_cast<string>(nSumGrpN);
	pMsg->info("wrote :"+ sSumGrpN +" block \"BoundaryEdge\"");

	return true;
}
//--
// BoundaryVolumeMesh : デフォルトはNeumann型で出力(理由:FrontISTRは境界値を使用しない)
//--
bool CFileWriter::WriteBndVolMesh(CAssyModel *pAssyModel, ofstream& ofs)
{
	if(!pAssyModel->existMesh()) return false;//Meshが無い

	CMessage* pMsg=CMessage::Instance();

	size_t nNMesh=pAssyModel->getNumOfMesh();
	size_t nSumGrpN=0;
	string sBndType;
	/*
	タグ：BoundaryVolumeMesh
	書式
	メッシュパーツID  境界条件数
	体積境界メッシュID  境界条件種類　名称  自由度数 自由度番号 自由度番号 ...(繰り返し) 数式数(DOF数) DOF 数式 DOF 数式 ...
	…
	*/
	for(size_t i=0; i < nNMesh; i++){
		CMesh *pMesh=pAssyModel->getMesh(i);
		nSumGrpN += pMesh->getNumOfEgrp();

		if(pMesh->getNumOfEgrp() > 0){
			ofs << "BoundaryVolumeMesh" << endl;
			ofs << format(" %lu %lu") %pMesh->getID() %pMesh->getNumOfEgrp() << endl;
			for(size_t ii=0; ii < pMesh->getNumOfEgrp(); ii++){
				CEgroup *pEgrp=pMesh->getEgrp(ii);

				switch(pEgrp->getBndType()){ //境界種類あり？
				case( BndType::Dirichlet ): sBndType = msDirichlet; break;
				case( BndType::Neumann ): sBndType = msNeumann; break;
				case( BndType::NotUse ): sBndType = msNeumann; break;
				default: 
					sBndType = msNeumann;
					break;
				}

				ofs << format(" %lu %s %s %lu %lu %lu") %pEgrp->getID() %sBndType %pEgrp->getGroupName() %1 %0 %0 << endl;//自由度数:1,自由度:0,数式数:0
			};
			ofs << "End\n" << endl;
		}
	};

	// regex
	string strRegex= FistrElemTypeS::regStandard();// MW3に存在する要素タイプ
	sregex regBasis = sregex::compile(strRegex);
	string strRegex2= FistrElemTypeS::regShellType();// シェル要素
	sregex regShell= sregex::compile(strRegex2);
	smatch res;

	for(size_t i=0; i < nNMesh; i++){
		CMesh *pMesh=pAssyModel->getMesh(i);
		size_t nNEgrp= pMesh->getNumOfEgrp();
		for(size_t ii=0; ii < nNEgrp; ii++){
			CEgroup *pEgrp=pMesh->getEgrp(ii);

			if(pEgrp->getBndType()!=BndType::Dirichlet){
				/*
				タグ：BoundaryVolume
				書式 ( ノイマン )
				体積境界メッシュID　境界条件種類  メッシュパーツID  境界節点数 要素数
				境界節点ID  節点ID  
				…
				形状タイプ   境界体積ID  要素ID  0  自由度番号 境界節点ID...(繰り返し)   境界値
				…
				*/
				ofs << "BoundaryVolume" << endl;
				ofs << format(" %lu %s %5lu %5lu %5lu") %pEgrp->getID() %msNeumann %pMesh->getID() %pEgrp->getNumOfNode() %pEgrp->getNumOfElement() << endl;
				// Node
				for(size_t inode=0; inode < pEgrp->getNumOfNode(); inode++){
					ofs << format(" %5lu %5lu") %inode %pEgrp->getNodeID(inode) << endl;
				};
				// Volume(Element)
				for(size_t ielem=0; ielem < pEgrp->getNumOfElement(); ielem++){

					CElement *pElem= pEgrp->getElement(ielem);
					string sType= pElem->getCType();
				
					if( regex_search(sType,res,regBasis) ){
						// 要素形状がMW３に存在.
						ofs << format(" %s %5lu %5lu %5lu %5lu") %sType %ielem %pElem->getID() %0 %0 ;
					}else if( regex_search(sType,res,regShell) ){
						// シェル要素⇒置き換え
						if(sType==FistrElemTypeS::TriShell()) sType=FistrElemTypeS::Triangle();
						if(sType==FistrElemTypeS::TriShell2())sType=FistrElemTypeS::Triangle2();
						if(sType==FistrElemTypeS::QuadShell())sType=FistrElemTypeS::Quad();
						if(sType==FistrElemTypeS::QuadShell2())sType=FistrElemTypeS::Quad2();
						ofs << format(" %s %5lu %5lu %5lu %5lu") %sType %ielem %pElem->getID() %0 %0 ;
					}else{
						// 規格外
						pMsg->error( res.str()+" not supported");//Error --- IFaceQuadも含む.
					}
					// 要素構成BNode
					for(size_t enode=0; enode < pElem->getNumOfNode(); enode++){
						size_t nID= pElem->getNodeID_Fistr2MW3(enode);
						size_t nIX= pEgrp->getBNodeIndex(nID);
						ofs << format(" %5lu") %nIX;
					};
					ofs << format("   %10.3e") %0.000 << endl;//境界値
				};
				ofs << "End\n" << endl;
			}else{
				/*
				タグ：BoundaryVolume
				書式 ( ディレクレ )
				体積境界メッシュID　境界条件種類  メッシュパーツID  境界節点数 要素数
				境界節点ID  節点ID  DOF数　 DOF番号 境界値 DOF番号 境界値 DOF番号 境界値…
				…
				形状タイプ   境界体積ID  要素ID  0  　境界節点ID...(繰り返し)   
				…
				*/
				ofs << "BoundaryVolume" << endl;
				ofs << format(" %lu %s %5lu %5lu %5lu") %pEgrp->getID() %msDirichlet %pMesh->getID() %pEgrp->getNumOfNode() %pEgrp->getNumOfElement() << endl;
				// Node
				for(size_t inode=0; inode < pEgrp->getNumOfNode(); inode++){
					ofs << format(" %5lu %5lu %lu %lu %lf") %inode %pEgrp->getNodeID(inode) %1 %0 %0.000 << endl;
				};
				// Volume(Element)
				for(size_t ielem=0; ielem < pEgrp->getNumOfElement(); ielem++){

					CElement *pElem= pEgrp->getElement(ielem);
					string sType= pElem->getCType();
				
					if( regex_search(sType,res,regBasis) ){
						// 要素形状がMW３に存在.
						ofs << format(" %s %5lu %5lu %5lu") %sType %ielem %pElem->getID() %0 ;
					}else if( regex_search(sType,res,regShell) ){
						// シェル要素⇒置き換え
						if(sType==FistrElemTypeS::TriShell()) sType=FistrElemTypeS::Triangle();
						if(sType==FistrElemTypeS::TriShell2())sType=FistrElemTypeS::Triangle2();
						if(sType==FistrElemTypeS::QuadShell())sType=FistrElemTypeS::Quad();
						if(sType==FistrElemTypeS::QuadShell2())sType=FistrElemTypeS::Quad2();
						ofs << format(" %s %5lu %5lu %5lu") %sType %ielem %pElem->getID() %0;
					}else{
						// 規格外
						pMsg->error( res.str()+" not supported");//Error --- IFaceQuadも含む.
					}
					// 要素構成BNode
					for(size_t enode=0; enode < pElem->getNumOfNode(); enode++){
						size_t nID= pElem->getNodeID_Fistr2MW3(enode);
						size_t nIX= pEgrp->getBNodeIndex(nID);
						ofs << format(" %5lu") %nIX;
					};
					ofs << endl;
				};
				ofs << "End\n" << endl;
			}//if(BndType)
		};
	};


	// メッセージ表示
	string sSumGrpN= lexical_cast<string>(nSumGrpN);
	pMsg->info("wrote :"+ sSumGrpN +" block \"BoundaryVolume\"");

	return true;
}
//--
// ElementGroup
//   ⇒ BoundaryVolumeMeshデータをElementGroupとして読み替えて出力
//--
bool CFileWriter::WriteElemGroup(CAssyModel *pAssyModel, ofstream& ofs)
{
	if(!pAssyModel->existMesh()) return false;//Meshが無い

	CMessage* pMsg=CMessage::Instance();

	size_t nNMesh=pAssyModel->getNumOfMesh();
	size_t nSumGrpN=0;

	/*
	ElementGroup
		グループ番号  グループ名  メッシュパーツ番号
		…
	End
	*/
	for(size_t i=0; i < nNMesh; i++){
		CMesh *pMesh=pAssyModel->getMesh(i);
		nSumGrpN += pMesh->getNumOfEgrp();
	};
	if(nSumGrpN > 0){
		ofs << "ElementGroup" << endl;
		for(size_t i=0; i < nNMesh; i++){
			CMesh *pMesh=pAssyModel->getMesh(i);
			for(size_t ii=0; ii < pMesh->getNumOfEgrp(); ii++){
				CEgroup *pEgrp=pMesh->getEgrp(ii);
				ofs << format(" %lu %s %lu") %pEgrp->getID() %pEgrp->getGroupName() %pMesh->getID() << endl;
			};//---for(グループ数)
		};//---for(メッシュ数)end
		ofs << "End\n" << endl;
	}

	/*
	ElementGroupEntity
		グループ番号  メッシュパーツ番号
		要素ID  要素ID  要素ID  要素ID  要素ID
		要素ID  要素ID  要素ID  要素ID  要素ID
		…
	End
	*/
	if(nSumGrpN > 0){
		for(size_t i=0; i < nNMesh; i++){
			CMesh *pMesh=pAssyModel->getMesh(i);
			for(size_t ii=0; ii < pMesh->getNumOfEgrp(); ii++){
				CEgroup *pEgrp=pMesh->getEgrp(ii);
				ofs << "ElementGroupEntity" << endl;
				ofs << format(" %lu %lu")  %pEgrp->getID() %pMesh->getID() << endl;

				for(size_t ielem=0; ielem < pEgrp->getNumOfElement(); ielem++){
					CElement *pElem= pEgrp->getElement(ielem);
					
					ofs << format(" %lu")  %pElem->getID();
					if(ielem!=0 && ielem%5==0) ofs << endl;
				};
				if((pEgrp->getNumOfElement()-1)%5!=0) ofs << endl;
				ofs << "End\n" << endl;
			};
		};
	}


	// メッセージ表示
	string sSumGrpN= lexical_cast<string>(nSumGrpN);
	pMsg->info("wrote :"+ sSumGrpN +" block \"ElementGroup\"");

	return true;
}

//--
// ContactMesh
//--
bool CFileWriter::WriteContactMesh(CAssyModel *pAssyModel, ofstream& ofs)
{
	if(!pAssyModel->existConPair() && !pAssyModel->existMPCPair()) return false;//MPCとContactの両方とも無い

	CMessage *pMsg=CMessage::Instance();
	
	/*
	タグ：ContactMesh
	書式
	接合メッシュ数 最大ID  最小ID
	接合メッシュID  自身Rank 　属性番号(0:MPC, 1:Contact, ...)
	接合節点数 最大ID  最小ID
	接合節点ID  X  Y  Z  メッシュパーツID  節点ID  rank  maslave  v 3 0 ("v 3 0"はダミー定数)
	…
	マスター接合面数  最大ID  最小ID
	マスター面ID メッシュパーツID  要素ID  (要素)面番号  形状  接合節点ID …  自身Rank
	…
	スレーブ接合面数  最大ID  最小ID
	スレーブ面ID メッシュパーツID  要素ID  (要素)面番号  形状  接合節点ID …  自身Rank
	…
	
	※ maslave：0:マスター、1:スレーブ
	*/

	size_t nNumMPC= pAssyModel->getNumOfMPCPair();
	
	//MPC全体 
	ofs << "ContactMesh" << endl;
	ofs << format(" %lu %lu %lu") %nNumMPC %pAssyModel->getMaxIDinMPCPair() %pAssyModel->getMinIDinMPCPair() << endl;
	
	for(size_t i=0; i < nNumMPC; i++){
		// MPC接合面(単位-接合面)
		ofs << format(" %lu %lu %lu") %i %0 %0 << endl;
		
		CGroupPair *pMPCPair= pAssyModel->getMPCPair(i);
		
		size_t nAllNNum= pMPCPair->getNumOfAllNode();   //総節点数
		size_t nMNNum  = pMPCPair->getNumOfMasterNode();//マスター節点数
		size_t nMMeshID= pMPCPair->getMasterMeshID();   //マスターパーツID
		size_t nSMeshID= pMPCPair->getSlaveMeshID();    //スレーブパーツID
		// 接合節点
		ofs << format(" %9lu %9lu %9lu") %nAllNNum %pMPCPair->getMaxNodeIX() %pMPCPair->getMinNodeIX() << endl;
		for(size_t inode=0; inode < nAllNNum; inode++){

			CNode* pNode= pMPCPair->getNode(inode);
			double X=pNode->getX(), Y=pNode->getY(), Z=pNode->getZ();
			size_t nNID=pNode->getID();
			if(inode < nMNNum){
				//マスター点
				ofs << format(" %9lu %23.15e %23.15e %23.15e %lu %9lu %lu %lu %s") %inode %X %Y %Z %nMMeshID %nNID %0 %0 %"v 3 0" << endl;
			}else{
				//スレーブ点
				ofs << format(" %9lu %23.15e %23.15e %23.15e %lu %9lu %lu %lu %s") %inode %X %Y %Z %nSMeshID %nNID %0 %1 %"v 3 0" << endl;
			}
		};

		//マスター面
		size_t nMFaceNum= pMPCPair->getNumOfMasterFace();
		size_t nMaxMFaceIX = pMPCPair->getMaxMFaceIX();
		size_t nMinMFaceIX = pMPCPair->getMinMFaceIX();
		ofs << format(" %9lu %9lu %9lu") %nMFaceNum %nMaxMFaceIX %nMinMFaceIX << endl;
		// マスター面ID メッシュパーツID  要素ID  (要素)面番号  形状  接合節点ID …  自身Rank
		for(size_t iface=0; iface < nMFaceNum; iface++){
			size_t nElemID= pMPCPair->getMasterElemID(iface);
			size_t nFaceN = pMPCPair->getMasterFaceN_MW3(iface);
			string sFaceTye=pMPCPair->getMasterFaceType(iface);
			ofs << format(" %9lu %9lu %9lu %9lu %s") %iface %nMMeshID %nElemID %nFaceN %sFaceTye;

			vector<size_t> vConIX= pMPCPair->getMFaceNodeIX(iface);
			for(size_t inode=0; inode < vConIX.size(); inode++){
				ofs << format(" %9lu") %vConIX[inode];
			};
			ofs << format(" %lu") %0 << endl;
		};

		//スレーブ面
		size_t nSFaceNum= pMPCPair->getNumOfSlaveFace();
		size_t nMaxSFaceIX= pMPCPair->getMaxSFaceIX();
		size_t nMinSFaceIX= pMPCPair->getMinSFaceIX();
		ofs << format(" %9lu %9lu %9lu") %nSFaceNum %nMaxSFaceIX %nMinSFaceIX << endl;
		// スレーブ面ID メッシュパーツID  要素ID  (要素)面番号  形状  接合節点ID …  自身Rank
		for(size_t iface=0; iface < nSFaceNum; iface++){
			size_t nElemID= pMPCPair->getSlaveElemID(iface);
			size_t nFaceN = pMPCPair->getSlaveFaceN_MW3(iface);
			string sFaceTye=pMPCPair->getSlaveFaceType(iface);
			ofs << format(" %9lu %9lu %9lu %9lu %s") %iface %nSMeshID %nElemID %nFaceN %sFaceTye;

			vector<size_t> vConIX= pMPCPair->getSFaceNodeIX(iface);
			for(size_t inode=0; inode < vConIX.size(); inode++){
				ofs << format(" %9lu") %vConIX[inode];
			};
			ofs << format(" %lu") %0 << endl;
		};

	};
	ofs << "End\n" << endl;

	string sNumMPC= lexical_cast<string>(nNumMPC);
	pMsg->info("wrote :"+sNumMPC+" block \"ContactMesh(MPC)\"" );
	
	return true;
}


//--
// MW3 メッシュ出力
//--
void CFileWriter::WriteMesh_MW3(CAssyModel *pAssyModel, string filename)
{
	CMessage *pMsg=CMessage::Instance();

	string strMW3std  = FileName::regMW3Mesh();
	string strMW3fistr= FileName::regMW3FistrMesh();
	sregex regMW3std  = sregex::compile(strMW3std);
	sregex regMW3fistr= sregex::compile(strMW3fistr);
	smatch res1, res2;
	//--
	// ファイル名チェック : rank番号"0"の付け方
	//--
	// 1.MW3標準          メッシュファイル名:".0.msh"
	// 2.MW3FrontISTR対応 メッシュファイル名:".msh.0"
	//
	if(!regex_search(filename, res1, regMW3std) && !regex_search(filename, res2, regMW3fistr) ){
		// Warn
		pMsg->warn("invalid output file extension.");

		string user_str;
		pMsg->user_in("select file extension kind: enter STD or FISTR => ", user_str);//----------------------- ユーザー選択
		transform(user_str.begin(), user_str.end(), user_str.begin(), ::tolower);//-- 全て小文字に変換 

		string strExt = FileName::regExt();
		sregex regExt = sregex::compile(strExt);
		smatch resExt;
		
		if( regex_search(filename, resExt, regExt) ){
			//--
			// 拡張子".msh"が付いているがrank番号は無い場合
			//--
			if(user_str=="std"){//MW3標準
				filename = regex_replace(filename, regExt, string(".0.msh") );//".msh"⇒".0.msh"に置換
			}else if(user_str=="fistr"){//MW3FrontISTR対応
				filename = regex_replace(filename, regExt, string(".msh.0") );//".msh"⇒".msh.0"に置換
			}else{
				pMsg->error("input error, change to std type.");
				filename = regex_replace(filename, regExt, string(".0.msh") );//".msh"⇒".0.msh"に置換
			}
		}else{
			//--
			// 拡張子に.mshすらない場合
			//--
			if(user_str=="std"){   filename += ".0.msh";}else
			if(user_str=="fistr"){ filename += ".msh.0";}
			else{ pMsg->error("input error, change to std type."); filename += ".0.msh";}
		}
	}//if(拡張子にrank番号の有無) 


	ofstream ofs;
	ofs.open(filename.c_str());// Open
	//--
	// ファイル出力
	//--
	if(ofs){
		pMsg->info("file open:\""+filename+"\"");

		WriteHeader(ofs);
		WriteAssyModel(pAssyModel, ofs);//---AssyModel
		WriteNode(pAssyModel, ofs);     //---Node
		WriteElement(pAssyModel, ofs);  //---Element
		WriteBndNodeMesh(pAssyModel, ofs);//-BoundaryNodeMesh
		WriteBndFaceMesh(pAssyModel, ofs);//-BoundaryFaceMesh
		WriteBndEdgeMesh(pAssyModel, ofs);//-BoundaryEdgeMesh
		////WriteBndVolMesh(pAssyModel, ofs); //-BoundaryVolumeMesh ← fistr v4 は不使用
		WriteElemGroup(pAssyModel, ofs);//---ElementGroup ← BoundaryVolumeMeshの読み替え
		//
		// # 通信界面はシングルなので無い : CommMesh2
		//
		WriteContactMesh(pAssyModel, ofs);//-ContactMesh

		pMsg->info("file close:\""+filename+"\"\n");
	}else{
		pMsg->error("\""+filename+"\" could not open.");
		exit(0);//------------------------------------------ exit(0)
	}

	ofs.close();// Close
}



