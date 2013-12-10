/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   FrontISTR_Type.h
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef E8F2AD0_FFFD_4ce7_B818_51DA3511C5DF
#define E8F2AD0_FFFD_4ce7_B818_51DA3511C5DF

//--
// 要素タイプ
//--
struct FistrElementType{
  enum{
		Beam=111,
    Beam2=112,
		Triangle=231,
    Triangle2=232,
    Quad=241,
    Quad2=242,
    Tetra=341,
    Tetra2=342,
    Prism=351,
    Prism2=352,
    Hexa=361,
    Hexa2=362,
		IFaceQuad=541,
		IFaceQuad2=542,
		TriShell=731,
    TriShell2=732,
    QuadShell=741,
		QuadShell2=742
  };
};

#include <string>
//--
// 要素タイプ(string)
//--
struct FistrElemTypeS{

	// MW3 に存在する要素タイプ
	static string Beam(){ return "Beam";}
	static string Beam2(){ return "Beam2";}
	static string Triangle(){ return "Triangle";}
	static string Triangle2(){ return "Triangle2";}
	static string Quad(){ return "Quad";}
	static string Quad2(){ return "Quad2";}
	static string Tetra(){ return "Tetra";}
	static string Tetra2(){ return "Tetra2";}
	static string Prism(){ return "Prism";}
	static string Prism2(){ return "Prism2";}
	static string Hexa(){ return "Hexa";}
	static string Hexa2(){ return "Hexa2";}
	// MW3 に存在する要素タイプ regex
	static string regStandard(){ return "Beam|Beam2|Triangle|Triangle2|Quad|Quad2|Tetra|Tetra2|Prism|Prism2|Hexa|Hexa2";}

	// MW3 に存在しない要素タイプ
	static string IFaceQuad(){ return "IFaceQuad";}
	static string IFaceQuad2(){ return "IFaceQuad2";}
	static string TriShell(){ return "TriShell";}
	static string TriShell2(){ return "TriShell2";}
	static string QuadShell(){ return "QuadShell";}
	static string QuadShell2(){ return "QuadShell2";}
	// MW3 に存在しない要素タイプ regex
	static string regNonStandard(){ return "IFaceQuad|IFaceQuad2|TriShell|TriShell2|QuadShell|QuadShell2";}

	// シェル要素タイプ
	static string regShellType(){ return "TriShell|TriShell2|QuadShell|QuadShell2";}
	// 全く対応不可なタイプ
	static string regNonSupported(){ return "IFaceQuad|IFaceQuad2";}

	// 平面・シェル
	static string regPlate(){ return "TriShell|TriShell2|QuadShell|QuadShell2|Triangle|Triangle2|Quad|Quad2";}
	// 線
	static string regLine(){ return "Beam|Beam2";}

	static string Unknown(){ return "Unknown";}
};

//--
// 出力ファイル名チェック
//--
struct FileName{
	// MW3 標準メッシュ・ファイルの後半 ".0.msh"であることを確認.
	static string regMW3Mesh(){ return "\\.0\\.msh$";}
	// MW3 FrontISTR対応メッシュ・ファイルの後半 ".msh.0"であることを確認.
	static string regMW3FistrMesh(){ return "\\.msh\\.0$";}

	// 拡張子 ".msh" の確認.
	static string regExt(){ return "\\.msh$";}
};
//--
// MW3 境界条件種類
//--
struct BndType{
	enum{
		Dirichlet,
		Neumann,
		Other,
		NotUse
	};
};

//--
// ファイル書式のタグ : FrontISTR ver4.1
//--
struct FistrTag{
	
	// ヘッダー 
	static string Header(){ return "!HEADER";}
	static string Node(){ return "!NODE";}
	static string Element(){ return "!ELEMENT";}
	static string Egroup(){ return "!EGROUP";}
	static string Sgroup(){ return "!SGROUP";}
	static string Ngroup(){ return "!NGROUP";}
	static string AssemblyPair(){ return "!ASSEMBLY_PAIR";}
	static string ContactPair(){ return "!CONTACT_PAIR";}
	static string End(){ return "!END";}
	static string Ex(){ return "!";} 

	// パラメータ 
	static string Ver(){ return "VER";}
	static string PartsName(){ return "PARTNAME";}
	static string Num(){ return "NUM";}
	static string Type(){ return "TYPE";}
	static string Egrp(){ return "EGRP";}
	static string Sgrp(){ return "SGRP";}
	static string Ngrp(){ return "NGRP";}
	static string Name(){ return "NAME";}
	static string Equal(){ return "=";}
	static string Generate(){ return "GENERATE";}



	// ヘッダー regex
	static string regHeader(){ return "^!HEADER";}
	static string regNode(){ return "^!NODE";}
	static string regElement(){ return "^!ELEMENT";}
	static string regEgroup(){ return "^!EGROUP";}
	static string regSgroup(){ return "^!SGROUP";}
	static string regNgroup(){ return "^!NGROUP";}
	static string regAssemblyPair(){ return "^!ASSEMBLY_PAIR";}
	static string regContactPair(){ return "^!CONTACT_PAIR";}
	static string regEnd(){ return "^!END";}

	static string regEx(){ return "^!.*";}
	static string regAllHeaders(){ return "!HEADER|!NODE|!ELEMENT|!EGROUP|!SGROUP|!NGROUP|!ASSEMBLY_PAIR|!CONTACT_PAIR|!END";}


	// パラメータ regex
	static string regVer(){ return "VER.*";}
	static string regPartsName(){ return "PARTNAME.*";}
	static string regNum(){ return "NUM.*";}
	static string regType(){ return "TYPE.*";}
	static string regEgrp(){ return "EGRP.*";}
	static string regSgrp(){ return "SGRP.*";}
	static string regNgrp(){ return "NGRP.*";}
	static string regName(){ return "NAME.*";}
	static string regEqual(){ return "=";}

	// バージョン番号(FrontISTRメッシュ書式バージョン)
	static string VerNumber(){ return "4";}

};

struct MW3_Tag{
	// バージョン番号(MW3 ファイル書式)
	static string VerNumber(){ return "0.4";}
};

#endif //include_guard





