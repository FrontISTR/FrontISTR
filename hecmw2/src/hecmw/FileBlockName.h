/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileBlockName.h
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
#ifndef FILE_BLOCK_NAME_7B85EDDD_CCF3_4c0b_B435_154AA1897A9C
#define FILE_BLOCK_NAME_7B85EDDD_CCF3_4c0b_B435_154AA1897A9C
#include "TypeDef.h"
#include <cstring>
using namespace std;
namespace FileIO
{
struct FileBlockName {
    static const char* StartMeshFileName() {
        return "MeshFileName";
    }
    static uiint MeshFileName_Len() {
        return strlen(StartMeshFileName());
    }
    static const char* EndMeshFileName() {
        return "End";
    }

    //--
    // メッシュ・バージョン
    //--
    static const char* StartVersion() {
        return "Version";
    }
    static const char* Version_Num() {
        return "0.4";   //----------mesh version 0.4
    }
    static const char* Default_Version_Num() {
        return "0.3";   //--バージョンブロックが存在しない場合のデフォルト:0.3
    }

    //--
    // メッシュデータ ブロック名
    //--
    static const char* StartNode() {
        return "Node";
    }
    static uiint Node_Len() {
        return strlen(StartNode());
    }
    static const char* EndNode() {
        return "End";
    }

    static const char* StartElement() {
        return "Element";
    }
    static uiint Element_Len() {
        return strlen(StartElement());
    }
    static const char* EndElement() {
        return "End";
    }

    static const char* StartBoundaryNodeMesh() {
        return "BoundaryNodeMesh";
    }
    static uiint BoundaryNodeMesh_Len() {
        return strlen(StartBoundaryNodeMesh());
    }
    static const char* EndBoundaryNodeMesh() {
        return "End";
    }

    static const char* StartBoundaryFaceMesh() {
        return "BoundaryFaceMesh";
    }
    static uiint BoundaryFaceMesh_Len() {
        return strlen(StartBoundaryFaceMesh());
    }
    static const char* EndBoundaryFaceMesh() {
        return "End";
    }

    static const char* StartBoundaryVolumeMesh() {
        return "BoundaryVolumeMesh";
    }
    static uiint BoundaryVolumeMesh_Len() {
        return strlen(StartBoundaryVolumeMesh());
    }
    static const char* EndBoundaryVolumeMesh() {
        return "End";
    }

    static const char* StartBoundaryEdgeMesh() {
        return "BoundaryEdgeMesh";
    }
    static uiint BoundaryEdgeMesh_Len() {
        return strlen(StartBoundaryEdgeMesh());
    }
    static const char* EndBoundaryEdgeMesh() {
        return "End";
    }

    static const char* StartBoundaryNode() {
        return "BoundaryNode";
    }
    static uiint BoundaryNode_Len() {
        return strlen(StartBoundaryNode());
    }
    static const char* EndBoundaryNode() {
        return "End";
    }

    static const char* StartBoundaryFace() {
        return "BoundaryFace";
    }
    static uiint BoundaryFace_Len() {
        return strlen(StartBoundaryFace());
    }
    static const char* EndBoundaryFace() {
        return "End";
    }

    static const char* StartBoundaryEdge() {
        return "BoundaryEdge";
    }
    static uiint BoundaryEdge_Len() {
        return strlen(StartBoundaryEdge());
    }
    static const char* EndBoundaryEdge() {
        return "End";
    }

    static const char* StartBoundaryVolume() {
        return "BoundaryVolume";
    }
    static uiint BoundaryVolume_Len() {
        return strlen(StartBoundaryVolume());
    }
    static const char* EndBoundaryVolume() {
        return "End";
    }

    static const char* StartMaterial() {
        return "Material";
    }
    static uiint Material_Len() {
        return strlen(StartMaterial());
    }
    static const char* EndMaterial() {
        return "End";
    }

    static const char* StartElementGroup() {
        return "ElementGroup";
    }
    static uiint ElementGroup_Len() {
        return strlen(StartElementGroup());
    }
    static const char* EndElementGroup() {
        return "End";
    }

    static const char* StartElementGroupEntity() {
        return "ElementGroupEntity";
    }
    static uiint ElementGroupEntity_Len() {
        return strlen(StartElementGroupEntity());
    }
    static const char* EndElementGroupEntity() {
        return "End";
    }

    static const char* StartCommMesh() {
        return "CommMesh";
    }
    static uiint CommMesh_Len() {
        return strlen(StartCommMesh());
    }
    static const char* EndCommMesh() {
        return "End";
    }

    static const char* StartCommNode() {
        return "CommNode";
    }
    static uiint CommNode_Len() {
        return strlen(StartCommNode());
    }
    static const char* EndCommNode() {
        return "End";
    }

    static const char* StartCommElement() {
        return "CommElement";
    }
    static uiint CommElement_Len() {
        return strlen(StartCommElement());
    }
    static const char* EndCommElement() {
        return "End";
    }

    static const char* StartCommMesh2() {
        return "CommMesh2";
    }
    static uiint CommMesh2_Len() {
        return strlen(StartCommMesh2());
    }
    static const char* EndCommMesh2() {
        return  "End";
    }

    static const char* StartCommNodeCM2() {
        return "CommNodeCM2";
    }
    static uiint CommNodeCM2_Len() {
        return strlen(StartCommNodeCM2());
    }
    static const char* EndCommNodeCM2() {
        return "End";
    }

    static const char* StartCommFace() {
        return "CommFace";
    }
    static uiint CommFace_Len() {
        return strlen(StartCommFace());
    }
    static const char* EndCommFace() {
        return "End";
    }

    static const char* StartRefine() {
        return "Refine";
    }
    static uiint Refine_Len() {
        return strlen(StartRefine());
    }
    static const char* EndRefine() {
        return "End";
    }

    static const char* StartAssyModel() {
        return "AssyModel";
    }
    static uiint AssyModel_Len() {
        return strlen(StartAssyModel());
    }
    static const char* EndAssyModel() {
        return "End";
    }

    static const char* StartMesh() {
        return "Mesh";
    }
    static uiint Mesh_Len() {
        return strlen(StartMesh());
    }
    static const char* EndMesh() {
        return "End";
    }

    static const char* StartContactMesh() {
        return "ContactMesh";
    }
    static uiint ContactMesh_Len() {
        return strlen(StartContactMesh());
    }
    static const char* EndContactMesh() {
        return "End";
    }


    //--
    // リスタートデータ ブロック名
    //--
    static const char* StartRes() {
        return "Restart";
    }
    static uiint Res_Len() {
        return strlen(StartRes());
    }
    static const char* EndRes() {
        return "End";
    }

    static const char* StartAlgebra() {
        return "Algebra";
    }
    static uiint Algebra_Len() {
        return strlen(StartAlgebra());
    }
    static const char* EndAlgebra() {
        return "End";
    }

    //--
    // リザルト ブロック名
    //--
    static const char* StartResult() {
        return "Result";
    }
    static uiint Result_Len() {
        return strlen(StartResult());
    }
    static const char* EndResult() {
        return "End";
    }

    //--
    // バイナリ・チェック ブロック名
    //--
    static const char* StartBinCheck() {
        return "BinCheck";
    }
    static uiint BinCheck_Len() {
        return strlen(StartBinCheck());
    }
    static const char* EndBinCheck() {
        return "End";
    }

    static const char* StartTitle() {
        return "Title";
    }
    static uiint Title_Len() {
        return strlen(StartTitle());
    }
    static const char* EndTitle() {
        return "End";
    }


    //--
    // End タグ(共通)
    //--
    static const char* End() {
        return "End";
    }
    static uiint End_Len() {
        return strlen(End());
    }
};
}
#endif
