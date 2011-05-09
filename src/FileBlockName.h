//
//  FileBlockName  struct
// 
//
//				2009.09.18
//				2008.12.10
//				k.Takeda
#ifndef FILE_BLOCK_NAME_7B85EDDD_CCF3_4c0b_B435_154AA1897A9C
#define FILE_BLOCK_NAME_7B85EDDD_CCF3_4c0b_B435_154AA1897A9C

/*
	<< HEC_MW2 Tag >>
	<<  FileTagName.h  >>

	<< cnt file  >>
	!CONTROL
	!MESH
	!MESH GROUP
	!RESTART
	!RESULT

	<< mesh file >>
	!HEADER
	!ZERO
	!NODE
	!ELEMENT
	!NGROUP
	!EGROUP
	!SGROUP
	!EQUATION
	!AMPLITUDE
	!SECTION
	!MATERIAL
	!INITIAL CONDITION
	!INCLUDE
	!CONNECTIVITY
	!END
*/

namespace FileIO{
struct FileBlockName{
    // Cnt
    // --
    static const char* StartMeshFileName(){ return "StartMeshFileName";}
    static const char* EndMeshFileName(){ return "EndMeshFileName";}


    // Mesh
    // ----
    // Node
    static const char* StartNode(){ return "StartNode";}
    static const char* EndNode(){ return "EndNode";}
    // Element
    static const char* StartElement(){ return "StartElement";}
    static const char* EndElement(){ return "EndElement";}
    // BoundaryNode
    static const char* StartBoundaryNode(){ return "StartBoundaryNode";}
    static const char* EndBoundaryNode(){ return "EndBoundaryNode";}
    // BoundaryFace
    static const char* StartBoundaryFace(){ return "StartBoundaryFace";}
    static const char* EndBoundaryFace(){ return "EndBoundaryFace";}
    // BoundaryVolume
    static const char* StartBoundaryVolume(){ return "StartBoundaryVolume";}
    static const char* EndBoundaryVolume(){ return "EndBoundaryVolume";}
    // Material
    static const char* StartMaterial(){ return "StartMaterial";}
    static const char* EndMaterial(){ return "EndMaterial";}





    // Communication
    // --
    // CommMesh
    static const char* StartCommMesh(){ return "StartCommMesh";}
    static const char* EndCommMesh(){ return "EndCommMesh";}
    // CommNode
    static const char* StartCommNode(){ return "StartCommNode";}
    static const char* EndCommNode(){ return "EndCommNode";}
    // CommElement
    static const char* StartCommElement(){ return "StartCommElement";}
    static const char* EndCommElement(){ return "EndCommElement";}




    // GMGModel
    //  MultiGrid Model
    // --
    static const char* StartRefine(){ return "StartRefine";}
    static const char* EndRefine(){ return "EndRefine";}
    // AssyModel(Assembly Model)
    static const char* StartAssyModel(){ return "StartAssyModel";}
    static const char* EndAssyModel(){ return "EndAssyModel";}
    // Mesh(Part)
    static const char* StartMesh(){ return "StartMesh";}
    static const char* EndMesh(){ return "EndMesh";}
    // ContactMesh(Interface)
    static const char* StartContactMesh(){ return "StartContactMesh";}
    static const char* EndContactMesh(){ return "EndContactMesh";}

};
}
#endif
