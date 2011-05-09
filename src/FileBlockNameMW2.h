//
//  old HEC_MW file_tag union
//  
//
//				2009.4.07
//				2009.3.24
//				k.Takeda
//
#ifndef FILE_TAG_NAME_CC0B8D86_67BD_493a_B552_FBF00C1F9B02
#define FILE_TAG_NAME_CC0B8D86_67BD_493a_B552_FBF00C1F9B02

namespace FileIO{
struct FileBlockNameMW2{
    // HEC_MW
    static const char* Comment(){ return "!!";}
    static const char* SingleExclamation(){ return "!";}// Exclamation(!)

    // HEC_MW mesh file
    static const char* Header(){ return "!HEADER";}
    static const char* Zero(){ return "!ZERO";}
    static const char* Node(){ return "!NODE";}
    static const char* Element(){ return "!ELEMENT";}
    static const char* NGroup(){ return "!NGROUP";}
    static const char* EGroup(){ return "!EGROUP";}
    static const char* SGroup(){ return "!SGROUP";}
    static const char* Equation(){ return "!EQUATION";}
    static const char* Amplitude(){ return "!AMPLITUDE";}
    static const char* Section(){ return "!SECTION";}
    static const char* Material(){ return "!MATERIAL";}
    static const char* Initial(){ return "!INITIAL CONDITION";}
    static const char* Include(){ return "!INCLUDE";}
    static const char* Connectivity(){ return "!CONNECTIVITY";}
    static const char* End(){ return "!END";}

    // HEC_MW cnt file
    static const char* Control(){return "!CONTROL";}
    static const char* Mesh(){ return "!MESH";}	
    static const char* MeshGroup(){ return "!MESH GROUP";}
    static const char* Restart(){return "!RESTART";}	  
    static const char* Result(){return "!RESULT";}
};
}
#endif
