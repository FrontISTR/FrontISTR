/* 
 * File:   FileWriterAVS.h
 * Author: ktakeda
 *
 * Created on 2011/04/09, 13:24
 */
#include "FileWriter.h"
#include "AssyVector.h"//解ベクトル

namespace FileIO{
#ifndef FILEWRITERAVS_H
#define	FILEWRITERAVS_H
class CFileWriterAVS:public CFileWriter{
private:
    CFileWriterAVS();
public:
    static CFileWriterAVS* Instance(){
        static CFileWriterAVS moWriterAVS;
        return &moWriterAVS;
    }
    virtual ~CFileWriterAVS();

private:
    // --
    // Meshパーツ別
    // --
    vector<vstring> mvLabel;//ラベル名の集合
    vector<map<string, string> >  mmUnit;//ラベル毎の単位名
    vector<map<string, uiint> >   mmDOF; //ラベル毎のDOF
    vector<map<string, vdouble> > mmVecVal;//ラベル毎の値の配列(節点×自由度)

    void WriteMesh(ofstream& ofs, const uiint& iMesh, const uiint& mgLevel);
public:
    virtual void WriteDebug(ofstream& ofs, const uiint& mgLevel);

    //変数登録
    void recAVS_Label(const uiint& iMesh, char* cLabel, char* cUnit, const uiint& nNumOfDOF);
    void recAVS_Variable(const uiint& iMesh, const uiint& nNumOfNode, char* cLabel, double* pvValue);

    //出力
    void WriteBasis(ofstream& ofs, const uiint& iMesh, const uiint& mgLevel);//基礎変数の出力
    void WriteFEM(ofstream& ofs, const uiint& iMesh, const uiint& mgLevel);  //登録変数の出力
};
#endif	/* FILEWRITERAVS_H */
}

