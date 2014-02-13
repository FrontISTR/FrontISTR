/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/IndexValuePair.h
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
#ifndef INDEXVALUEPAIR_H_
#define INDEXVALUEPAIR_H_
namespace pmw
{
template<int NDOF>
class CIndexValuePair
{
public:
    CIndexValuePair(int i, int j, double *val[NDOF])
        : m_i(i), m_j(j) {
        for (int k = 0; k < NDOF; k++)
            for (int l = 0; l < NDOF; l++)
                m_val[k][l] = val[k][l];
    }
    virtual ~CIndexValuePair() {}
    int operator <(CIndexValuePair& other) {
        if (m_i < other.m_i) return -1;
        if (m_i > other.m_i) return 1;
        if (m_j < other.m_j) return -1;
        if (m_j > other.m_j) return 1;
        return 0;
    }
public:
    int m_i;
    int m_j;
    double m_val[NDOF][NDOF];
};
}
#endif /* INDEXVALUEPAIR_H_ */
