/* 
 * File:   Film.h
 * Author: ktakeda
 *
 * 熱伝達率を表現する薄膜
 * 
 *
 * Created on 2009/05/19, 14:16
 */

#ifndef _FILM_H_ad201417_2f36_4da7_8847_40ffded43b0c
#define	_FILM_H_ad201417_2f36_4da7_8847_40ffded43b0c

#include "CommonStd.h"
#include "TypeDef.h"

#include "Element.h"

namespace pmw{
class CFilm:public CElement{
public:
    CFilm();
    virtual ~CFilm();

private:
    double mHeatTrans;//Heat Transfer Rate
    
    static uiint mnElemType;
    static uiint mnBaseType;

public:
    void setHeatTrans(const double& heat_trans);
    double& getHeatTrans();
};
}
#endif	/* _FILM_H */

