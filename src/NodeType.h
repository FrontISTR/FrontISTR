/* 
 * File:   NodeType.h
 * Author: ktakeda
 *
 * Created on 2009/05/26, 16:52
 */

#ifndef _NODETYPE_H_0ca01a33_f12a_4ae4_9dde_3da45f62b934
#define	_NODETYPE_H_0ca01a33_f12a_4ae4_9dde_3da45f62b934

namespace pmw{
union NodeType{
    enum{
        Scalar,
        Vector,
        ScalarVector
    };
};
}
#endif	/* _NODETYPE_H */

