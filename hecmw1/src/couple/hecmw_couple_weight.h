/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 1.00                                              *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : Coupling Interface                                *
 *                                                                     *
 *            Written by Shin'ichi Ezure (RIST)                        *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using Hight End Computing Middleware (HEC-MW)"     *
 *                                                                     *
 *=====================================================================*/




#ifndef INC_HECMW_COUPLE_WEIGHT
#define INC_HECMW_COUPLE_WEIGHT


struct hecmw_couple_weight {
	int n;				/* アイテム数 (被補間節点)					*/
	int type;			/* アイテムタイプ (補間アイテム)			*/
	int *index;			/* インデックス配列 (補間アイテム)			*/
	int *id;			/* 連成境界情報におけるID (補間アイテム)	*/
	double *weight;		/* 重み係数 (補間アイテム)					*/
};


struct hecmw_couple_weight_list {
	struct hecmw_couple_weight *info;		/* 補間係数情報						*/
	struct hecmw_couple_weight_list *next;	/* 次の補間係数リストへのポインタ	*/
};


extern struct hecmw_couple_weight *
HECMW_couple_alloc_weight(void);
extern void
HECMW_couple_free_weight(struct hecmw_couple_weight *p);

extern struct hecmw_couple_weight_list *
HECMW_couple_alloc_weight_list(void);
extern void
HECMW_couple_free_weight_list(struct hecmw_couple_weight_list *p);

#endif	/* INC_HECMW_COUPLE_WEIGHT */
