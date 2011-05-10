#ifndef HECMW_VARRAY_INT_INCLUDED
#define HECMW_VARRAY_INT_INCLUDED

struct hecmw_varray_int {
	int n_val;
	int max_val;

	int *vals;
};


extern int HECMW_varray_int_init(struct hecmw_varray_int *varray);

extern void HECMW_varray_int_finalize(struct hecmw_varray_int *varray);


extern int HECMW_varray_int_nval(const struct hecmw_varray_int *varray);

extern int HECMW_varray_int_append(struct hecmw_varray_int *varray, int value);

extern int HECMW_varray_int_get(const struct hecmw_varray_int *varray, int index);

extern int HECMW_varray_int_cat(struct hecmw_varray_int *varray,
				const struct hecmw_varray_int *varray2);

extern int HECMW_varray_int_sort(struct hecmw_varray_int *varray);

extern int HECMW_varray_int_uniq(struct hecmw_varray_int *varray);


extern int HECMW_varray_int_resize(struct hecmw_varray_int *varray, int len);

extern int *HECMW_varray_int_get_v(struct hecmw_varray_int *varray);


extern int HECMW_varray_int_copy(const struct hecmw_varray_int *varray,
				 struct hecmw_varray_int *varray2);

extern int HECMW_varray_int_rmdup(struct hecmw_varray_int *varray);

#endif /* HECMW_VARRAY_INT_INCLUDED */
