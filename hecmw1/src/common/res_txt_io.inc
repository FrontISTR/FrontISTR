/*---------------------------------------------------------------------------*/
/* TEXT MODE I/O ---- output_result                                          */
/*---------------------------------------------------------------------------*/


static int
output_result_header(FILE *fp)
{
	int rc;

	/* header */
	rc = fprintf(fp, "%s\n", head);
	if(rc < 0) {
		HECMW_set_error(HECMW_UTIL_E0205, "");
		return -1;
	}

	/* n_node, n_elem */
	rc = fprintf(fp, "%d %d\n", nnode, nelem);
	if(rc < 0) {
		HECMW_set_error(HECMW_UTIL_E0205, "");
		return -1;
	}

	/* nn_component, ne_component */
	rc = fprintf(fp, "%d %d\n", count_nn_comp(), count_ne_comp());
	if(rc < 0) {
		HECMW_set_error(HECMW_UTIL_E0205, "");
		return -1;
	}

	return 0;
}


static int
output_result_node(FILE *fp)
{
	int i,j,k,n,rc,nn_comp;
	struct result_list *p,**data;

	/* nn_dof */
	n = 0;
	for(p=node_list; p; p=p->next) {
		rc = fprintf(fp, "%d%c", p->n_dof, (n+1)%COL_INT ? ' ' : '\n');
		if(rc < 0) {
			HECMW_set_error(HECMW_UTIL_E0205, "");
			return -1;
		}
		n++;
	}
	if(n % COL_INT) {
		rc = fprintf(fp, "\n");
		if(rc < 0) {
			HECMW_set_error(HECMW_UTIL_E0205, "");
			return -1;
		}
	}

	/* node_label */
	for(p=node_list; p; p=p->next) {
		rc = fprintf(fp, "%s\n", p->label);
		if(rc < 0) {
			HECMW_set_error(HECMW_UTIL_E0205, "");
			return -1;
		}
	}

	/* node_val_item */
	nn_comp = count_nn_comp();
	if(nn_comp == 0) return 0;
	data = HECMW_malloc(sizeof(*data) * nn_comp);
	if(data == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	i = 0;
	for(p=node_list; p; p=p->next) {
		data[i++] = p;
	}
	for(i=0; i < nnode; i++) {
        rc = fprintf(fp, "%d \n", node_global_ID[i]);
		if(rc < 0) {
			HECMW_set_error(HECMW_UTIL_E0205, "");
			return -1;
		}
		n = 0;
		for(j=0; j < nn_comp; j++) {
			p = data[j];
			for(k=0; k < p->n_dof; k++) {
				rc = fprintf(fp, "%.16E%c", p->ptr[i*p->n_dof+k], (n+1)%COL_DOUBLE ? ' ' : '\n');
				n++;
			}
		}
		if(n % COL_DOUBLE) {
			rc = fprintf(fp, "\n");
			if(rc < 0) {
				HECMW_set_error(HECMW_UTIL_E0205, "");
				return -1;
			}
		}
	}
	HECMW_free(data);

	return 0;
}


static int
output_result_elem(FILE *fp)
{
	int i,j,k,n,rc,ne_comp;
	struct result_list *p,**data;

	/* ne_dof */
	n = 0;
	for(p=elem_list; p; p=p->next) {
		rc = fprintf(fp, "%d%c", p->n_dof, (n+1)%COL_INT ? ' ' : '\n');
		if(rc < 0) {
			HECMW_set_error(HECMW_UTIL_E0205, "");
			return -1;
		}
		n++;
	}
	if(n % COL_INT) {
		rc = fprintf(fp, "\n");
		if(rc < 0) {
			HECMW_set_error(HECMW_UTIL_E0205, "");
			return -1;
		}
	}

	/* elem_label */
	for(p=elem_list; p; p=p->next) {
		rc = fprintf(fp, "%s\n", p->label);
		if(rc < 0) {
			HECMW_set_error(HECMW_UTIL_E0205, "");
			return -1;
		}
	}

	/* elem_val_item */
	ne_comp = count_ne_comp();
	if(ne_comp == 0) return 0;
	data = HECMW_malloc(sizeof(*data) * ne_comp);
	if(data == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	i = 0;
	for(p=elem_list; p; p=p->next) {
		data[i++] = p;
	}
	for(i=0; i < nelem; i++) {
        rc = fprintf(fp, "%d\n", elem_global_ID[i]);
		if(rc < 0) {
			HECMW_set_error(HECMW_UTIL_E0205, "");
			return -1;
		}
		n = 0;
		for(j=0; j < ne_comp; j++) {
			p = data[j];
			for(k=0; k < p->n_dof; k++) {
				rc = fprintf(fp, "%.16E%c", p->ptr[i*p->n_dof+k], (n+1)%COL_DOUBLE ? ' ' : '\n');
				n++;
			}
		}
		if(n % COL_DOUBLE) {
			rc = fprintf(fp, "\n");
			if(rc < 0) {
				HECMW_set_error(HECMW_UTIL_E0205, "");
				return -1;
			}
		}
	}
	HECMW_free(data);

	return 0;
}


static int
output_result_data(FILE *fp)
{
	HECMW_assert(fp);

	if(output_result_header(fp)) {
		return -1;
	}
	if(output_result_node(fp)) {
		return -1;
	}
	if(output_result_elem(fp)) {
		return -1;
	}

	return 0;
}


/*---------------------------------------------------------------------------*/
/* TEXT MODE I/O ---- output_result_ST                                       */
/*---------------------------------------------------------------------------*/


static int
output_result_header_ST(struct hecmwST_result_data *result, char *header,
		int n_node, int n_elem, FILE *fp)
{
	int rc,len;
	char *p,*q;
	char head[HECMW_HEADER_LEN+1];

	if(header == NULL) {
		head[0] = '\0';
	} else {
		len = 0;
		p = header;
		q = head;
		while(len < sizeof(head)-1 && *p && *p != '\n') {
			*q++ = *p++;
			len++;
		}
		*++q = '\0';
	}

	/* header */
	rc = fprintf(fp, "%s\n", header);
	if(rc < 0) {
		HECMW_set_error(HECMW_UTIL_E0205, "");
		return -1;
	}

	/* n_node, n_elem */
	rc = fprintf(fp, "%d %d\n", n_node, n_elem);
	if(rc < 0) {
		HECMW_set_error(HECMW_UTIL_E0205, "");
		return -1;
	}

	/* nn_component, ne_component */
	rc = fprintf(fp, "%d %d\n", result->nn_component, result->ne_component);
	if(rc < 0) {
		HECMW_set_error(HECMW_UTIL_E0205, "");
		return -1;
	}

	return 0;
}


static int
output_result_node_ST(struct hecmwST_result_data *result, int n_node, FILE *fp)
{
	int i,j,k,n,m,rc;

	/* nn_dof */
	n = 0;
	for(i=0; i < result->nn_component; i++) {
		rc = fprintf(fp, "%d%c", result->nn_dof[i], (n+1)%COL_INT ? ' ' : '\n');
		if(rc < 0) {
			HECMW_set_error(HECMW_UTIL_E0205, "");
			return -1;
		}
		n++;
	}
	if(n % COL_INT) {
		rc = fprintf(fp, "\n");
		if(rc < 0) {
			HECMW_set_error(HECMW_UTIL_E0205, "");
			return -1;
		}
	}

	/* node_label */
	for(i=0; i < result->nn_component; i++) {
		rc = fprintf(fp, "%s\n", result->node_label[i]);
		if(rc < 0) {
			HECMW_set_error(HECMW_UTIL_E0205, "");
			return -1;
		}
	}

	/* node_val_item */
	if(result->nn_component == 0) return 0;
	m = 0;
	for(i=0; i < n_node; i++) {
        rc = fprintf(fp, "%d \n", node_global_ID[i]);
		if(rc < 0) {
			HECMW_set_error(HECMW_UTIL_E0205, "");
			return -1;
		}
		n = 0;
		for(j=0; j < result->nn_component; j++) {
			for(k=0; k < result->nn_dof[j]; k++) {
				rc = fprintf(fp, "%.16E%c", result->node_val_item[m], (n+1)%COL_DOUBLE ? ' ' : '\n');
				n++;
				m++;
			}
		}
		if(n % COL_DOUBLE) {
			rc = fprintf(fp, "\n");
			if(rc < 0) {
				HECMW_set_error(HECMW_UTIL_E0205, "");
				return -1;
			}
		}
	}

	return 0;
}


static int
output_result_elem_ST(struct hecmwST_result_data *result, int n_elem, FILE *fp)
{
	int i,j,k,n,m,rc;

	/* ne_dof */
	n = 0;
	for(i=0; i < result->ne_component; i++) {
		rc = fprintf(fp, "%d%c", result->ne_dof[i], (n+1)%COL_INT ? ' ' : '\n');
		if(rc < 0) {
			HECMW_set_error(HECMW_UTIL_E0205, "");
			return -1;
		}
		n++;
	}
	if(n % COL_INT) {
		rc = fprintf(fp, "\n");
		if(rc < 0) {
			HECMW_set_error(HECMW_UTIL_E0205, "");
			return -1;
		}
	}

	/* elem_label */
	for(i=0; i < result->ne_component; i++) {
		rc = fprintf(fp, "%s\n", result->elem_label[i]);
		if(rc < 0) {
			HECMW_set_error(HECMW_UTIL_E0205, "");
			return -1;
		}
	}

	/* elem_val_item */
	if(result->ne_component == 0) return 0;
	m = 0;
	for(i=0; i < n_elem; i++) {
        rc = fprintf(fp, "%d\n", elem_global_ID[i]);
		if(rc < 0) {
			HECMW_set_error(HECMW_UTIL_E0205, "");
			return -1;
		}
		n = 0;
		for(j=0; j < result->ne_component; j++) {
			for(k=0; k < result->ne_dof[j]; k++) {
				rc = fprintf(fp, "%.16E%c", result->elem_val_item[m], (n+1)%COL_DOUBLE ? ' ' : '\n');
				n++;
				m++;
			}
		}
		if(n % COL_DOUBLE) {
			rc = fprintf(fp, "\n");
			if(rc < 0) {
				HECMW_set_error(HECMW_UTIL_E0205, "");
				return -1;
			}
		}
	}

	return 0;
}


static int
output_result_data_ST(struct hecmwST_result_data *result, int n_node, int n_elem, char *header, FILE *fp)
{
	HECMW_assert(fp);

	if(output_result_header_ST(result, header, n_node, n_elem, fp)) {
		return -1;
	}
	if(output_result_node_ST(result, n_node, fp)) {
		return -1;
	}
	if(output_result_elem_ST(result, n_elem, fp)) {
		return -1;
	}

	return 0;
}


/*---------------------------------------------------------------------------*/
/* TEXT MODE I/O ---- input_result                                           */
/*---------------------------------------------------------------------------*/


static int
get_line(char *buf, int bufsize, FILE *fp)
{
	if(fgets(buf, bufsize, fp) == NULL) {
		HECMW_set_error(HECMW_UTIL_E0205, "");
		return -1;
	}
	return strlen(line_buf);
}


static int
input_result_header(struct hecmwST_result_data *result,
		int *n_node, int *n_elem, FILE *fp)
{
	/* header */
	if(get_line(line_buf, sizeof(line_buf), fp) < 0) {
		return -1;
	}
	line_buf[ strlen(line_buf)-1 ] = 0;/* remove CR/LF*/
	strcpy( head, line_buf );

	/* n_node, n_elem */
	if(get_line(line_buf, sizeof(line_buf), fp) < 0) {
		return -1;
	}
	if(sscanf(line_buf, "%d%d", n_node, n_elem) != 2) {
		HECMW_set_error(HECMW_UTIL_E0205, "");
		return -1;
	}

	/* nn_component, ne_component */
	if(get_line(line_buf, sizeof(line_buf), fp) < 0) {
		return -1;
	}
	if(sscanf(line_buf, "%d%d", &result->nn_component, &result->ne_component) != 2) {
		HECMW_set_error(HECMW_UTIL_E0205, "");
		return -1;
	}

	return 0;
}


static int
input_result_node(struct hecmwST_result_data *result, int n_node, FILE *fp)
{
#define DELIM " \n"
	int i,rc,n;
	int label_counter;
	char *p,*buf;

	if(result->nn_component <= 0) {
		return 0;
	}

	/* nn_dof */
	result->nn_dof = HECMW_malloc(sizeof(*result->nn_dof)*result->nn_component);
	if(result->nn_dof == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}

	if(get_line(line_buf, sizeof(line_buf), fp) < 0) {
		return -1;
	}
	i = n = 0;
	buf = line_buf;
	while(i < result->nn_component) {
		p = strtok(buf, DELIM);
		if(p == NULL) {
			if(get_line(line_buf, sizeof(line_buf), fp) < 0) {
				return -1;
			}
			buf = line_buf;
			continue;
		}
		buf = NULL;
		rc = sscanf(p, "%d", &result->nn_dof[i]);
		if(rc == EOF) {
			HECMW_set_error(HECMW_UTIL_E0204, "");
			return -1;
		}
		if(rc != 1) {
			HECMW_set_error(HECMW_UTIL_E0205, "");
			return -1;
		}
		n += result->nn_dof[i];
		i++;
	}

	/* node_label */
	result->node_label = HECMW_malloc(sizeof(*result->node_label)*result->nn_component);
	if(result->node_label == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}

	for(i=0; i < result->nn_component; i++) {
		char label[HECMW_NAME_LEN+1];
		if(get_line(line_buf, sizeof(line_buf), fp) < 0) {
			return -1;
		}
		rc = sscanf(line_buf, "%s", label);
		if(rc == EOF) {
			HECMW_set_error(HECMW_UTIL_E0204, "");
			return -1;
		}
		if(rc != 1) {
			HECMW_set_error(HECMW_UTIL_E0205, "");
			return -1;
		}
		result->node_label[i] = HECMW_strdup(label);
		if(result->node_label[i] == NULL) {
			HECMW_set_error(errno, "");
			return -1;
		}
	}
	/* node_val_item */
	node_global_ID = HECMW_malloc(sizeof(*node_global_ID)*n_node);
	if(node_global_ID == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	result->node_val_item = HECMW_malloc(sizeof(*result->node_val_item)*n*n_node);
	if(result->node_val_item == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}

	if(get_line(line_buf, sizeof(line_buf), fp) < 0) {
		return -1;
	}
	i = 0;
	label_counter = 0;
	n++;         /**** For global node ID ****/
	buf = line_buf;
	while(i < n*n_node) {
		p = strtok(buf, DELIM);
		if(p == NULL) {
			if(get_line(line_buf, sizeof(line_buf), fp) < 0) {
				return -1;
			}
			buf = line_buf;
			continue;
		}
		buf = NULL;
                if ( (i%n) == 0 ) {
           		rc = sscanf(p, "%d", &node_global_ID[label_counter]);
                label_counter++;
                } else {
        		rc = sscanf(p, "%lf", &result->node_val_item[i-label_counter]);
                }
		if(rc == EOF) {
			HECMW_set_error(HECMW_UTIL_E0204, "");
			return -1;
		}
		if(rc != 1) {
			HECMW_set_error(HECMW_UTIL_E0205, "");
			return -1;
		}
		i++;
	}

	return 0;
}


static int
input_result_elem(struct hecmwST_result_data *result, int n_elem, FILE *fp)
{
#define DELIM " \n"
	int i,rc,n;
	int label_counter;
	char *p,*buf;

	if(result->ne_component <= 0) {
		return 0;
	}

	/* ne_dof */
	result->ne_dof = HECMW_malloc(sizeof(*result->ne_dof)*result->ne_component);
	if(result->ne_dof == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}

	if(get_line(line_buf, sizeof(line_buf), fp) < 0) {
		return -1;
	}
	i = n = 0;
	buf = line_buf;
	while(i < result->ne_component) {

		p = strtok(buf, DELIM);
		if(p == NULL) {
			if(get_line(line_buf, sizeof(line_buf), fp) < 0) {
				return -1;
			}
			buf = line_buf;
			continue;
		}
		buf = NULL;
		rc = sscanf(p, "%d", &result->ne_dof[i]);
		if(rc == EOF) {
			HECMW_set_error(HECMW_UTIL_E0204, "");
			return -1;
		}
		if(rc != 1) {
			HECMW_set_error(HECMW_UTIL_E0205, "");
			return -1;
		}
		n += result->ne_dof[i];
		i++;
	}

	/* elem_label */
	result->elem_label = HECMW_malloc(sizeof(*result->elem_label)*result->ne_component);
	if(result->elem_label == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}

	for(i=0; i < result->ne_component; i++) {
		char label[HECMW_NAME_LEN+1];
		if(get_line(line_buf, sizeof(line_buf), fp) < 0) {
			return -1;
		}
		rc = sscanf(line_buf, "%s", label);
		if(rc == EOF) {
			HECMW_set_error(HECMW_UTIL_E0204, "");
			return -1;
		}
		if(rc != 1) {
			HECMW_set_error(HECMW_UTIL_E0205, "");
			return -1;
		}
		result->elem_label[i] = HECMW_strdup(label);
		if(result->elem_label[i] == NULL) {
			HECMW_set_error(errno, "");
			return -1;
		}
	}

	/* elem_val_item */
	elem_global_ID = HECMW_malloc(sizeof(*elem_global_ID)*n_elem);
	if(elem_global_ID == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	result->elem_val_item = HECMW_malloc(sizeof(*result->elem_val_item)*n*n_elem);
	if(result->elem_val_item == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}

	if(get_line(line_buf, sizeof(line_buf), fp) < 0) {
		return -1;
	}
	i = 0;
    label_counter = 0;
	n++;          /**** For global element ID ****/
	buf = line_buf;
	while(i < n*n_elem) {
		p = strtok(buf, DELIM);
		if(p == NULL) {
			if(get_line(line_buf, sizeof(line_buf), fp) < 0) {
				return -1;
			}
			buf = line_buf;
			continue;
		}
		buf = NULL;
                if ( (i%n) == 0 ) {
           		rc = sscanf(p, "%d", &elem_global_ID[label_counter]);
                label_counter++;
                } else {
         		rc = sscanf(p, "%lf", &result->elem_val_item[i-label_counter]);
                }
		if(rc == EOF) {
			HECMW_set_error(HECMW_UTIL_E0204, "");
			return -1;
		}
		if(rc != 1) {
			HECMW_set_error(HECMW_UTIL_E0205, "");
			return -1;
		}
		i++;
	}

	return 0;
}


static struct hecmwST_result_data *
input_result_data(FILE *fp)
{
	int n_node, n_elem;
	struct hecmwST_result_data *result;

	HECMW_assert(fp);

	result = HECMW_calloc(1, sizeof(*result));
	if(result == NULL) {
		HECMW_set_error(errno, "");
		return NULL;
	}
	if(input_result_header(result, &n_node, &n_elem, fp)) {
		return NULL;
	}
	nnode = n_node;
	nelem = n_elem;
	if(input_result_node(result, n_node, fp)) {
		return NULL;
	}
	if(input_result_elem(result, n_elem, fp)) {
		return NULL;
	}

	return result;
}

/*---------------------------------------------------------------------------*/

struct hecmwST_result_data *
result_read_txt_by_fname(char *filename)
{
	FILE *fp;
	struct hecmwST_result_data *result;

	if((fp = fopen(filename, "r")) == NULL) {
		HECMW_set_error(HECMW_UTIL_E0201, "File: %s, %s", filename, HECMW_strmsg(errno));
		return NULL;
	}

	result = input_result_data(fp);
	if(result == NULL) {
		return NULL;
	}

	if(fclose(fp)) {
		HECMW_set_error(HECMW_UTIL_E0202, "");
		return NULL;
	}

	return result;
}

