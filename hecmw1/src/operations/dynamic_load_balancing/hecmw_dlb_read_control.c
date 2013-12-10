/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 1.00                                              *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : Dynamic Load Balancing                            *
 *                                                                     *
 *            Written by Li Chen (Univ. of Tokyo)                      *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using Hight End Computing Middleware (HEC-MW)"     *
 *                                                                     *
 *=====================================================================*/

#include "hecmw_repart.h"

int HECMW_dlb_is_blank_line(char *buf){
   int i,j;
   int flag;

   i=0;
   flag=1;
   while(buf[i]!='\n') {
	   if(buf[i]!=' ') {
		   flag=0;
		   break;
	   }
	   i++;
   }
   return(flag);
}

int HECMW_dlb_is_comment_line(char *buf) {
	int  i,j;
	int  flag;
    flag=0;
	if(buf[0]=='#')
		flag=1;
	else if((buf[0]=='!') && (buf[1]=='!'))
		flag=1;
	return(flag);
}


void HECMW_dlb_get_string_item(char *para, char *buf, int *start_location, char para2[128])
{
   int value;
   int  i,j;

   i=*start_location;
   while((buf[i]==',') || (buf[i]==' ') || (buf[i]=='='))
	   i++;
   if(buf[i]=='\n') {
	   fprintf(stderr, "No string value for %s\n", para);
	   HECMW_dlb_print_exit("ERROR: HEC-MW-VIS-E0002: The control parameter format error: should start from !");
   }
   j=0;
    	while((buf[i]!=' ') && (buf[i]!=',') && (buf[i]!='\n')) {
	    	para2[j]=buf[i];
    		i++;
    		j++;
		}
		para2[j]='\0';
		*start_location=i;
		return;
}

int HECMW_dlb_get_int_item(char *para, char *buf, int *start_location)
{
   int value;
   int  i,j;
   char para2[128];

   i=*start_location;
   while((buf[i]==',') || (buf[i]==' ') || (buf[i]=='='))
	   i++;
   if(buf[i]=='\n') {
	   HECMW_dlb_print_exit("ERROR: HEC-MW-VIS-E0003:The control parameter format error:No integer value for %s");
   }
   j=0;
    	while((buf[i]!=' ') && (buf[i]!=',') && (buf[i]!='\n')) {
	    	para2[j]=buf[i];
    		i++;
    		j++;
		}
		para2[j]='\0';
		  if((isdigit(para2[0])==0) && (para2[0]!='+') && (para2[0]!='-')) {
			  fprintf(stderr, "ERROR: HEC-MW-VIS-E0004: %s should be integer \n", para);
			  HECMW_dlb_print_exit("Please re-input and run again!");
		  }
		value=atoi(para2);
		*start_location=i;
		return(value);
}
       

double HECMW_dlb_get_double_item(char *para, char *buf, int *start_location)
{
   double value;
   int  i,j;
   char para2[128];

   i=*start_location;
   while((buf[i]==',') || (buf[i]==' ') || (buf[i]=='='))
	   i++;
   if(buf[i]=='\n') {
	   fprintf(stderr, "No integer value for %s\n", para);
	   HECMW_dlb_print_exit("The control parameter format error:!");
   }
   j=0;
    	while((buf[i]!=' ') && (buf[i]!=',') && (buf[i]!='\n')) {
	    	para2[j]=buf[i];
    		i++;
    		j++;
		}
		para2[j]='\0';
		  if((isdigit(para2[0])==0) && (para2[0]!='+') && (para2[0]!='-')) {
			  fprintf(stderr, "ERROR: HEC-MW-VIS-E0005:%s should be a real \n", para);
			  HECMW_dlb_print_exit("Please re-input and run again!");
		  }
		value=atof(para2);
		*start_location=i;
		return(value);
}
       
int HECMW_dlb_get_keyword_item(char *buf, char *para) {
	int i,j;
	int flag;
    i=0;
	while(buf[i]==' ')
		i++;
	if(buf[i]!='!') {
		fprintf(stderr, "Please check the line %s\n", buf);
		HECMW_dlb_print_exit("The control parameter format error:!");
	}
	i++;
	j=0;
	while((buf[i]!=' ') && (buf[i]!='=') && (buf[i]!=',') && (buf[i]!='\n')) {
		para[j]=buf[i];
		i++;
		j++;
	}
	para[j]='\0';
	return (i);
}

int HECMW_dlb_get_keyword_repart(char *buf, Control_para *ctl_para) {
	int i,j;
	int flag;
	char para[128], para2[128];

    i=0;
	while(buf[i]==' ')
		i++;
	if(buf[i]!='!') 
		HECMW_dlb_print_exit("ERROR: HEC-MW-DLB-E0002: The control parameter format error: should start from !");
	i=i+1;  j=0;
	while((buf[i]!=' ') && (buf[i]!=',') && (buf[i]!='\n')) {
		para[j]=buf[i];
		i++;
		j++;
	}
	flag=1;
	if((strncmp(para, "DLB_CTRL", 8)!=0) && (strncmp(para, "dlb_ctrl", 8)!=0)) {
		flag=0;
		return (flag);
	}
/*	while((buf[i]!='=') && (buf[i]!='\n'))
		i++;
	if(buf[i]=='\n')
		HECMW_dlb_print_exit("ERROR: HEC-MW-DLB-E0006:The control parameter format error: DLB_METHOD");
	i++;
	while((buf[i]==' ') && (buf[i]!='\n'))
		i++;
	if(buf[i]=='\n')
		HECMW_dlb_print_exit("ERROR: HEC-MW-VIS-E0006:The control parameter format error: visual ID");
    j=0;
	while((buf[i]!=' ') && (buf[i]!=',') && (buf[i]!='=') && (buf[i]!='\n')) {
		para[j]=buf[i];
		i++;
		j++;
	}
	if((strncmp(para, "PSR", 3)==0) || (strncmp(para, "psr", 3)==0)) {
		flag=1;
		return (flag);
	}

	if((strncmp(para, "PVR", 3)==0) || (strncmp(para, "pvr", 3)==0)) {
		flag=2;
		return (flag);
	}

	HECMW_dlb_print_exit("ERROR: HEC-MW-VIS-E0007:The control parameter format error: method only can be PSR or PVR");
	*/
	strncpy(ctl_para->adaptive_repartition,"off", 3);

	while((buf[i]==',') || (buf[i]==' '))
		  i++;
	while(buf[i]!='\n') {
         j=0; 
    	while((buf[i]!=' ') && (buf[i]!=',') && (buf[i]!='=') && (buf[i]!='\n')) {
	    	para[j]=buf[i];
    		i++;
    		j++;
		}
		para[j]='\0';
		if((strncmp(para, "method", 6)==0) || (strncmp(para, "METHOD", 6)==0)) {
			HECMW_dlb_get_string_item(para, buf, &i, para2);
			if((strncmp(para, "ADAPT", 5)==0) || (strncmp(para, "adapt", 5)==0))
	           strncpy(ctl_para->adaptive_repartition,"on", 2);
		}
				

	  	  while((buf[i]==',') || (buf[i]==' '))
		  i++;
	  }
	return flag;

}



void hecmw_dlb_read_control(char *contfile, Control_para *ctl_para, int stat_para[NUM_CONTROL_PARAS], int pesize)
{
  int i, j, k, offset;
  char		buf[MAX_LINE_LEN];
  char*     parameters[NUM_CONTROL_PARAS];
   int      len_para[NUM_CONTROL_PARAS];

   char     para[128], para1[128];
   int      hit;
   int      surface_num;
   int      len_str; 
   int     location, visual_method;
   int flag, flag_surface;
   int    cont_flag;

   FILE     *fp;
   
   fp=fopen(contfile, "r");
   if(fp==NULL) 
	   HECMW_dlb_print_exit("Cannot find the control input file");
   for(i=0;i<NUM_CONTROL_PARAS;i++) {
	   parameters[i]=(char *)calloc(128, sizeof(char));
	   if(parameters[i]==NULL) {
		   fprintf(stderr, "There is no enough memory for parameters\n");
		   exit(0);
	   }
   }

  parameters[0]="adaptive_repartition";
  len_para[0]=11;
  parameters[1]="num_of_criteria";
  len_para[1]=12;
  parameters[2]="balance_rate";
  len_para[2]=10;
  parameters[3]="num_of_repartition";
  len_para[3]=15;
  parameters[4]="itr_rate";
  len_para[4]=8;
  parameters[5]="wgtflag";
  len_para[5]=7;
  parameters[6]="vwgt_filename";
  len_para[6]=9;
  parameters[7]="adjwgt_filename";
  len_para[7]=12;
  parameters[8]="machine_wgt";
  len_para[8]=11;



  for(i=0;i<NUM_CONTROL_PARAS;i++) {
	  stat_para[i]=0;
  }
  ctl_para->num_criteria=1;

  offset=0;
  cont_flag=1;
  while (cont_flag) {
	  if(fgets(buf, MAX_LINE_LEN, fp) != NULL) {
        if ((HECMW_dlb_is_blank_line(buf)==0) && (HECMW_dlb_is_comment_line(buf)==0)) break;
	  }
	  else 
		  cont_flag=0;
  }
  
  hit=0;
  if(cont_flag==1)
  cont_flag=HECMW_dlb_get_keyword_repart(buf,ctl_para);
  while (cont_flag==1) {
	  if(fgets(buf, MAX_LINE_LEN, fp) != NULL) {
        if ((HECMW_dlb_is_blank_line(buf)==0) && (HECMW_dlb_is_comment_line(buf)==0)) break;
	  }
	  else
		  cont_flag=0;
  }

  while(cont_flag==1) {
	  hit=-1;
	  location=HECMW_dlb_get_keyword_item(buf, para);
	  for(i=1;i<NUM_CONTROL_PARAS;i++) {
		  if((strncmp(para, parameters[i], len_para[i]))==0) {
			  hit=i;
			  stat_para[i]=1;
			  break;
		  }
	  }
/*	 fprintf(stderr, "para=%s hit=%d\n", para, hit); 
*/
         if((hit>=0) && (hit<NUM_CONTROL_PARAS)) {
         switch (hit) {
/*	  case 0:
		  sf[k].surface_style=get_int_item(para, buf, &location);
		  break;
		  */
	  case 1:
		  ctl_para->num_criteria=HECMW_dlb_get_int_item(para, buf, &location);
		  break;
	  case 2:
		  ctl_para->balance_rate=(float *)calloc(ctl_para->num_criteria, sizeof(float));
		  for(i=0;i<ctl_para->num_criteria;i++) {
		  
		  ctl_para->balance_rate[i]=(float)HECMW_dlb_get_double_item(para, buf, &location);
		  }
		  break;
	  case 3:
          ctl_para->num_repartition=HECMW_dlb_get_int_item(para, buf, &location);
		  break;
	  case 4:
		  ctl_para->itr_rate=(float)HECMW_dlb_get_double_item(para, buf, &location);
              break;
	  case 5:
		  ctl_para->wgtflag=HECMW_dlb_get_int_item(para, buf, &location);
		  break;
	  case 6:
		  HECMW_dlb_get_string_item(para, buf, &location, ctl_para->vwgt_filename);
		  break;
	  case 7:
		  HECMW_dlb_get_string_item(para, buf, &location, ctl_para->adjwgt_filename);
		  break;
	  case 8:
		  ctl_para->machine_wgt=(float *)calloc(pesize, sizeof(float));
		  for(i=0;i<pesize;i++) {
		  ctl_para->machine_wgt[i]=(float)HECMW_dlb_get_double_item(para, buf, &location);
		  }
		  break;
/*	  case 9:
		  fscanf(contfp, "%s", ctl_para->output_filename);
		  */

		  break;


		  }
              }
         while (cont_flag) {
	        if(fgets(buf, MAX_LINE_LEN, fp) != NULL) {
                if ((HECMW_dlb_is_blank_line(buf)==0) && (HECMW_dlb_is_comment_line(buf)==0)) break;
			}
			else 
		       cont_flag=0;
		 }
		if(cont_flag==0)
		 break;
        }
		  /* check the parameters */
/*		  if(stat_para[0]==0) {
			  strncpy(ctl_para->adaptive_repartition, "off", 3);
		  }
		  
		  if((strncmp(ctl_para->adaptive_repartition, "on", 2)!=0) && (strncmp(ctl_para->adaptive_repartition, "off", 3)!=0)){
			  fprintf(stderr, "adaptive_repartition should be on or off\n");
			  fprintf(stderr, "Please re-input again\n");
			  exit(0);
		  }
*/
          if(stat_para[1]==0) {
			  ctl_para->num_criteria=1;
		  }
		  
		  if(ctl_para->num_criteria<=0) {
			  fprintf(stderr, "#### HEC-MW-DLB-E1001: num_of_criteria should be greater than 0\n");
			  fprintf(stderr, "Please re-input a correct one\n");
			  exit(0);
		  }
		  if(stat_para[2]==0) {
			  ctl_para->balance_rate=(float *)calloc(ctl_para->num_criteria, sizeof(float));
			  for(i=0;i<ctl_para->num_criteria;i++)
				  ctl_para->balance_rate[i]=1.05;
		  }
		  for(i=0;i<ctl_para->num_criteria;i++) {
			  if(ctl_para->balance_rate[i]<1.0) {
				  fprintf(stderr, "#### HEC-MW-DLB-E1002: The balance rate should be >=1.0\n");
				  fprintf(stderr, "Please input again\n");
				  exit(0);
			  }
		  }

		  if(stat_para[3]==0) {
			  ctl_para->num_repartition=pesize;
		  }
		  if(ctl_para->num_repartition<1) {
			  fprintf(stderr, "#### HEC-MW-DLB-E1003: The num_of_repartition cannot less than 1\n");
			  fprintf(stderr, "Please input again\n");
			  exit(0);
		  }
		  if(stat_para[4]==0) {
		     ctl_para->itr_rate=10000.0;
		  }
		  if(ctl_para->itr_rate<0.0) {
			  fprintf(stderr, "#### HEC-MW-DLB-E1004:itr_rate cannot be less than 0.0\n");
			  exit(0);
		  }

		  if(stat_para[5]==0) {
			  ctl_para->wgtflag=0;
		  }
		  if((ctl_para->wgtflag<0) || (ctl_para->wgtflag>3)) {
			  fprintf(stderr, "#### HEC-MW-DLB-E1005:wgtflag only can be in 0--3\n");
			  exit(0);
		  }
		  if(stat_para[8]==0) {
			  ctl_para->machine_wgt=(float *)calloc(pesize, sizeof(float));
			  for(i=0;i<pesize;i++)
				  ctl_para->machine_wgt[i]=1.0/(float)pesize;
		  }

		  return;

}

void hecmw_dlb_set_default_control(Control_para *ctl_para, int stat_para[NUM_CONTROL_PARAS], int pesize)
{
	int i;
			  
	strncpy(ctl_para->adaptive_repartition,"on", 3);
	ctl_para->num_criteria=1;
		  
	ctl_para->balance_rate=(float *)calloc(1, sizeof(float));
	ctl_para->balance_rate[0]=1.02;		
	ctl_para->num_repartition=pesize;
	ctl_para->itr_rate=10000.0;
    ctl_para->wgtflag=0;
    for(i=6;i<NUM_CONTROL_PARAS;i++)
		stat_para[i]=0;
	stat_para[8]=1;
	ctl_para->machine_wgt=(float *)calloc(pesize, sizeof(float));
	for(i=0;i<pesize;i++)
		ctl_para->machine_wgt[i]=1.0/(float)pesize;

	return;

}






