/* ----------------------------------------------------------------  */
/* readFoam                                        7.07.2005 Wittig  */
/* ----------------------------------------------------------------  */


#include <cgx.h>
#include <sys/dir.h>

#define kappa    1.4

#define MAXINTEGER 2147483647
#define MAXFLOAT   1.e32

#define TEST     0

#define MAX_ENTITIES      100
#define CHECK_INPUT       0

extern Alias     *alias;
extern Sets      *set;
extern SumGeo    anzGeo[1];
extern SumAsci   sumAsci[1];


typedef struct {
  int anz_f;
  int *face;
  int *ori;
}Elemf;

typedef struct {
  int p;
  int u;
  int t;
  int rho;
}Ds;


int get_values_in_string( char *handle, double *fval)
{
  int i=0, pos=0;
  static int parenthesis=0;
  char byte;
  char string[MAX_LINE_LENGTH];


  /* count open and closed parenthesis */
  do
  {
    byte = handle[pos++];
    if(byte=='(') parenthesis++;
    else if(byte==')')
    {
      parenthesis--;

      /* return if the first open '(' is closed */
      if(parenthesis==0) break;

      /* if one "(" is open get the values and return */
      if(parenthesis==1)
      {
        string[i++]=0;
        sscanf(string, "%lf %lf %lf", &fval[0], &fval[1], &fval[2]);
        return(pos);
      }
    }
    else if (parenthesis==2)
    {
      string[i++]=byte;
    }

  }while(byte!=0);
  fval[0]=fval[1]=fval[2]=0.;
  return(-parenthesis);
}

int defineElemNodes(Summen *anz, Faces *face, Elements *elem, int i, int *eface, int anz_f)
{
  int f, l, n, j, k, ni[8];
  int fbuf[6], slavnod[2];

  if (anz_f==6) /* HEXA8  */
  {
    /* mark the faces as used */
    for(f=0; f<anz_f; f++) face[eface[f]].group++;
  
    elem[i].type = 1;
    anz->etype[elem[i].type]++;
  
    /* bottom */
    /* owner */
    if(face[eface[0]].group==1)
    {
      elem[i].nod[0] = face[eface[0]].nod[0];         
      elem[i].nod[4] = face[eface[0]].nod[1];         
      elem[i].nod[7] = face[eface[0]].nod[2];         
      elem[i].nod[3] = face[eface[0]].nod[3];
    }
    /* neighbour */
    else if(face[eface[0]].group==2)
    {
      elem[i].nod[0] = face[eface[0]].nod[0];         
      elem[i].nod[4] = face[eface[0]].nod[3];         
      elem[i].nod[7] = face[eface[0]].nod[2];         
      elem[i].nod[3] = face[eface[0]].nod[1];
    }
    else { printf("ERROR in group:%d\n",face[eface[0]].group); exit(1); }
  
    /* search face on node 0 and 1 from face0 */
    for(f=2; f<anz_f; f++)
    { 
      for(l=0; l<4; l++)
      {
        if(face[eface[f]].nod[l]==face[eface[0]].nod[0])
        {
          slavnod[0]=l;
          goto found_face;
        }
      }
    }
    printf("ERROR: found no slavface for elem:%d\n",i ); exit(1);
  found_face:;
    for(l=0; l<4; l++)
    {
      for(n=1; n<4; n++)
      {
        if(face[eface[f]].nod[l]==face[eface[0]].nod[n])
        {
          if((l==0)&&(slavnod[0]==1)) slavnod[1]=2;
        else if((l==1)&&(slavnod[0]==2)) slavnod[1]=3;
        else if((l==2)&&(slavnod[0]==3)) slavnod[1]=0;
        else if((l==3)&&(slavnod[0]==0)) slavnod[1]=1;
  
        else if((l==2)&&(slavnod[0]==1)) slavnod[1]=0;
        else if((l==3)&&(slavnod[0]==2)) slavnod[1]=1;
        else if((l==0)&&(slavnod[0]==3)) slavnod[1]=2;
        else if((l==1)&&(slavnod[0]==0)) slavnod[1]=3;
        else { printf("ERROR: slavnod nod found l:%d slavnod[0]:%d\n",l,slavnod[0]); exit(1); }
          goto found_nodn;
        }
      }
    }
    printf("ERROR: found no slavnod[1] for elem:%d\n",i ); exit(1);
  found_nodn:;
  
    /* determine the node from face1 at the slavnod[1] */
    for(l=0; l<4; l++)
    {
      if(face[eface[f]].nod[slavnod[1]]==face[eface[1]].nod[l])
      {
        /* top */
        if(face[eface[1]].group==1)
        {
          if(l==0)
          {
            elem[i].nod[1] = face[eface[1]].nod[0];         
            elem[i].nod[2] = face[eface[1]].nod[1];         
            elem[i].nod[6] = face[eface[1]].nod[2];         
            elem[i].nod[5] = face[eface[1]].nod[3];
          }        
          else if(l==1)
          {
            elem[i].nod[1] = face[eface[1]].nod[1];         
            elem[i].nod[2] = face[eface[1]].nod[2];         
            elem[i].nod[6] = face[eface[1]].nod[3];         
            elem[i].nod[5] = face[eface[1]].nod[0];
          }        
          else if(l==2)
          {
            elem[i].nod[1] = face[eface[1]].nod[2];         
            elem[i].nod[2] = face[eface[1]].nod[3];         
            elem[i].nod[6] = face[eface[1]].nod[0];         
            elem[i].nod[5] = face[eface[1]].nod[1];
          }        
          else if(l==3)
          {
            elem[i].nod[1] = face[eface[1]].nod[3];         
            elem[i].nod[2] = face[eface[1]].nod[0];         
            elem[i].nod[6] = face[eface[1]].nod[1];         
            elem[i].nod[5] = face[eface[1]].nod[2];
          }        
          else { printf("ERROR: l nod found\n"); exit(1); }
        }
        else if(face[eface[1]].group==2)
        {
          if(l==0)
          {
            elem[i].nod[1] = face[eface[1]].nod[0];         
            elem[i].nod[2] = face[eface[1]].nod[3];         
            elem[i].nod[6] = face[eface[1]].nod[2];         
            elem[i].nod[5] = face[eface[1]].nod[1];
          }         
          else if(l==1)
          {
            elem[i].nod[1] = face[eface[1]].nod[1];         
            elem[i].nod[2] = face[eface[1]].nod[0];         
            elem[i].nod[6] = face[eface[1]].nod[3];         
            elem[i].nod[5] = face[eface[1]].nod[2];
          }         
          else if(l==2)
          {
            elem[i].nod[1] = face[eface[1]].nod[2];         
            elem[i].nod[2] = face[eface[1]].nod[1];         
            elem[i].nod[6] = face[eface[1]].nod[0];         
            elem[i].nod[5] = face[eface[1]].nod[3];
          }         
          else if(l==3)
          {
            elem[i].nod[1] = face[eface[1]].nod[3];         
            elem[i].nod[2] = face[eface[1]].nod[2];         
            elem[i].nod[6] = face[eface[1]].nod[1];         
            elem[i].nod[5] = face[eface[1]].nod[0];
          }         
          else { printf("ERROR: l nod found\n"); exit(1); }
        }
        else { printf("ERROR in group:%d\n",face[eface[1]].group); exit(1); }
        break;
      }
    }
  
    for(f=0; f<anz_f; f++)
    {
      face[eface[f]].elem_nr=elem[i].nr;
  
      n=0;
      for(j=0; j<8; j++) 
        for(k=0; k<4; k++)
          if(elem[i].nod[j]==face[eface[f]].nod[k])
            ni[n++]=j+1;
      if (n==4)
      {
        /* check which sides are involved */
        for (j=0; j<6; j++) fbuf[j]=0;
        for (j=0; j<4; j++)
        {
          if ((ni[j]==1)||(ni[j]==2)||(ni[j]==3)||(ni[j]==4)) fbuf[0]++;
          if ((ni[j]==5)||(ni[j]==8)||(ni[j]==7)||(ni[j]==6)) fbuf[1]++;
            if ((ni[j]==1)||(ni[j]==5)||(ni[j]==6)||(ni[j]==2)) fbuf[2]++;
          if ((ni[j]==2)||(ni[j]==6)||(ni[j]==7)||(ni[j]==3)) fbuf[3]++;
          if ((ni[j]==3)||(ni[j]==7)||(ni[j]==8)||(ni[j]==4)) fbuf[4]++;
          if ((ni[j]==4)||(ni[j]==8)||(ni[j]==5)||(ni[j]==1)) fbuf[5]++;
        }
      }
      for (j=0; j<6; j++)
      {
        /* printf(" j:%d fbuf[j]:%d\n", j, fbuf[j]); */
        if(fbuf[j]==4)
        {
          face[eface[f]].nr=j;
          if( anz->emax>anz->f) 
          {
            if ( (face = (Faces *)realloc((Faces *)face, (anz->emax+2) * sizeof(Faces))) == NULL )
            {  printf("\n\n ERROR: malloc failed, face\n\n") ; exit(1); }
	  }
          face[elem[i].nr].indx[j]=eface[f];
	}
      }
    }
  }
  else if (anz_f==5) /* PE6    */
  {
    elem[i].type = 2; printf("elem-type:pe6 not installed\n"); return(0);
  }
  else if (anz_f==4) /* TET4   */
  {
    elem[i].type = 3; printf("elem-type:tet4 not installed\n"); return(0); 
  }
  else
  {
    printf("\n\n ERROR: number of faces:%d of element:%d do not match available types of elements\n", anz_f, i);
    return(0);
  }
  
#if TEST
          printf (" e:%d t:%d f:",elem[i].nr, elem[i].type);
          for(j=0; j<anz_f; j++)  printf (" %d",eface[j]);
          printf ("\nnode:");
          for(j=0; j<8; j++)  printf (" %d",elem[i].nod[j]);
          printf ("\n");
#endif
   return(1);
}





int readFoam(char *datin, Summen *anz, Nodes **nptr, Elements **eptr, Datasets **lptr  )
{
  int i,j,k=0,l,n,f,f1=0,f2, anz_f, length, pos;
  int eface[8], sum_dir, ncomps;
  FILE *handle;
  DIR *dirp;
  struct direct *dp;
  Rsort *rsort=NULL;
  Ds ds;

  double cp, Rg;
  double ps, us, ts, pt, tt, rho, ma;


  char file[MAX_LINE_LENGTH], rec_str[MAX_LINE_LENGTH];
  char dat[2][MAX_LINE_LENGTH];
  char resultname[MAX_LINE_LENGTH];
  char buffer[MAX_LINE_LENGTH];
  char **resultdir=NULL;
  char **resultdir2=NULL;

  Nodes     *node = NULL;
  Faces     *face = NULL;
  Elements  *elem = NULL;
  Datasets  *lcase= NULL;
  Elemf     *elemf= NULL;

  int anz_dat=0, ipuf=0, *vpn=NULL;
  double  fval[10];

  int startFace, nFaces, setNr=-1, formFlag=0, foundFlag;

  anz->n=-1;
  anz->f=-1;
  anz->e=-1;
  anz->l=-1;

  /* Open the files and check to see that it was opened correctly */

  /* mesh */
  strcpy(anz->model, "foam");
  printf (" MODEL NAME:  %s\n", anz->model);

  anz->n=anz->f=anz->e=0;
  anz->nmax=-MAXINTEGER;  anz->nmin=MAXINTEGER;
  anz->emax=-MAXINTEGER;  anz->emin=MAXINTEGER;

  printf ("\n Read the Mesh  \n");

  /* nodes */
  /* located in "datin"/constant/polymesh/points */
  sprintf( file, "%s/constant/polyMesh/points", datin);
  handle = fopen (file, "r");
  if ( handle== NULL )
  {
    printf ("\nThe input file \"%s\" could not be opened.\n", file);
  }
  else printf ("\n%s opened\n",file);
  while(1)
  {
    length = frecord( handle, rec_str);
    if (rec_str[length] == (char)EOF) break;

    /* read data if a block starts */
    if(rec_str[0]=='(')
    {
      if ( (node = (Nodes *)malloc( (anz->n+1) * sizeof(Nodes))) == NULL )
      {  printf("\n\n ERROR: malloc failed, node\n\n") ; exit(1); }
      for(i=0; i<anz->n; i++)
      {
        node[i].nr=i+1;
        node[node[i].nr].indx=i;
        length = frecord( handle, rec_str);
        sscanf(&rec_str[1], "%lf %lf %lf", &node[node[i].nr].nx, &node[node[i].nr].ny, &node[node[i].nr].nz);
#if TEST
        printf (" n=%d x=%lf y=%lf z=%lf \n", node[i].nr, node[node[i].nr].nx, node[node[i].nr].ny, node[node[i].nr].nz); 
#endif
      }
      break;
    }
    else if(length>0) anz->n=atoi(rec_str);
  }
  fclose(handle);
  anz->nmax=anz->n;
  anz->nmin=1;

  /* faces */
  /* located in "datin"/constant/polymesh/faces */
  sprintf( file, "%s/constant/polyMesh/faces", datin);
  handle = fopen (file, "r");
  if ( handle== NULL )
  {
    printf ("\nThe input file \"%s\" could not be opened.\n", file);
  }
  else printf ("\n%s opened\n",file);
  while(1)
  {
    length = frecord( handle, rec_str);
    if (rec_str[length] == (char)EOF) break;

    /* read data if a block starts */
    if(rec_str[0]=='(')
    {
      if ( (face = (Faces *)malloc( (anz->f+2) * sizeof(Faces))) == NULL )
      {  printf("\n\n ERROR: malloc failed, face\n\n") ; exit(1); }
      for(i=0; i<anz->f; i++)
      {
        face[i].group=0;
        face[i].mat=-1;
        face[i].side=NULL;
        length = frecord( handle, rec_str);
        sscanf(rec_str, "%d",&n);
        sscanf(&rec_str[2], "%d %d %d %d %d %d %d %d",&face[i].nod[0],&face[i].nod[1],&face[i].nod[2],&face[i].nod[3],&face[i].nod[4],&face[i].nod[5],&face[i].nod[6],&face[i].nod[7]);
        for(j=0; j<n; j++)  face[i].nod[j]++; /* nodenr is one higher! */

        if(n==3)      face[i].type = 7;  /* TRI3  */ 
        else if(n==6) face[i].type = 8;  /* TRI6  */  
        else if(n==4) face[i].type = 9;  /* QUAD4 */  
        else if(n==8) face[i].type = 10; /* QUAD8 */
        else
        {
          printf("\n\n ERROR: number of nodes:%d of face:%d do not match available types of elements\n",n,i);
          exit(1);
        }

#if TEST
        printf (" f:%d t:%d",i, face[i].type);
        for(j=0; j<n; j++)  printf (" %d",face[i].nod[j]);
        printf ("\n");
#endif
      }
      break;
    }
    else if(length>0) anz->f=atoi(rec_str);
  }
  fclose(handle);

  /* elements */
  /* located in "datin"/constant/polymesh/cells */
  sprintf( file, "%s/constant/polyMesh/cells", datin);
  handle = fopen (file, "r");
  if ( handle== NULL )
  {
    //printf ("\nThe input file \"%s\" could not be opened.\n", file);

    /* no cells found, try to read owner and neighbour */
    sprintf( file, "%s/constant/polyMesh/owner", datin);
    handle = fopen (file, "r");
    if ( handle== NULL )
    {
      printf("ERROR: found no owner file\n"); exit(1);
    }

    /* create a list of faces per element */
    printf ("\n%s opened\n",file);
    while(1)
    {
      length = frecord( handle, rec_str);
      if (rec_str[length] == (char)EOF) break;
  
      if(rec_str[0]=='(') break;
    }
    /* read data if a block starts */
    i=0;
    while(1)
    {
      length = frecord( handle, rec_str);
      if (rec_str[length] == (char)EOF) break;
      if(rec_str[0]==')') break;
      sscanf(rec_str, "%d",&face[i].elem_nr);
      if(face[i].elem_nr>anz->e) anz->e=face[i].elem_nr;
      i++;
    }
    fclose(handle);
    anz->e++;

    sprintf( file, "%s/constant/polyMesh/neighbour", datin);
    handle = fopen (file, "r");
    if ( handle== NULL )
    {
      printf("ERROR: found no owner file\n"); exit(1);
    }

    /* create a list of faces per element */
    printf ("\n%s opened\n",file);
    while(1)
    {
      length = frecord( handle, rec_str);
      if (rec_str[length] == (char)EOF) break;
  
      if(rec_str[0]=='(') break;
    }
    /* read data if a block starts */
    i=0;
    while(1)
    {
      length = frecord( handle, rec_str);
      if (rec_str[length] == (char)EOF) break;
      if(rec_str[0]==')') break;
      sscanf(rec_str, "%d",&face[i++].mat);  
    }
    fclose(handle);

    /* define the cells */
    if ( (elemf = (Elemf *)malloc( (anz->e+1) * sizeof(Elemf))) == NULL )
    {  printf("\n\n ERROR: malloc failed, elemf\n\n") ; exit(1); }
    for(i=0; i<anz->e; i++) elemf[i].face=NULL;
    for(i=0; i<anz->e; i++) elemf[i].anz_f=0;
    for(i=0; i<anz->f; i++)
    {
      if ( (elemf[face[i].elem_nr].face = (int *)realloc( (int *)elemf[face[i].elem_nr].face, (elemf[face[i].elem_nr].anz_f+1) * sizeof(int))) == NULL )
      {  printf("\n\n ERROR: malloc failed, elemf\n\n") ; exit(1); }

      elemf[face[i].elem_nr].face[elemf[face[i].elem_nr].anz_f]=i;
      elemf[face[i].elem_nr].anz_f++;
      if(face[i].mat>-1)
      {
        if ( (elemf[face[i].mat].face = (int *)realloc( (int *)elemf[face[i].mat].face, (elemf[face[i].mat].anz_f+1) * sizeof(int))) == NULL )
        {  printf("\n\n ERROR: malloc failed, elemf\n\n") ; exit(1); }
        elemf[face[i].mat].face[elemf[face[i].mat].anz_f]=i;
        elemf[face[i].mat].anz_f++;
      }
    }

    /* define the elements */
    if ( (elem = (Elements *)malloc( (anz->e+1) * sizeof(Elements))) == NULL )
    {  printf("\n\n ERROR: malloc failed, elem\n\n") ; exit(1); }
    for(i=0; i<anz->e; i++)
    {
      /* search two disjunct faces */
      eface[0]=elemf[i].face[0];
      f2=2;
      if     (face[eface[0]].type == 7) f1=3;  /* TRI3  */ 
      else if(face[eface[0]].type == 8) f1=6;  /* TRI6  */  
      else if(face[eface[0]].type == 9) f1=4;  /* QUAD4 */  
      else if(face[eface[0]].type == 10) f1=8; /* QUAD8 */
      for(f=1; f<elemf[i].anz_f; f++)
      {
	if     (face[elemf[i].face[f]].type == 7) k=3;  /* TRI3  */ 
        else if(face[elemf[i].face[f]].type == 8) k=6;  /* TRI6  */  
        else if(face[elemf[i].face[f]].type == 9) k=4;  /* QUAD4 */  
        else if(face[elemf[i].face[f]].type == 10) k=8; /* QUAD8 */
        for(n=0; n<k; n++)
	{
	  for(j=0; j<f1; j++)
          {
	    if(face[elemf[i].face[f]].nod[n]==face[eface[0]].nod[j])
            {
              eface[f2++]=elemf[i].face[f];
              goto isconnected;
            }
          }
        }
        eface[1]=elemf[i].face[f];
      isconnected:;
      }
#if TEST
      printf("face order (%d): ", elemf[i].anz_f);
      for(f=0; f<elemf[i].anz_f; f++) printf(" %d", eface[f]);
      printf("\n");
#endif

      elem[i].nr = i+1;
      elem[i].group= 0;
      elem[i].mat = 1;

      defineElemNodes(anz, face, elem, i, eface, elemf[i].anz_f);
    }
  }
  else
  {
    printf ("\n%s opened\n",file);
    while(1)
    {
      length = frecord( handle, rec_str);
      if (rec_str[length] == (char)EOF) break;
  
      /* read data if a block starts */
      if(rec_str[0]=='(')
      {
        if ( (elem = (Elements *)malloc( (anz->e+1) * sizeof(Elements))) == NULL )
        {  printf("\n\n ERROR: malloc failed, elem\n\n") ; exit(1); }
        for(i=0; i<anz->e; i++)
        {
    	  elem[i].nr = i+1;
          elem[i].group= 0;
          elem[i].mat = 1;
          length = frecord( handle, rec_str);
          sscanf(rec_str, "%d",&anz_f);
          sscanf(&rec_str[2], "%d %d %d %d %d %d",&eface[0],&eface[1],&eface[2],&eface[3],&eface[4],&eface[5]);
          defineElemNodes(anz, face, elem, i, eface, anz_f);  
        }
        break;
      }
      else if(length>0) anz->e=atoi(rec_str);
    }
    fclose(handle);
  }
  anz->emax=anz->e;
  anz->emin=1;


  /* boundaries */
  /* located in "datin"/constant/polymesh/boundary */
  sprintf( file, "%s/constant/polyMesh/boundary", datin);
  handle = fopen (file, "r");
  if ( handle== NULL )
  {
    printf ("\nThe input file \"%s\" could not be opened.\n", file);
  }
  else printf ("\n%s opened\n",file);
  while(1)
  {
    length = frecord( handle, rec_str);
    if (rec_str[length] == (char)EOF) break;

    /* read data if a block starts */
    if(rec_str[0]=='(')
    {
      while(1)
      {
        length = frecord( handle, rec_str);
        if (rec_str[length] == (char)EOF) break;
        length = strfind(rec_str, ")");
        if ( length!=-1 ) break;

        length = strfind(rec_str, "{");
        if ( length!=-1 )
        {
          printf("\n   setname:%s\n",  buffer);
          setNr=getSetNr(buffer);
          if(setNr<0) setNr=pre_seta( buffer, "i", 0);
        }
        else sscanf(rec_str, "%s", buffer);

        length = strfind(rec_str, "nFaces");
        if ( length!=-1 )
        {
          sscanf(rec_str,"%*s %d", &nFaces);
          printf("   nFaces:%d\n", nFaces );
        }

        length = strfind(rec_str, "startFace");
        if ( length!=-1 )
        {
          sscanf(rec_str,"%*s %d", &startFace);
          printf("   startFace:%d\n", startFace );
        }

        /* store the faces in the set if the defining block was closed */
        length = strfind(rec_str, "}");
        if ( length!=-1 )
        {
          for(j=0; j<nFaces; j++)
          {
            if(face[j+startFace].nr>-1) 
            {
              seta(setNr,"f", j+startFace);
              i=seta(setNr,"j",0);
              if(i>-1)
              {
                set[setNr].elf[i].e=face[j+startFace].elem_nr;
                set[setNr].elf[i].f=face[j+startFace].nr;
              }
	    }
          }
        }
      }
    }
  }
  fclose(handle);

  /* scan over all project directories and search for results */
  sum_dir=0;
  dirp = opendir(datin);
  printf("seach for results in %s\n", datin );
  if (dirp != NULL) while ((dp = readdir(dirp)) != NULL)
  {

    /* if the dir-name starts with a number, its a result dir */
    if((dp->d_name[0]>=48)&&(dp->d_name[0]<=57))
    {
      if ( (resultdir2 = (char **)realloc( resultdir2, (sum_dir+2) * sizeof(char *))) == NULL )
        printf("\n\n ERROR: realloc failure\n\n" );
      if ( (resultdir2[sum_dir] = (char *)malloc( (MAX_LINE_LENGTH) * sizeof(char))) == NULL )
        printf("\n\n ERROR: realloc failure\n\n" );
      sprintf( resultdir2[sum_dir], "%s/%s", datin, dp->d_name );

      if ( (rsort = (Rsort *)realloc( rsort, (sum_dir+1) * sizeof(Rsort))) == NULL )
        printf("ERROR: realloc failed: Rsort\n\n" ); 
      rsort[sum_dir].r=atof(dp->d_name);
      rsort[sum_dir].i=sum_dir;
      sum_dir++;
    }
  }
  closedir(dirp);

  /* sort the resultdirs according to their name-value (its the time-step value) */
  qsort( rsort, sum_dir, sizeof(Rsort), (void *)compareRsort );
  if ( (resultdir = (char **)realloc( resultdir, (sum_dir+1) * sizeof(char *))) == NULL )
    printf("\n\n ERROR: realloc failure\n\n" );
  for( i=0; i<sum_dir; i++)
  {
    if ( (resultdir[i] = (char *)malloc( (MAX_LINE_LENGTH) * sizeof(char))) == NULL )
      printf("\n\n ERROR: realloc failure\n\n" );
    strcpy( resultdir[i], resultdir2[rsort[i].i]);
    printf(" found results in dir:%s\n",resultdir[i]);
  }


  if ( (vpn = (int *)malloc( (anz->nmax+1) * sizeof(int))) == NULL )
    printf("\n\n ERROR: malloc failure\n\n" );	               
  for(i=0; i<=anz->nmax; i++) vpn[i]=0;

  for(l=0; l<sum_dir; l++)
  {
    /* search for files in dir */
    dirp = opendir(resultdir[l]);
    if (dirp != NULL)
    {
     ds.p=ds.t=ds.u=ds.rho=-1;

     /* scan all files */
     while ((dp = readdir(dirp)) != NULL)
     {
      if ((strlen(dp->d_name) >0 )&&(dp->d_name[strlen(dp->d_name)-1]!='~')&&(dp->d_name[0]!='.'))
      {
        /* open the dir/file */
        sprintf(file, "%s/%s",resultdir[l],dp->d_name); 
        handle = fopen (file, "r");
        if ( handle== NULL )
        {
          printf ("\nThe input file \"%s\" could not be opened.\n", file);
          goto noFile;
        }
        else printf("file:%s opened\n", file);
 
        /* determine the type of data and create the entity names */
        ncomps=0;
        while(1)
        {
          length = frecord( handle, rec_str);
          // printf ("record: %s", rec_str); 
          if (rec_str[length] == (char)EOF) break;
  
          /* read data when the header starts */
          else if(rec_str[0]=='{')
          {
            while(1)
            {
              length = frecord( handle, rec_str);
              /* printf ("record: %s", rec_str); */
              if (rec_str[length] == (char)EOF) { printf("ERROR\n"); exit(1); }
              if (rec_str[0] == '}') goto foundEntities;

              sscanf(rec_str, "%s %s", dat[0], dat[1]);
              // printf ("dat:%s| dat:%s|\n",dat[0], dat[1]);

              if(compareStrings(dat[0], "class")>0)
              {
                if( compare(dat[1], "volScalarField", 14)==14) ncomps=1;
                else if( compare(dat[1], "volVectorField", 14)==14) ncomps=3;
                else { printf("Warning; class:%s not known\n",dat[1]); ncomps=0; }
              }
              if(compareStrings(dat[0], "object")>0)
              { strcpy(resultname,dat[1]); resultname[strlen(dat[1])]=0; }
            }
          }
        } foundEntities:;
        printf ("  resultname:%s entities:%d\n", resultname, ncomps);

        if(ncomps<1) goto foundNoEntities;
    
        /* define the lcases */
        anz->l++;
    
        /* store Datasets */
        if ( (lcase = (Datasets *)realloc(lcase, (anz->l+1) * sizeof(Datasets))) == NULL )
          printf("\n\n ERROR: malloc failed, lcase\n\n") ;

        i=strlen(resultdir[l]); while((i)&&(resultdir[l][i]!='/')) i--;
        lcase[anz->l].value=atof(&resultdir[l][i+1]);
        if(compareStrings(&resultdir[l][i+1],"0")==1) sprintf( lcase[anz->l].dataset_name,"%s", "initial");
        else sprintf( lcase[anz->l].dataset_name,"%s    ", "result");

        lcase[anz->l].irtype = 1;
        lcase[anz->l].npheader = 0;
        strcpy(lcase[anz->l].analysis_name,"");
        strcpy(lcase[anz->l].name, dp->d_name );
        strcpy(lcase[anz->l].dataset_text,"");
        lcase[anz->l].step_number=l+1;
        lcase[anz->l].analysis_type=1;

        printf ("lcase.name[%d]= %s value:%f\n", anz->l, lcase[anz->l].name, lcase[anz->l].value);
        if(lcase[anz->l].name[0]=='p') ds.p=anz->l;
        if(lcase[anz->l].name[0]=='U') ds.u=anz->l;
        if(lcase[anz->l].name[0]=='T') ds.t=anz->l;
        if(lcase[anz->l].name[0]=='r') ds.rho=anz->l;
        lcase[anz->l].ncomps  = ncomps;
    

        if ( (lcase[anz->l].nmax = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
          printf("\n\n ERROR: malloc failure\n\n" );
        if ( (lcase[anz->l].nmin = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
          printf("\n\n ERROR: malloc failure\n\n" );
        if ( (lcase[anz->l].max = (float *)malloc( lcase[anz->l].ncomps * sizeof(float))) == NULL )
          printf("\n\n ERROR: malloc failure\n\n" );
        if ( (lcase[anz->l].min = (float *)malloc( lcase[anz->l].ncomps * sizeof(float))) == NULL )
          printf("\n\n ERROR: malloc failure\n\n" );
        if ( (lcase[anz->l].dat = (float **)malloc( lcase[anz->l].ncomps * sizeof(float *))) == NULL )
          printf("\n\n ERROR: malloc failure\n\n" );
        if ( (lcase[anz->l].icname = (char **)malloc( lcase[anz->l].ncomps * sizeof(char *))) == NULL )
          printf("\n\n ERROR: malloc failure\n\n" );
        if ( (lcase[anz->l].compName = (char **)malloc( lcase[anz->l].ncomps * sizeof(char *))) == NULL )
          printf("\n\n ERROR: malloc failure\n\n" );
        for(i=0; i<lcase[anz->l].ncomps; i++)
        {
          if ( (lcase[anz->l].compName[i] = (char *)malloc( MAX_LINE_LENGTH * sizeof(char))) == NULL )
            printf("\n\n ERROR: malloc failed\n\n" );
          if ( (lcase[anz->l].dat[i] = (float *)malloc( (anz->nmax+1) * sizeof(float))) == NULL )
            printf("\n\n ERROR: malloc failed\n\n" );	               
          if ( (lcase[anz->l].icname[i] = (char *)malloc( MAX_LINE_LENGTH * sizeof(char))) == NULL )
             printf("\n\n ERROR: malloc failed\n\n" );
          lcase[anz->l].max[i]=-MAXFLOAT;
          lcase[anz->l].min[i]=MAXFLOAT;
        }
        if ( (lcase[anz->l].menu = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
          printf("\n\n ERROR: malloc failure\n\n" );
        if ( (lcase[anz->l].ictype = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
          printf("\n\n ERROR: malloc failure\n\n" );
        if ( (lcase[anz->l].icind1 = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
          printf("\n\n ERROR: malloc failure\n\n" );
        if ( (lcase[anz->l].icind2 = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
          printf("\n\n ERROR: malloc failure\n\n" );
        if ( (lcase[anz->l].iexist = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
          printf("\n\n ERROR: malloc failure\n\n" );
    
        for(i=0; i<lcase[anz->l].ncomps; i++)
        {
          sprintf(lcase[anz->l].compName[i],"%s%d ",resultname, i+1);
          lcase[anz->l].menu[i] = 1;
          if(ncomps==1)
	  {
            lcase[anz->l].ictype[i] = 1;
            lcase[anz->l].icind1[i] = 0;
	  }
          else if(ncomps==3)
          {
            lcase[anz->l].ictype[i]=2;
            lcase[anz->l].icind1[i]=i+1;
          }
          lcase[anz->l].icind2[i] = 0;
          lcase[anz->l].iexist[i] = 0;
          for(j=0; j<=anz->nmax; j++) lcase[anz->l].dat[i][j]=0.;
        }

        while(1)
        {
          length = frecord( handle, rec_str);
          if (rec_str[length] == (char)EOF) break;
      
          /* read "internalField" data */
          if ( strfind(rec_str, "internalField") !=-1 )
          { 
            sscanf(rec_str, "%*s %s", buffer);
            if(compareStrings(buffer,"uniform")>0)
            {
              formFlag=1;
              if     (ncomps==1) sscanf(rec_str, "%*s %*s %lf", &fval[0]); 
              else if(ncomps==3)
              {
                length = strfind(rec_str, "(");
                 sscanf(&rec_str[length+1], "%lf %lf %lf", &fval[0], &fval[1], &fval[2]);
	      }
              else { printf("ERROR1; class not known\n"); exit(1); }
		//printf("%s vals:%lf %lf %lf\n", &rec_str[length], fval[0], fval[1], fval[2]);
	    }
            else
            {
              formFlag=0;
              length = frecord( handle, rec_str);
              anz_dat=atoi(rec_str);
              if(anz->e!=anz_dat)
              { printf("ERROR: sum-data:%d do not match sum-elem:%d\n",anz_dat,anz->e); exit(1); }
              length = frecord( handle, rec_str);
	    }
    
            for(i=0; i<anz->e; i++)
            {
	      if(formFlag==0)
	      {
                length = frecord( handle, rec_str);
                if     (ncomps==1) sscanf(rec_str, "%lf", &fval[0]); 
                else if(ncomps==3)
                {
                  sscanf(&rec_str[1], "%lf %lf %lf", &fval[0], &fval[1], &fval[2]);
	        }
                else { printf("ERROR2; class not known (ncomps:%d)\n",ncomps); exit(1); }
	      }
              if (elem[i].type == 1) ipuf = 8;  /* HEXA8 */
              if (elem[i].type == 2) ipuf = 6;  /* PENTA6 */
              if (elem[i].type == 3) ipuf = 4;  /* TET4 */
              for(n=0; n<lcase[anz->l].ncomps; n++)
              {
                for(j=0; j<ipuf; j++)
                {
                  lcase[anz->l].dat[n][elem[i].nod[j]]+=fval[n];
                  if((vpn[0]==0)&&(n==0))
                  {
                    vpn[elem[i].nod[j]]++; /* count the hits */
	          }
                }
              }
            }
	    vpn[0]=1;
          }
	 
          /* read "boundaryField" data */
          if ( strfind(rec_str, "boundaryField") !=-1 )
          {
            length = frecord( handle, rec_str);
            while(1)
            {
              length = frecord( handle, rec_str);
              if (rec_str[length] == (char)EOF) break;
            
              length = strfind(rec_str, "{");
              if ( length!=-1 )
              {
                printf("\n   setname:%s\n",  buffer);

                setNr=getSetNr(buffer);
                if(setNr<0) { printf("ERROR; Set:%s in boundaryField not known\n",buffer); exit(1); }

		while(1)
                {
                  length = frecord( handle, rec_str);
                  if (rec_str[length] == (char)EOF) break;
                  if ( strfind(rec_str, "}")!=-1) break;


                  if ( strfind(rec_str, "fixedValue")!=-1)
                  {
		    do
                    {
                      length = frecord( handle, rec_str); if (rec_str[length] == (char)EOF) break;
                    }while(strfind(rec_str, "value")==-1);

                    /* reset all affected nodes to 0. */
                    for (j=0; j<set[setNr].anz_f; j++)
                    {
                      f=set[setNr].face[j];
                      if (face[f].type == 7) ipuf = 3;  /* TRI3 */
                      else if (face[f].type== 8) ipuf = 6;  /* TRI6  */  
                      else if (face[f].type == 9) ipuf = 4;  /* QU4 */
                      else if (face[f].type == 10) ipuf = 8;  /* QU4 */
                      else
                      {
                        printf("\n\n ERROR: face.type:%d not known\n",face[f].type);
                        exit(1);
                      }
                      for(n=0; n<lcase[anz->l].ncomps; n++)
                      {
                        for(k=0; k<ipuf; k++)
                        {
    	                  lcase[anz->l].dat[n][face[f].nod[k]]=0.;
                          vpn[face[f].nod[k]]=0;
                        }
                      }

    		    }

                    if ( strfind(rec_str, "nonuniform")!=-1)
                    {
		      f=-1;
                      do
		      {
  		        pos=0;
                        do
		        {
		          length=get_values_in_string( &rec_str[pos], &fval[0]);
                          if (length<1) break;
                          pos+=length;

                          /* go over all faces and look if the face is a member of the set */
                          do{
                            f++;
                            for (j=0; j<set[setNr].anz_elf; j++)
                            {
                              if(face[set[setNr].elf[j].e].indx[set[setNr].elf[j].f]==f) goto found_f;
			    }
                          }while(f<anz->f);

			found_f:;
                          if (face[f].type == 7) ipuf = 3;  /* TRI3 */
                          else if (face[f].type== 8) ipuf = 6;  /* TRI6  */  
                          else if (face[f].type == 9) ipuf = 4;  /* QU4 */
                          else if (face[f].type == 10) ipuf = 8;  /* QU4 */
                          else
                          {
                            printf("\n\n ERROR: face.type:%d not known\n",face[f].type);
                            exit(1);
                          }
                          for(n=0; n<lcase[anz->l].ncomps; n++)
                          {
                            for(j=0; j<ipuf; j++)
                            {
                              lcase[anz->l].dat[n][face[f].nod[j]]+=fval[n];
                              if(n==0) vpn[face[f].nod[j]]++; /* count the hits */
                            }
                          }
/*
			  if(ncomps==1) printf("face:%d e:%d val:%lf \n", f, face[f].elem_nr, fval[0]);
			  if(ncomps==3) printf("face:%d e:%d vals:%lf %lf %lf\n", f, face[f].elem_nr, fval[0], fval[1], fval[2]);
*/
		        }while(1);

                        if((length<0)||(f==-1)) { length = frecord( handle, rec_str); if (rec_str[length] == (char)EOF) break; }
                        else break;
		      }while(1); /* as long as a parenthesis is open */

		    }
                    else
		    {
                      if     (ncomps==1) sscanf(rec_str, "%*s %*s %lf", &fval[0]); 
                      else if(ncomps==3)
                      {
                        length = strfind(rec_str, "(");
                        sscanf(&rec_str[length+1], "%lf %lf %lf", &fval[0], &fval[1], &fval[2]);
  	              }
                      else { printf("ERROR3; class not known\n"); exit(1); }
  		      // printf("%s vals:%lf %lf %lf\n", &rec_str[length], fval[0], fval[1], fval[2]);
     
                      for(f=0; f<anz->f; f++)
                      {
                        foundFlag=0;
                        for (j=0; j<set[setNr].anz_elf; j++)
                        {
                          if(face[set[setNr].elf[j].e].indx[set[setNr].elf[j].f]==f) { foundFlag=1; break; }
	      	        }
                        if(foundFlag)
                        { 
                          if (face[f].type == 7) ipuf = 3;  /* TRI3 */
                          else if (face[f].type== 8) ipuf = 6;  /* TRI6  */  
                          else if (face[f].type == 9) ipuf = 4;  /* QU4 */
                          else if (face[f].type == 10) ipuf = 8;  /* QU4 */
                          else
                          {
                            printf("\n\n ERROR: face.type:%d not known\n",face[f].type);
                            exit(1);
                          }
                          for(n=0; n<lcase[anz->l].ncomps; n++)
                          {
                            for(j=0; j<ipuf; j++)
                            {
			      lcase[anz->l].dat[n][face[f].nod[j]]+=fval[n];
                              vpn[face[f].nod[j]]++; /* count the hits */
                            }
                          }
                        }
    		      }
		    }
		  }


                }
              }
              else sscanf(rec_str, "%s", buffer);
                
            }

          }
        }
 
        /* divide the nodal values by the amount of contributor elements */
        for(i=0; i<anz->n; i++)
        {
          for(n=0; n<lcase[anz->l].ncomps; n++)
          {
            lcase[anz->l].dat[n][node[i].nr]/=vpn[node[i].nr];
    
            if(lcase[anz->l].dat[n][node[i].nr] > lcase[anz->l].max[n])
            {
              lcase[anz->l].max[n]=lcase[anz->l].dat[n][node[i].nr];
              lcase[anz->l].nmax[n]=node[i].nr;
            }
            if(lcase[anz->l].dat[n][node[i].nr] < lcase[anz->l].min[n])
            {
              lcase[anz->l].min[n]=lcase[anz->l].dat[n][node[i].nr];
              lcase[anz->l].nmin[n]=node[i].nr;
            }
          }
        }

        foundNoEntities:;
        fclose(handle);
        noFile:;
      }
     }
     closedir(dirp);

     /* if data were found calc derived data like Ma etc. (based on pre-defined kappa ) */
     printf("p:%d u:%d t:%d rho:%d\n",ds.p, ds.u, ds.t, ds.rho) ;
     if((ds.p!=-1)&&(ds.u!=-1)&&(ds.t!=-1)&&(ds.rho!=-1))
     {
       anz->l++;
   
       /* store Datasets */
       if ( (lcase = (Datasets *)realloc(lcase, (anz->l+1) * sizeof(Datasets))) == NULL )
         printf("\n\n ERROR: malloc failed, lcase\n\n") ;

       lcase[anz->l].value=lcase[anz->l-1].value;
       sprintf( lcase[anz->l].dataset_name,"%s", lcase[anz->l-1].dataset_name);

       lcase[anz->l].irtype = 1;
       lcase[anz->l].npheader = 0;
       strcpy(lcase[anz->l].analysis_name,"");
       strcpy(lcase[anz->l].name, "derived");
       strcpy(lcase[anz->l].dataset_text,"");
       lcase[anz->l].step_number=lcase[anz->l-1].step_number;
       lcase[anz->l].analysis_type=1;

       printf ("lcase.name[%d]= %s value:%f\n", anz->l, lcase[anz->l].name, lcase[anz->l].value);
       lcase[anz->l].ncomps  = 3;

       if ( (lcase[anz->l].nmax = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
         printf("\n\n ERROR: malloc failure\n\n" );
       if ( (lcase[anz->l].nmin = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
         printf("\n\n ERROR: malloc failure\n\n" );
       if ( (lcase[anz->l].max = (float *)malloc( lcase[anz->l].ncomps * sizeof(float))) == NULL )
         printf("\n\n ERROR: malloc failure\n\n" );
       if ( (lcase[anz->l].min = (float *)malloc( lcase[anz->l].ncomps * sizeof(float))) == NULL )
         printf("\n\n ERROR: malloc failure\n\n" );
       if ( (lcase[anz->l].dat = (float **)malloc( lcase[anz->l].ncomps * sizeof(float *))) == NULL )
         printf("\n\n ERROR: malloc failure\n\n" );
       if ( (lcase[anz->l].icname = (char **)malloc( lcase[anz->l].ncomps * sizeof(char *))) == NULL )
         printf("\n\n ERROR: malloc failure\n\n" );
       if ( (lcase[anz->l].compName = (char **)malloc( lcase[anz->l].ncomps * sizeof(char *))) == NULL )
         printf("\n\n ERROR: malloc failure\n\n" );

       for(i=0; i<lcase[anz->l].ncomps; i++)
       {
         if ( (lcase[anz->l].compName[i] = (char *)malloc( MAX_LINE_LENGTH * sizeof(char))) == NULL )
           printf("\n\n ERROR: malloc failed\n\n" );
         if ( (lcase[anz->l].dat[i]      = (float *)malloc( (anz->nmax+1) * sizeof(float))) == NULL )
           printf("\n\n ERROR: malloc failed\n\n" );	               
         if ( (lcase[anz->l].icname[i]   = (char *)malloc( MAX_LINE_LENGTH * sizeof(char))) == NULL )
           printf("\n\n ERROR: malloc failed\n\n" );
         lcase[anz->l].max[i]=-MAXFLOAT;
         lcase[anz->l].min[i]=MAXFLOAT;
       }     
       if ( (lcase[anz->l].menu = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
         printf("\n\n ERROR: malloc failure\n\n" );
       if ( (lcase[anz->l].ictype = (int *)malloc(lcase[anz->l].ncomps * sizeof(int))) == NULL )
         printf("\n\n ERROR: malloc failure\n\n" );
       if ( (lcase[anz->l].icind1 = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
         printf("\n\n ERROR: malloc failure\n\n" );
       if ( (lcase[anz->l].icind2 = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
         printf("\n\n ERROR: malloc failure\n\n" );
       if ( (lcase[anz->l].iexist = (int *)malloc( lcase[anz->l].ncomps * sizeof(int))) == NULL )
         printf("\n\n ERROR: malloc failure\n\n" );
       
       for(i=0; i<lcase[anz->l].ncomps; i++)
       {
         lcase[anz->l].menu[i] = 1;
         lcase[anz->l].ictype[i] = 1;
         lcase[anz->l].icind1[i] = 0;
         lcase[anz->l].icind2[i] = 0;
         lcase[anz->l].iexist[i] = 0;
       }

       sprintf(lcase[anz->l].compName[0],"TT ");
       sprintf(lcase[anz->l].compName[1],"PT ");
       sprintf(lcase[anz->l].compName[2],"Ma ");

       for(i=0; i<anz->n; i++)
       {
         ps=lcase[ds.p].dat[0][node[i].nr];
         us=sqrt(lcase[ds.u].dat[0][node[i].nr]*lcase[ds.u].dat[0][node[i].nr]+lcase[ds.u].dat[1][node[i].nr]*lcase[ds.u].dat[1][node[i].nr]+lcase[ds.u].dat[2][node[i].nr]*lcase[ds.u].dat[2][node[i].nr]);
         ts=lcase[ds.t].dat[0][node[i].nr];
         rho=lcase[ds.rho].dat[0][node[i].nr];
         Rg=ps/rho/ts;
         cp=Rg/(1-1/kappa);
	 //printf("R:%f cp:%f\n", Rg, cp);
         tt=us*us*.5/cp + ts;
	 pt = ps/pow( (ts/tt), (kappa/(kappa-1)) );
         ma=us/sqrt(kappa*ps/rho);
         lcase[anz->l].dat[0][node[i].nr]=tt;
         lcase[anz->l].dat[1][node[i].nr]=pt;
         lcase[anz->l].dat[2][node[i].nr]=ma;

         for(j=0; j<lcase[anz->l].ncomps; j++)
	 {
           if(lcase[anz->l].dat[j][node[i].nr] > lcase[anz->l].max[j])
           {
             lcase[anz->l].max[j]=lcase[anz->l].dat[j][node[i].nr];
             lcase[anz->l].nmax[j]=node[i].nr;
           }
           if(lcase[anz->l].dat[j][node[i].nr] < lcase[anz->l].min[j])
           {
             lcase[anz->l].min[j]=lcase[anz->l].dat[j][node[i].nr];
             lcase[anz->l].nmin[j]=node[i].nr;
           }
         }
       }
     }

    }
  }
  free(vpn);
  anz->l++;
  elemChecker( anz->e, node, elem);


  /* free the temporary faces */
  free(face);
  anz->f=0;

  /* set all .loaded to 1 to indicate that the data are available */
  for (i=0; i<anz->l; i++) { lcase[i].fileptr=0; lcase[i].loaded=1; }

  /* delete all faces in sets */
  for(i=0; i<anz->sets; i++) set[i].anz_f=0;

 *nptr =  node; *eptr = elem; *lptr = lcase;
  return (1);
}
