/* --------------------------------------------------------------------  */
/*                          CALCULIX                                     */
/*                   - GRAPHICAL INTERFACE -                             */
/*                                                                       */
/*     A 3-dimensional pre- and post-processor for finite elements       */
/*              Copyright (C) 1996 Klaus Wittig                          */
/*                                                                       */
/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation; version 2 of           */
/*     the License.                                                      */
/*                                                                       */
/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */
/*                                                                       */
/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */
/* --------------------------------------------------------------------  */

#include <extUtil.h>

#define TEST     0

#define INI_FIELD_SIZE 100000


int readNG( char *datin, Summen *anz, Nodes **nptr, Elements **eptr, Datasets **lptr )
{
  FILE *handle;
  int i=0, elemflag=0;

  char rec_str[MAX_LINE_LENGTH];
  int  node_field_size, elem_field_size;
  int  e_nmax=1, e_nmin=1;
  int  length, sum,n;

  Nodes     *node=NULL;
  Elements  *elem=NULL;


  node_field_size=INI_FIELD_SIZE;
  do
  {
    if ( (node = (Nodes *)realloc( (Nodes *)node, (node_field_size+1) * sizeof(Nodes))) == NULL )
    {
      printf("WARNING: in readfrd() is INI_FIELD_SIZE:%d to large and is reduced\n", node_field_size );
      node_field_size/=2;
    }
    if(node_field_size<100)
    {
      printf("\n\n ERROR: not enough memory in readfrd()\n\n");
      exit(-1);
    }
  }while(!node);

  elem_field_size=INI_FIELD_SIZE;
  do
  {
    if((elem = (Elements *)realloc( (Elements *)elem, (elem_field_size+1) * sizeof(Elements))) == NULL )
    {
      printf("WARNING: in readfrd() is INI_FIELD_SIZE:%d to large and is reduced\n", elem_field_size );
      elem_field_size/=2;
    }
    if(elem_field_size<100)
    {
      printf("\n\n ERROR: not enough memory in readfrd()\n\n");
      exit(-1);
    }
  }while(!elem);

  anz->n=anz->e=anz->l=-1;

  /* Open the files and check to see that it was opened correctly */
  handle = fopen (datin, "r");
  if ( handle== NULL )  { printf ("ERROR: The input file \"%s\" could not be opened.\n\n", datin); return(-1); }
  else  printf (" file:%s opened\n", datin);

  printf (" reading NG format\n");

  strcpy(anz->model, datin);
  printf (" MODEL NAME:  %s", anz->model);

  anz->emax=-MAX_INTEGER;  anz->emin=MAX_INTEGER; elemflag=1;
  anz->nmax=-MAX_INTEGER;  anz->nmin= MAX_INTEGER;

  do
  {
    length = frecord( handle, rec_str);
    if(compare(rec_str,"surfaceelements",6)==6) break;
    if (rec_str[length] == (char)EOF) break;
  }while(1);

  /* jump the nr of elems */
  length = frecord( handle, rec_str);
  sscanf(rec_str, "%d", &sum);
  for(i=0; i<sum; i++)
  {
    do { length = frecord( handle, rec_str); }while(!length); 
    //printf("surf:%s\n",rec_str);
    if (rec_str[length] == (char)EOF) break;

    anz->e++;
    if (anz->e>=elem_field_size)
    {
      elem_field_size=anz->e+100;
      if((elem=(Elements *)realloc((Elements *)elem,(elem_field_size+1)*sizeof(Elements)))==NULL)
      {
        printf("\n\n ERROR: realloc failed, elem-index:%d\n\n", anz->e);
        return(-1);
      }
    }
    elem[anz->e].nr    = anz->e+1;
    elem[anz->e].group = 1;
    elem[anz->e].type  = 0;
    if (elem[anz->e].nr >  anz->emax)  anz->emax=elem[anz->e].nr;
    if (elem[anz->e].nr <  anz->emin)  anz->emin=elem[anz->e].nr;
    sscanf(rec_str, "%*s %*s %*s %*s %d", &n );
    if(n==3)
    {
      elem[anz->e].type  = 7;
      elem[anz->e].mat   = 7;
      sscanf(rec_str, "%*s %d %*s %*s %*s %d %d %d", &elem[anz->e].mat, &elem[anz->e].nod[0], &elem[anz->e].nod[1], &elem[anz->e].nod[2] );
    }
    else if(n==6)
    {
      elem[anz->e].type  = 8;
      elem[anz->e].mat   = 8;
      sscanf(rec_str, "%*s %d %*s %*s %*s %d %d %d %d %d %d", &elem[anz->e].mat, &elem[anz->e].nod[0], &elem[anz->e].nod[1], &elem[anz->e].nod[2], &elem[anz->e].nod[4], &elem[anz->e].nod[5], &elem[anz->e].nod[3] );
    }
    else printf("elem-type not known, nr of nodes:%d\n", n);
    anz->etype[elem[anz->e].type]++;
  }

  do
  {
    length = frecord( handle, rec_str);
    if(compare(rec_str,"volumeelements",6)==6) break;
    if (rec_str[length] == (char)EOF) break;
  }while(1);

  /* jump the nr of elems */
  length = frecord( handle, rec_str);
  sscanf(rec_str, "%d", &sum);
  for(i=0; i<sum; i++)
  {
    length = frecord( handle, rec_str);
    //printf("volu:%s\n",rec_str);
    if (rec_str[length] == (char)EOF) break;

    anz->e++;
    if (anz->e>=elem_field_size)
    {
      elem_field_size=anz->e+100;
      if((elem=(Elements *)realloc((Elements *)elem,(elem_field_size+1)*sizeof(Elements)))==NULL)
      {
        printf("\n\n ERROR: realloc failed, elem-index:%d\n\n", anz->e);
        return(-1);
      }
    }
    elem[anz->e].nr    = anz->e+1;
    elem[anz->e].group = 1;
    elem[anz->e].type  = 0;
    if (elem[anz->e].nr >  anz->emax)  anz->emax=elem[anz->e].nr;
    if (elem[anz->e].nr <  anz->emin)  anz->emin=elem[anz->e].nr;
    sscanf(rec_str, "%*s %d", &n );
    if(n==4)
    {
      elem[anz->e].mat   = 3;
      elem[anz->e].type  = 3;
      sscanf(rec_str, "%d %*s %d %d %d %d", &elem[anz->e].mat, &elem[anz->e].nod[0], &elem[anz->e].nod[2], &elem[anz->e].nod[1], &elem[anz->e].nod[3] );
    }
    else if(n==10)
    {
      elem[anz->e].mat   = 6;
      elem[anz->e].type  = 6;
      sscanf(rec_str, "%d %*s  %d %d %d %d %d %d %d %d %d %d", &elem[anz->e].mat, &elem[anz->e].nod[0], &elem[anz->e].nod[2], &elem[anz->e].nod[1],
      &elem[anz->e].nod[3], &elem[anz->e].nod[6], &elem[anz->e].nod[4], &elem[anz->e].nod[7], &elem[anz->e].nod[5], &elem[anz->e].nod[9], &elem[anz->e].nod[8] );
    }
    else printf("elem-type not known, nr of nodes:%d\n", n);
    anz->etype[elem[anz->e].type]++;
  }


  do
  {
    length = frecord( handle, rec_str);
    if(compare(rec_str,"edgesegments",6)==6) break;
    if (rec_str[length] == (char)EOF) break;
  }while(1);

  /* jump the nr of elems */
  length = frecord( handle, rec_str);
  sscanf(rec_str, "%d", &sum);
  for(i=0; i<sum; i++)
  {
    length = frecord( handle, rec_str);
    //printf("edge:%s\n",rec_str);
    if (rec_str[length] == (char)EOF) break;

    anz->e++;
    if (anz->e>=elem_field_size)
    {
      elem_field_size=anz->e+100;
      if((elem=(Elements *)realloc((Elements *)elem,(elem_field_size+1)*sizeof(Elements)))==NULL)
      {
        printf("\n\n ERROR: realloc failed, elem-index:%d\n\n", anz->e);
        return(-1);
      }
    }
    elem[anz->e].nr    = anz->e+1;
    elem[anz->e].group = 1;
    elem[anz->e].mat   = 1;
    elem[anz->e].type  = 11;
    if (elem[anz->e].nr >  anz->emax)  anz->emax=elem[anz->e].nr;
    if (elem[anz->e].nr <  anz->emin)  anz->emin=elem[anz->e].nr;
    sscanf(rec_str, "%d %*s %d %d", &elem[anz->e].mat, &elem[anz->e].nod[0], &elem[anz->e].nod[1] );
    anz->etype[elem[anz->e].type]++;
  }

  do
  {
    length = frecord( handle, rec_str);
    if (rec_str[length] == (char)EOF) break;
    if(compare(rec_str,"points",4)==4) break;
  }while(1);

  /* jump the nr of points */
  length = frecord( handle, rec_str);
  sscanf(rec_str, "%d", &sum);
  for(i=0; i<sum; i++)
  {
    length = frecord( handle, rec_str);
    if (rec_str[length] == (char)EOF) break;
    //printf("node:%s\n",rec_str);
    anz->n++;
    //node[anz->n].nr = anz->n;
    node[anz->n].nr = anz->n+1;
    if (node[anz->n].nr>=node_field_size)
    {
      node_field_size=node[anz->n].nr+100;
      if ( (node = (Nodes *)realloc((Nodes *)node, (node_field_size+1) * sizeof(Nodes))) == NULL )
      {
        printf("\n\n ERROR: realloc failed, nodenr:%d\n\n", node[anz->n].nr) ;
        return(-1);
      }
    }
    node[node[anz->n].nr].indx=anz->n;

    sscanf(rec_str,"%lf %lf %lf", &node[node[anz->n].nr].nx,&node[node[anz->n].nr].ny,&node[node[anz->n].nr].nz);
    if (node[anz->n].nr >  anz->nmax)  anz->nmax=node[anz->n].nr;
    if (node[anz->n].nr <  anz->nmin)  anz->nmin=node[anz->n].nr;
  }

  anz->n++;
  anz->e++;
  anz->l++;
  fclose(handle);

  node_field_size=anz->nmax+1;
  if((node =(Nodes *)realloc( (Nodes *)node, node_field_size*sizeof(Nodes)))==NULL)
    printf("\n\n ERROR: realloc failed\n\n") ;
  else
    printf ("\n %d nodes reallocated \n",anz->nmax);

  elem_field_size=anz->e+1;
  if ( (elem = (Elements *)realloc((Elements *)elem, elem_field_size * sizeof(Elements))) == NULL )
    printf("\n\n ERROR: in readfrd realloc failed\n\n") ;
  else
    printf ("\n %d elements reallocated \n", anz->e);

  if ( e_nmax > (anz->nmax) )
  {
    printf ("\nWARNING: element requestes a nodename higher than allocated\n\n");
    printf (" e_nmax=%d e_nmin=%d\n", e_nmax, e_nmin );
  }
  if ( e_nmin < 1 )
  {
    printf ("\nWARNING: element requestes a nodename lower than allocated\n\n");
    printf (" e_nmax=%d e_nmin=%d\n", e_nmax, e_nmin );
  }
  elemChecker( anz->e, node, elem);

  *nptr = node; *eptr = elem;
  return(1);
}

