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

/* TO DO
 TOL_AMBIG could be iterated to improve the yeld
*/

#include "uselibSNL.h"

#include <iostream>

#include <snlSurface.h>
#include <snlUtil.h>
#include <snlNurbsCommon.h>

#define ITER_TOL   1.e-10
#define COS_TOL    1.e-10
#define MAX_PASS   50
#define MAX_RETURNED_PNTS_PER_PNT 10
#define MIN_DIST   1.e-5
#define MAXVALUE   1.e99

#define     TEST            0  /* debugging */
#define     TEST2           0  /* debugging */
#define     PRINT_SNL       0  /* debugging */

#ifdef WIN32
extern "C" char  printFlag;                     /* printf 1:on 0:off */
#else
extern char  printFlag;                     /* printf 1:on 0:off */
#endif



typedef struct {
  int n;
  double paramU[MAX_RETURNED_PNTS_PER_PNT];
  double paramV[MAX_RETURNED_PNTS_PER_PNT];
  double dist[MAX_RETURNED_PNTS_PER_PNT];
}PntProj;


void vl_result( double *A, double *B, double *C )
/**********************************************************/
/*    Vektorbetrag: C =  Vektor(B)-Vektor(A) == Vector(AB)*/
/**********************************************************/
{
register int i;

for (i=0; i<3; i++)
         C[i]=B[i]-A[i];
}
void vl_prod( double *A, double *B, double *C )
/*********************************************************/
/*    Vektormultiplikation: C = A x B                    */
/*********************************************************/
{
      C[0]=A[1]*B[2]-A[2]*B[1];
      C[1]=A[2]*B[0]-A[0]*B[2];
      C[2]=A[0]*B[1]-A[1]*B[0];
}
double vl_norm( double *A, double *C )
/*********************************************************/
/*         Vector(B)=Vektor(A)/|A|                       */
/*********************************************************/
{
  int i;
  double    B;

      B=0.;
      for ( i=0; i<3; i++)
       B=B+A[i]*A[i];

      B = sqrt(B);
      if(B==0.) {C[0]=C[1]=C[2]=0.; return(B);} 
      for ( i=0; i<3; i++)
       C[i]=A[i]/B;
      return(B);
}
double vl_betrag(double *a)
/* ********************************************************* */
/*      laenge von Vektor a                                  */
/* ********************************************************* */
{
  register int i;
  double b;
  b=0.;

  for (i=0; i<3; i++)
    {
    b=b+a[i]*a[i];
    }
  b=sqrt(b);
  return(b);
}
void vl_scal( double *A, double *B, double *C )
/**********************************************************/
/* Vektormultiplikation: vektor(C) =  scalar(A)*Vektor(B) */
/**********************************************************/
{
  register int i;

  for (i=0; i<3; i++){
         C[i]= *A * B[i];
        }
}
void vl_add( double *A, double *B, double *C )
/**********************************************************/
/*    Vektoraddition: C =  Vektor(B)+Vektor(A)            */
/**********************************************************/
{
  register int i;

  for (i=0; i<3; i++){
         C[i]=B[i]+A[i];
        }
}


#ifdef __cplusplus
extern "C" {
#endif
   SumGeo    anzGeo[1];
   SumAsci   sumAsci[1];
#ifdef __cplusplus
}
#endif

int createBlendedNurbs(int nr, Points **pntpntr, Lines *line, Lcmb *lcmb, Gsur *surf, Nurbs **nurbspntr )
{
  int i,j,k,l,cl,n,p, nip, flag, np[4], sum, index;
  static double *px[4], *py[4], *pz[4];
  char name[MAX_LINE_LENGTH], buffer[MAX_LINE_LENGTH];
  char sname[MAX_LINE_LENGTH];

  int e;
  double p0[3],p1[3],p0p1[3],ep0p1[3],vl,dvl,p0pn[3],pn[3];

  Points *point;
  Nurbs *nurbs;

  nurbs=*nurbspntr;
  point=*pntpntr;

  // printf("surf:%s\n", surf[nr].name);
  if(surf[nr].nl!=4)
  {
    printf("ERROR: nr of edges must be 4 but is %d for surf:%s to create a blended nurbs\n\n",surf[nr].nl,surf[nr].name);
    return(-1);
  }
  /* create an array of points for each edge */
  for(j=0; j<surf[nr].nl; j++)
  {
    p=0;
    nip=0; 
    if(surf[nr].typ[j]=='l')
    {
      // printf("line:%s\n", line[surf[nr].l[j]].name);
      l=surf[nr].l[j];
      nip+=line[l].nip;

      if( (px[j] = (double *)realloc( (double *)px[j], (nip/3+1)*sizeof(double ) ))==NULL )
      { printf(" ERROR: realloc failure in createBlendedNurbs()\n\n");
        return(-1); }
      if( (py[j] = (double *)realloc( (double *)py[j], (nip/3+1)*sizeof(double ) ))==NULL )
      { printf(" ERROR: realloc failure in createBlendedNurbs()\n\n");
        return(-1); }
      if( (pz[j] = (double *)realloc( (double *)pz[j], (nip/3+1)*sizeof(double ) ))==NULL )
      { printf(" ERROR: realloc failure in createBlendedNurbs()\n\n");
        return(-1); }

      if(surf[nr].o[j]=='+')
      {
        n=0; flag=line[l].nip; 
        do
        {
          px[j][p]=line[l].ip[n++];
          py[j][p]=line[l].ip[n++];
          pz[j][p]=line[l].ip[n++];
          p++;
        }while(n<flag);
      }
      else
      {             
        n=line[l].nip; flag=0;
        while(n>flag)
        {
          pz[j][p]=line[l].ip[--n];
          py[j][p]=line[l].ip[--n];
          px[j][p]=line[l].ip[--n];
          p++;
        }
      }
    }
    else
    {
      // printf("lcmb:%s\n", lcmb[surf[nr].l[j]].name);
      if(surf[nr].o[j]=='+')
      {
       for(cl=0; cl<lcmb[surf[nr].l[j]].nl; cl++)
       {
        l=lcmb[surf[nr].l[j]].l[cl];

        nip+=line[l].nip;
        if( (px[j] = (double *)realloc( (double *)px[j], (nip/3+1)*sizeof(double ) ))==NULL )
        { printf(" ERROR: realloc failure in createBlendedNurbs()\n\n");
          return(-1); }
        if( (py[j] = (double *)realloc( (double *)py[j], (nip/3+1)*sizeof(double ) ))==NULL )
        { printf(" ERROR: realloc failure in createBlendedNurbs()\n\n");
          return(-1); }
        if( (pz[j] = (double *)realloc( (double *)pz[j], (nip/3+1)*sizeof(double ) ))==NULL )
        { printf(" ERROR: realloc failure in createBlendedNurbs()\n\n");
          return(-1); }

        if(lcmb[surf[nr].l[j]].o[cl]=='-') flag=-1;
        else flag=1;

        if(flag==1)
        {
                       
          if(cl==0)                     { n=0; flag=line[l].nip-3; }
          else if(cl==lcmb[surf[nr].l[j]].nl-1)    { n=0; flag=line[l].nip; }
          else                         { n=0; flag=line[l].nip-3; }
          
          do
          {
            px[j][p]=line[l].ip[n++];
            py[j][p]=line[l].ip[n++];
            pz[j][p]=line[l].ip[n++];
            p++;
          }while(n<flag);
        }
        else
        {
                      
          if(cl==0)                     { n=line[l].nip; flag=3; }
          else if(cl==lcmb[surf[nr].l[j]].nl-1)    { n=line[l].nip; flag=0; }
          else                         { n=line[l].nip; flag=3; }

          while(n>flag)
          {
            pz[j][p]=line[l].ip[--n];
            py[j][p]=line[l].ip[--n];
            px[j][p]=line[l].ip[--n];
            p++;
          }
        }
       }
      }
      else
      {
       for(cl=lcmb[surf[nr].l[j]].nl-1; cl>=0; cl--)
       {
        l=lcmb[surf[nr].l[j]].l[cl];

        nip+=line[l].nip;
        if( (px[j] = (double *)realloc( (double *)px[j], (nip/3+1)*sizeof(double ) ))==NULL )
        { printf(" ERROR: realloc failure in createBlendedNurbs()\n\n");
          return(-1); }
        if( (py[j] = (double *)realloc( (double *)py[j], (nip/3+1)*sizeof(double ) ))==NULL )
        { printf(" ERROR: realloc failure in createBlendedNurbs()\n\n");
          return(-1); }
        if( (pz[j] = (double *)realloc( (double *)pz[j], (nip/3+1)*sizeof(double ) ))==NULL )
        { printf(" ERROR: realloc failure in createBlendedNurbs()\n\n");
          return(-1); }

        if(lcmb[surf[nr].l[j]].o[cl]=='-') flag=1;
        else flag=-1;

        if(flag==1)
        {
                       
          if(cl==0)                     { n=0; flag=line[l].nip; }
          else if(cl==lcmb[surf[nr].l[j]].nl-1)    { n=0; flag=line[l].nip-3; }
          else                         { n=0; flag=line[l].nip-3; }
          
          do
          {
            px[j][p]=line[l].ip[n++];
            py[j][p]=line[l].ip[n++];
            pz[j][p]=line[l].ip[n++];
            p++;
          }while(n<flag);
        }
        else
        {
                      
          if(cl==0)                     { n=line[l].nip; flag=0; }
          else if(cl==lcmb[surf[nr].l[j]].nl-1)    { n=line[l].nip; flag=3; }
          else                         { n=line[l].nip; flag=3; }

          while(n>flag)
          {
            pz[j][p]=line[l].ip[--n];
            py[j][p]=line[l].ip[--n];
            px[j][p]=line[l].ip[--n];
            p++;
          }
        }
       }
      }
    }
    np[j]=p;
  }

  /* each edge must have at least 4 points */
  for (e=0; e<4; e++) if(np[e]<4)
  {
    if( (px[e] = (double *)realloc( (double *)px[e], (5)*sizeof(double ) ))==NULL )
    { printf(" ERROR: realloc failure in createBlendedNurbs()\n\n");
      return(-1); }
    if( (py[e] = (double *)realloc( (double *)py[e], (5)*sizeof(double ) ))==NULL )
    { printf(" ERROR: realloc failure in createBlendedNurbs()\n\n");
      return(-1); }
    if( (pz[e] = (double *)realloc( (double *)pz[e], (5)*sizeof(double ) ))==NULL )
    { printf(" ERROR: realloc failure in createBlendedNurbs()\n\n");
      return(-1); }
    p0[0]=px[e][np[e]-2];
    p0[1]=py[e][np[e]-2];
    p0[2]=pz[e][np[e]-2];
    p1[0]=px[e][np[e]-1];
    p1[1]=py[e][np[e]-1];
    p1[2]=pz[e][np[e]-1];
    vl_result(p0,p1,p0p1);
    vl_norm(p0p1,ep0p1); 
    vl=vl_betrag(p0p1);
    nip=4-np[e];
    dvl=vl/(nip+1);
    vl=0.;
    for(i=0; i<=nip; i++)
    {
      vl+=dvl;
      vl_scal(&vl,ep0p1, p0pn);
      vl_add(p0,p0pn,pn);
      // printf("PNT ! %f %f %f\n",  pn[0],pn[1],pn[2]);
      px[e][np[e]-1+i]=pn[0];
      py[e][np[e]-1+i]=pn[1];
      pz[e][np[e]-1+i]=pn[2];
    }
    np[e]+=nip;
  }
  /*
  for(j=0; j<4; j++)
  {
    printf("# pnts:%d\n", np[j]);
    for ( i= 0; i < np[j] ; i ++ ) printf(" PNT ! %f %f %f\n",  px[j][i],  py[j][i],  pz[j][i]); 
  }
  */

  // Create array of points to be passed to snlCurve constructor.
  // printf("  Create array of points\n"); 
  snlPoint* curvePointsEdge1 = new snlPoint [ np[0] ];
  snlPoint* curvePointsEdge2 = new snlPoint [ np[1] ];
  snlPoint* curvePointsEdge3 = new snlPoint [ np[2] ];
  snlPoint* curvePointsEdge4 = new snlPoint [ np[3] ];
  
  j=0; for ( i= 0; i < np[j] ; i ++ ) curvePointsEdge1 [ i ].components ( px[j][i], py[j][i], pz[j][i] );
  j=1; for ( i= 0; i < np[j] ; i ++ ) curvePointsEdge2 [ i ].components ( px[j][i], py[j][i], pz[j][i] );
  j=2; for ( i= np[j]-1; i >= 0 ; i -- ) curvePointsEdge3 [ np[j]-1-i ].components ( px[j][i], py[j][i], pz[j][i] );
  j=3; for ( i= np[j]-1; i >= 0 ; i -- ) curvePointsEdge4 [ np[j]-1-i ].components ( px[j][i], py[j][i], pz[j][i] );
  
  // Create curves. All curves must have the same degree. "degree" is an integer.
  int degree=3;
  snlCurve* curveEdgeU1 = new snlCurve ( 	curvePointsEdge1,
  								np[0],
  								snlCurve::SNL_GLOBAL_INTERP_CENTRIFUGAL,
  								degree 	);
  snlCurve* curveEdgeV2 = new snlCurve ( 	curvePointsEdge2,
  								np[1],
  								snlCurve::SNL_GLOBAL_INTERP_CENTRIFUGAL,
  								degree 	);
  snlCurve* curveEdgeU2 = new snlCurve ( 	curvePointsEdge3,
  								np[2],
  								snlCurve::SNL_GLOBAL_INTERP_CENTRIFUGAL,
  								degree 	);
  snlCurve* curveEdgeV1 = new snlCurve ( 	curvePointsEdge4,
  								np[3],
  								snlCurve::SNL_GLOBAL_INTERP_CENTRIFUGAL,
  								degree 	);
  
  // Create bilinear Coons patch. The orientation of the curves given to this
  // function is very important. There are two curves in the U direction
  // and two in the V direction. They should be oriented as follows:
  //
  //		V1
  //		-------->
  //		|         |
  //	U1	|         |  U2
  //		|         |
  //		V	 V
  //		 -------->
  //		V2
  //
  // The arrow heads are the end of the curve. If you want to reverse
  // the direction of the curve call it's reverseEvalDirection() function.
  
  // You will have to correspond curveEdge1 etc to one of the directions in
  // the following constructor.
  
  // printf(" Create bilinear Coons patch\n"); 
  snlSurface * surface = new snlSurface ( curveEdgeU1,curveEdgeU2,curveEdgeV1,curveEdgeV2 );
  
  // You now have a new NURBS Coons Patch.
  
  // Clean up allocated memory.
  
  delete[] curvePointsEdge1;
  delete[] curvePointsEdge2;
  delete[] curvePointsEdge3;
  delete[] curvePointsEdge4;
  
  delete curveEdgeU1;
  delete curveEdgeU2;
  delete curveEdgeV1;
  delete curveEdgeV2;


  // (5)   create the cgx nurbs
  // printf(" create blended nurbs\n"); 
  buffer[0]='S';
  buffer[1]=(char)NULL;
  getNewName( name, buffer );
  for (i=0; i<MAX_LINE_LENGTH; i++) sname[i]=name[i];
  if(printFlag) printf ("store NURS Nr:%d Name:%s\n", anzGeo->nurs+1, name);

  if ((nurbs = (Nurbs *)realloc( (Nurbs *)nurbs, (anzGeo->nurs+1)*sizeof(Nurbs)) ) == NULL )
  { printf("\n\nERROR: realloc failure in Nurs, nurbs:%s not installed\n\n", name); return(-1); }

  nr=anzGeo->nurs;
  sum=hashNurs( sumAsci, name, nr );
  anzGeo->nurs++;
  k=strlen(sname);
  if((nurbs[nr].name= (char *)malloc((k+1)*sizeof(char))) == NULL )
  { printf("ERROR: malloc failed\n\n" ); return(-1); }
  for (i=0; i<=k; i++) nurbs[nr].name[i]=name[i];

  /* the nurbs fits exactly, no trimming is necessary */
  nurbs[nr].trimFlag=0;
  nurbs[nr].u_exp=surface->degreeU();
  nurbs[nr].v_exp=surface->degreeV();
  nurbs[nr].u_npnt=surface->sizeU();
  nurbs[nr].v_npnt=surface->sizeV();
  nurbs[nr].u_nknt=1+nurbs[nr].u_npnt+ nurbs[nr].u_exp;
  nurbs[nr].v_nknt=1+nurbs[nr].v_npnt+ nurbs[nr].v_exp;
  nurbs[nr].v_stride=4;
  nurbs[nr].u_stride=4* nurbs[nr].v_npnt;

  nurbs[nr].np=NULL;
  nurbs[nr].uv=NULL;
  nurbs[nr].xyz=NULL;

  nurbs[nr].nc=NULL;
  nurbs[nr].vmax=NULL;
  nurbs[nr].umax=NULL;
  nurbs[nr].ctlarray=NULL; 

  nurbs[nr].patches=0;
  nurbs[nr].endFlag=1;       
  nurbs[nr].type=GL_MAP2_VERTEX_4;       
  nurbs[nr].Nurb = (GLUnurbsObj *)gluNewNurbsRenderer();

  // knots
  if ( (nurbs[nr].uknt = (float *)malloc((nurbs[nr].u_nknt+1) * sizeof(float))) == NULL )
    printf("\n\n ERROR: malloc failed uknt\n\n");
  if ( (nurbs[nr].vknt = (float *)malloc((nurbs[nr].v_nknt+1) * sizeof(float))) == NULL )
    printf("\n\n ERROR: malloc failed vknt\n\n");
  for(i=0; i<nurbs[nr].u_nknt; i++) { nurbs[nr].uknt[i]=surface->knotsU()[i]; }
  for(i=0; i<nurbs[nr].v_nknt; i++) { nurbs[nr].vknt[i]=surface->knotsV()[i]; }

  // discrete controll points and weights
  if ( (nurbs[nr].ctlpnt = (int **)malloc(  (nurbs[nr].u_npnt+1) * sizeof(int *))) == NULL )
    printf("\n\n ERROR: malloc failed ctlpnt\n\n");
  for (i=0; i<nurbs[nr].u_npnt; i++)
  {
    if ( (nurbs[nr].ctlpnt[i] = (int *)malloc(  (nurbs[nr].v_npnt+1) * sizeof( int ))) == NULL )
      printf("\n\n ERROR: malloc failed ctlpnt[i]\n\n");
  }
  if ( (nurbs[nr].weight = (GLfloat **)malloc(  (nurbs[nr].u_npnt+1) * sizeof(GLfloat *))) == NULL )
    printf("\n\n ERROR: malloc failed weight\n\n");
  for (i=0; i<nurbs[nr].u_npnt; i++)
  {
    if ( (nurbs[nr].weight[i] = (GLfloat *)malloc(  (nurbs[nr].v_npnt+1) * sizeof(GLfloat))) == NULL )
      printf("\n\n ERROR: malloc failed weight[i]\n\n");
  }

  index=0;
  for (i=0; i<nurbs[nr].u_npnt; i++)
  {
    for (j=0; j<nurbs[nr].v_npnt; j++)
    {
      buffer[0]='p';
      buffer[1]=(char)NULL;
      getNewName( name, buffer );
      nurbs[nr].ctlpnt[i][j]  = pnt( name, surface->controlPoints()[ index ].x(), surface->controlPoints()[ index ].y(), surface->controlPoints()[ index ].z(), 0);
      point=*pntpntr;  /* update pointer */
#if TEST
      printf("p:%s\n", point[nurbs[nr].ctlpnt[i][j]].name);
#endif
      nurbs[nr].weight[i][j]  = surface->controlPoints()[ index ].w();
      index++;
    }
  }

#if TEST
  printf("out: degUV %d %d sizeUV %u %u\n", surface->degreeU(), surface->degreeV(), surface->sizeU(), surface->sizeV());
  printf("knotu:%d\n", nurbs[nr].u_nknt);
  printf("knotv:%d\n", nurbs[nr].v_nknt);
  for(i=0; i<nurbs[nr].u_nknt; i++) { printf("knu:%f\n", nurbs[nr].uknt[i]); }
  for(i=0; i<nurbs[nr].v_nknt; i++) { printf("knv:%f\n", nurbs[nr].vknt[i]); }
  index=0;
  for (i=0; i<nurbs[nr].u_npnt; i++)
  {
    for (j=0; j<nurbs[nr].v_npnt; j++)
    {
      printf("pnt: %f %f %f %f \n", surface->controlPoints()[ index ].x(),surface->controlPoints()[ index ].y(),surface->controlPoints()[ index ].z(),surface->controlPoints()[ index ].w());
      index++;
    }
  }
#endif

  delete surface;   // Release surface object.
  *nurbspntr=nurbs;
  *pntpntr=point;
  return(nr);
}



// uses double reduceDegree( int dir, unsigned numDeg, double tolerance ) 
// returns the error during reduction. "dir" is the direction you wish to degree
// reduce, it is one of the constants:
//
//        snlSurface::SNL_U_DIR
//        snlSurface::SNL_V_DIR
//
// for the U and V directions respectively. "numDeg" is the number of degrees you w
// ish the surface to be reduced by and "tolerance" is the maximum error allowed du
// ring processing. If tolerance is exceeded then the surface is restored to origin
// al. If you don't want a tolerance and will take any error then use a value of 0.
// This function will increase the size of the data describing the surface consider
// ably. I am yet to implement a function that remedies that, but in the mean time
// it will work just the same.
double repairNurbs( Nurbs *nurbs, int nr, int deg, int dir )
{
  int i,j,index=0;
  int numberOfControlPoints =nurbs[nr].u_npnt*nurbs[nr].v_npnt;
  snlCtrlPoint* controlPoints = new snlCtrlPoint [ numberOfControlPoints ];
  double error;


  // (1)      Create an array of control points from your cgx control points:
  for (i=0; i<nurbs[nr].u_npnt; i++)
  {
    for (j=0; j<nurbs[nr].v_npnt; j++)
    {
      controlPoints [ index++ ].components (
	  nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+0],
	  nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+1],
	  nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+2],
	  nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+3] );
    }
  }

  // (2)      Create an array of knots for each knot vector in your cgx surface:
  knot* knotsU = new knot [ nurbs[nr].u_nknt ];
  knot* knotsV = new knot [ nurbs[nr].v_nknt ];

  for(i=0; i<nurbs[nr].u_nknt; i++)  knotsU[i]=nurbs[nr].uknt[i];
  for(i=0; i<nurbs[nr].v_nknt; i++)  knotsV[i]=nurbs[nr].vknt[i];

  // (3)      Create an snlSurface:
  // You need to supply degreeU, degreeV, sizeU and sizeV.
  snlSurface* surface = new snlSurface ( nurbs[nr].u_exp, nurbs[nr].v_exp, nurbs[nr].u_npnt, nurbs[nr].v_npnt, controlPoints, knotsU, knotsV );

#if TEST
  // print the nurbs-parameter before reduction 
  printf("in: degUV %d %d sizeUV %u %u\n", surface->degreeU(), surface->degreeV(), surface->sizeU(), surface->sizeV());
  for(i=0; i<nurbs[nr].u_nknt; i++) {  printf("knu:%f\n", nurbs[nr].uknt[i]); }
  for(i=0; i<nurbs[nr].v_nknt; i++) {  printf("knv:%f\n", nurbs[nr].vknt[i]); }
  index=0;
  for (i=0; i<nurbs[nr].u_npnt; i++)
  {
    for (j=0; j<nurbs[nr].v_npnt; j++)
    {
      printf("pnt: %f %f %f %f \n", controlPoints [ index ].x(),controlPoints [ index ].y(),controlPoints [ index ].z(),controlPoints [ index ].w());
      index++;
    }
  }
#endif

  // (4) reduce degree
  if(dir==0) error=surface -> reduceDegree ( snlSurface::SNL_U_DIR, deg, 0.0 );
  else error=surface -> reduceDegree ( snlSurface::SNL_V_DIR, deg, 0.0 );

  // (5)   change the cgx nurbs
  nurbs[nr].u_exp=surface->degreeU();
  nurbs[nr].v_exp=surface->degreeV();
  nurbs[nr].u_npnt=surface->sizeU();
  nurbs[nr].v_npnt=surface->sizeV();
  nurbs[nr].u_nknt=1+nurbs[nr].u_npnt+ nurbs[nr].u_exp;
  nurbs[nr].v_nknt=1+nurbs[nr].v_npnt+ nurbs[nr].v_exp;
  nurbs[nr].u_stride=4* nurbs[nr].v_npnt;

  // knots
  if ( (nurbs[nr].uknt = (float *)realloc(  (float *)nurbs[nr].uknt, (nurbs[nr].u_nknt+1) * sizeof(float))) == NULL )
    printf("\n\n ERROR: realloc failed uknt\n\n");
  if ( (nurbs[nr].vknt = (float *)realloc(  (float *)nurbs[nr].vknt, (nurbs[nr].v_nknt+1) * sizeof(float))) == NULL )
    printf("\n\n ERROR: realloc failed vknt\n\n");
  for(i=0; i<nurbs[nr].u_nknt; i++) { nurbs[nr].uknt[i]=surface->knotsU()[i]; }
  for(i=0; i<nurbs[nr].v_nknt; i++) { nurbs[nr].vknt[i]=surface->knotsV()[i]; }

  // controll-array, discrete controll points and weights
  if( (nurbs[nr].ctlarray = (GLfloat *)realloc( (GLfloat *)nurbs[nr].ctlarray, (nurbs[nr].u_npnt*nurbs[nr].v_npnt*nurbs[nr].v_stride+5)*sizeof(GLfloat) )) == NULL )
  { printf(" ERROR: realloc failure in repairNurbs(), nurbs:%s can not be shaped\n\n", nurbs[nr].name);
      return(-1); }
  index=0;
  for (i=0; i<nurbs[nr].u_npnt; i++)
  {
    for (j=0; j<nurbs[nr].v_npnt; j++)
    {
      nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+0]=surface->controlPoints()[ index ].x();
      nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+1]=surface->controlPoints()[ index ].y();
      nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+2]=surface->controlPoints()[ index ].z();
      nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+3]=surface->controlPoints()[ index ].w();
      index++;
    }
  }

#if TEST
  // print 
  printf("out: degUV %d %d sizeUV %u %u\n", surface->degreeU(), surface->degreeV(), surface->sizeU(), surface->sizeV());
  printf("knotu:%d\n", nurbs[nr].u_nknt);
  printf("knotv:%d\n", nurbs[nr].v_nknt);
  for(i=0; i<nurbs[nr].u_nknt; i++) { printf("knu:%f\n", nurbs[nr].uknt[i]); }
  for(i=0; i<nurbs[nr].v_nknt; i++) { printf("knv:%f\n", nurbs[nr].vknt[i]); }
  index=0;
  for (i=0; i<nurbs[nr].u_npnt; i++)
  {
    for (j=0; j<nurbs[nr].v_npnt; j++)
    {
      printf("pnt: %f %f %f %f \n", surface->controlPoints()[ index ].x(),surface->controlPoints()[ index ].y(),surface->controlPoints()[ index ].z(),surface->controlPoints()[ index ].w());
      index++;
    }
  }
#endif

  delete surface;   // Release surface object.

  // The rest of the objects are deleted by libSNL.

  return(error);
}


int evalNurbs( Nurbs *nurbs, int nr, int sum_p, double *pnt_u, double *pnt_v, Points *pnt)
{ 
  int i,j,index=0;
  snlPoint snlPnt;
  int numberOfControlPoints =nurbs[nr].u_npnt*nurbs[nr].v_npnt;
  snlCtrlPoint* controlPoints = new snlCtrlPoint [ numberOfControlPoints ];


  // (1)      Create an array of control points from your cgx control points:
  for (i=0; i<nurbs[nr].u_npnt; i++)
  {
    for (j=0; j<nurbs[nr].v_npnt; j++)
    {
      controlPoints [ index++ ].components (
	  nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+0],
	  nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+1],
	  nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+2],
	  nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+3] );
    }
  }

  // (2)      Create an array of knots for each knot vector in your cgx surface:
  knot* knotsU = new knot [ nurbs[nr].u_nknt ];
  knot* knotsV = new knot [ nurbs[nr].v_nknt ];

  for(i=0; i<nurbs[nr].u_nknt; i++)  knotsU[i]=nurbs[nr].uknt[i];
  for(i=0; i<nurbs[nr].v_nknt; i++)  knotsV[i]=nurbs[nr].vknt[i];

  // (3)      Create an snlSurface:
  // You need to supply degreeU, degreeV, sizeU and sizeV.
  snlSurface* surface = new snlSurface ( nurbs[nr].u_exp, nurbs[nr].v_exp, nurbs[nr].u_npnt, nurbs[nr].v_npnt, controlPoints, knotsU, knotsV );

  // (4) get coordinates
  for(i=0; i<sum_p; i++)
  {
    // negative values lead to a segfault
    if(pnt_u[i]<0.) pnt_u[i]=0.;
    if(pnt_v[i]<0.) pnt_v[i]=0.;
    snlPnt= surface -> eval ( pnt_u[i], pnt_v[i] );
    pnt[i].px=snlPnt.x();
    pnt[i].py=snlPnt.y();
    pnt[i].pz=snlPnt.z();
  }

  delete surface;   // Release surface object.

  // The rest of the objects are deleted by libSNL.
  return(1);
}



void projSurfToNurbs( Nurbs *nurbs, int nr, Gsur *surf, int snr, Nodes *node )
{
  int i,j,n,k,index=0;
  int sum_p=0, sum_inverted=0, returnedPntsPerPnt;

  double convergTol,  normTol;
  int maxPass; 
  snlSurfLocn* inverted;
  int *isort=NULL;
  int numberOfControlPoints =nurbs[nr].u_npnt*nurbs[nr].v_npnt;
  snlCtrlPoint* controlPoints = new snlCtrlPoint [ numberOfControlPoints ];

  normTol=COS_TOL;
  convergTol=ITER_TOL;
  maxPass=MAX_PASS;

  if(node) sum_p=surf[snr].nn;
  else     sum_p=surf[snr].npgn;

  if(!sum_p) return;

  // (1)      Create an array of control points from your cgx control points:
  for (i=0; i<nurbs[nr].u_npnt; i++)
  {
    for (j=0; j<nurbs[nr].v_npnt; j++)
    {
      controlPoints [ index++ ].components (
	  nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+0],
	  nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+1],
	  nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+2],
	  nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+3] );
    }
  }

  // (2)      Create an array of knots for each knot vector in your cgx surface:
  knot* knotsU = new knot [ nurbs[nr].u_nknt ];
  knot* knotsV = new knot [ nurbs[nr].v_nknt ];

  for(i=0; i<nurbs[nr].u_nknt; i++)  knotsU[i]=nurbs[nr].uknt[i];
  for(i=0; i<nurbs[nr].v_nknt; i++)  knotsV[i]=nurbs[nr].vknt[i];

  // (3)      Create an snlSurface:
  // You need to supply degreeU, degreeV, sizeU and sizeV.
  snlSurface* surface = new snlSurface ( nurbs[nr].u_exp, nurbs[nr].v_exp, nurbs[nr].u_npnt, nurbs[nr].v_npnt,
                                               controlPoints, knotsU, knotsV );
  // (4)      Project some points to the surface:
  // Create an array of snlPoints to project onto the surface.
  snlPoint* toProject = new snlPoint [ sum_p ]; // You supply numPointsToProject.

  // Fill the snlPoint array with data from cgx
  if(node)
  {
    for(i=0; i<sum_p; i++)
    {
      toProject[i].components (node[surf[snr].nod[i]].nx, node[surf[snr].nod[i]].ny, node[surf[snr].nod[i]].nz);
    }
  }
  else
  {
    sum_p=n=0;
    while((surf[snr].npgn-n))
    {
      /* create temporary nodes */
      n++;
      j=(int)surf[snr].pgn[n++];
      n+=3;
      for(k=0; k<j; k++)
      {
        toProject[sum_p].components (surf[snr].pgn[n], surf[snr].pgn[n+1], surf[snr].pgn[n+2]);
        sum_p++;
        n+=3;
      }
    }
  }

  // Project points to surface.
  returnedPntsPerPnt =MAX_RETURNED_PNTS_PER_PNT/5;
  isort = new int [ sum_p ];

 finerProjection:;
  inverted =surface -> fastProject( toProject, sum_p, &sum_inverted, convergTol, normTol, maxPass, 1, returnedPntsPerPnt );

  //printf("sum_inverted:%d\n", sum_inverted);
  //for(i=0; i<sum_inverted; i++) printf("%d node:%d dist:%f\n", inverted[i].origPtIndex, surf[snr].nod[inverted[i].origPtIndex],  inverted[i].dist);

  /* store the index of the clostest entity */
  for(i=0; i<sum_p; i++) isort[i]=-1;
  for(i=0; i<sum_inverted; i++)
  {
    n=inverted[i].origPtIndex;
    if(isort[n]==-1) isort[n]=i;
    else if(inverted[i].dist<inverted[isort[n]].dist) isort[n]=i;
  }
  //printf("sum:%d\n", sum_p);
  //for(i=0; i<sum_p; i++) printf("%d node:%d dist:%f\n", isort[i], surf[snr].nod[inverted[isort[i]].origPtIndex],  inverted[isort[i]].dist);

  if(node)
  {
    for(i=0; i<sum_p; i++)
    {
      // check if the distance is "close"
      if((inverted[isort[i]].dist>MIN_DIST)&&(returnedPntsPerPnt<MAX_RETURNED_PNTS_PER_PNT))
      { returnedPntsPerPnt=MAX_RETURNED_PNTS_PER_PNT; delete[] inverted; goto finerProjection; }

      node[surf[snr].nod[i]].nx = inverted[isort[i]].pt.x();
      node[surf[snr].nod[i]].ny = inverted[isort[i]].pt.y();
      node[surf[snr].nod[i]].nz = inverted[isort[i]].pt.z();
    }
  }
  else
  {
    i=n=0;
    while((surf[snr].npgn-n))
    {
      /* create temporary nodes */
      n++;
      j=(int)surf[snr].pgn[n++];
      n+=3;
      for(k=0; k<j; k++)
      {
        // check if the distance is "close"
        if((inverted[isort[i]].dist>MIN_DIST)&&(returnedPntsPerPnt<MAX_RETURNED_PNTS_PER_PNT))
        { returnedPntsPerPnt=MAX_RETURNED_PNTS_PER_PNT; delete[] inverted; goto finerProjection; }

        surf[snr].pgn[n]   = inverted [ isort[i]].pt.x();
        surf[snr].pgn[n+1] = inverted [ isort[i]].pt.y();
        surf[snr].pgn[n+2] = inverted [ isort[i]].pt.z();
        i++;
        n+=3;
      }
    }
  }
  delete[] isort; 
  delete surface;  // Release surface object.
  delete[] inverted;  // Release snlVertex array returned from projection function.
  delete[] toProject;  // Release points that were projected onto surface.

  // The rest of the objects are deleted by libSNL.
}


void projPntsToNurbs( Nurbs *nurbs, int nr, Sets *set, int setNr, Points *pnt)
{
  int i,j,n,index=0;

  double convergTol,  normTol;
  int maxPass,sum_p, sum_inverted=0, returnedPntsPerPnt; 
  snlSurfLocn* inverted;
  int numberOfControlPoints =nurbs[nr].u_npnt*nurbs[nr].v_npnt;
  snlCtrlPoint* controlPoints = new snlCtrlPoint [ numberOfControlPoints ];

  normTol=COS_TOL;
  convergTol=ITER_TOL;
  maxPass=MAX_PASS;

 
  // (1)      Create an array of control points from your cgx control points:
  for (i=0; i<nurbs[nr].u_npnt; i++)
  {
    for (j=0; j<nurbs[nr].v_npnt; j++)
    {
      controlPoints [ index++ ].components (
	  nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+0],
	  nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+1],
	  nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+2],
	  nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+3] );
    }
  }

  // (2)      Create an array of knots for each knot vector in your cgx surface:

  knot* knotsU = new knot [ nurbs[nr].u_nknt ];
  knot* knotsV = new knot [ nurbs[nr].v_nknt ];

  for(i=0; i<nurbs[nr].u_nknt; i++)  knotsU[i]=nurbs[nr].uknt[i];
  for(i=0; i<nurbs[nr].v_nknt; i++)  knotsV[i]=nurbs[nr].vknt[i];

  // (3)      Create an snlSurface:
  // You need to supply degreeU, degreeV, sizeU and sizeV.
  snlSurface* surface = new snlSurface ( nurbs[nr].u_exp, nurbs[nr].v_exp, nurbs[nr].u_npnt, nurbs[nr].v_npnt,
                                               controlPoints, knotsU, knotsV );
  // (4)      Project some points to the surface:
  // Create an array of snlPoints to project onto the surface.
  snlPoint* toProject = new snlPoint [ set[setNr].anz_p ]; // You supply numPointsToProject.

  // Fill the snlPoint array with data from cgx
  for(i=0; i<set[setNr].anz_p; i++)
  {
    toProject[i].components (pnt[set[setNr].pnt[i]].px, pnt[set[setNr].pnt[i]].py, pnt[set[setNr].pnt[i]].pz);
#if TEST
    printf("1xyz: %f %f %f\n",pnt[set[setNr].pnt[i]].px, pnt[set[setNr].pnt[i]].py, pnt[set[setNr].pnt[i]].pz);
#endif
  }

  // Project points to surface.
  sum_p=set[setNr].anz_p;
  returnedPntsPerPnt =MAX_RETURNED_PNTS_PER_PNT/5;
  int* isort = new int [ sum_p ];

 finerProjection:;
  inverted =surface -> fastProject( toProject, sum_p, &sum_inverted, convergTol, normTol, maxPass, 1, returnedPntsPerPnt );

  //printf("sum_inverted:%d\n", sum_inverted);
  //for(i=0; i<sum_inverted; i++) printf("%d node:%d dist:%f\n", inverted[i].origPtIndex, surf[snr].nod[inverted[i].origPtIndex],  inverted[i].dist);

  /* store the index of the clostest entity */
  for(i=0; i<sum_p; i++) isort[i]=-1;
  for(i=0; i<sum_inverted; i++)
  {
    n=inverted[i].origPtIndex;
    if(isort[n]==-1) isort[n]=i;
    else if(inverted[i].dist<inverted[isort[n]].dist) isort[n]=i;
  }
  //printf("sum:%d\n", sum_p);
  //for(i=0; i<sum_p; i++) printf("%d node:%d dist:%f\n", isort[i], surf[snr].nod[inverted[isort[i]].origPtIndex],  inverted[isort[i]].dist);

  for(i=0; i<set[setNr].anz_p; i++)
  {
    // check if the distance is "close"
    if((inverted[isort[i]].dist>MIN_DIST)&&(returnedPntsPerPnt<MAX_RETURNED_PNTS_PER_PNT))
    { returnedPntsPerPnt=MAX_RETURNED_PNTS_PER_PNT; delete[] inverted; goto finerProjection; }

    pnt[set[setNr].pnt[i]].px = inverted [ isort[i]].pt.x();
    pnt[set[setNr].pnt[i]].py = inverted [ isort[i]].pt.y();
    pnt[set[setNr].pnt[i]].pz = inverted [ isort[i]].pt.z();

#if TEST
    printf("2xyz: %f %f %f\n",pnt[set[setNr].pnt[i]].px, pnt[set[setNr].pnt[i]].py, pnt[set[setNr].pnt[i]].pz);
#endif
  }

  delete[]  isort; 
  delete surface;  // Release surface object.
  delete[] inverted;  // Release snlVertex array returned from projection function.
  delete[] toProject;  // Release points that were projected onto surface.

  // The rest of the objects are deleted by libSNL.
}



/* returns the type ( -1error, 0plate, 1cyl, 2torus, 3ball, 4half-ball-bottom, 5half-ball-top) */
/* ( half-ball-top: closed at top)  */
/* axis=1 ->u (ambig point has either vmax or vmin) */
/* axis=2 ->v (ambig point has either umax or umin) */
/* axis=3 ->uv (ambig point has either umax or umin and either vmax or vmin) */
/* ambig[0]==1 u-ambig at vmin  axis=+2 */
/* ambig[1]==1 u-ambig at vmax  axis=+2 */
/* ambig[2]==1 v-ambig at umin  axis=1 */
/* ambig[3]==1 v-ambig at umax  axis=1 */
/* collapsed[0]==1  u-collapsed at vmin  axis=2 */
/* collapsed[1]==1  u-collapsed at vmax  axis=2 */
/* collapsed[2]==1  v-collapsed at umin  axis=1 */
/* collapsed[3]==1  v-collapsed at umax  axis=1 */
int getNurbsType( Nurbs *nurbs, int nr, int *axis)
{
  int i,j,k,n;
  int collapsed[4], ambig[4];
  
  double diff1;

   *axis=0;
  for(k=0; k<4; k++)  collapsed[k]=ambig[k]=0;


  /* check if one edge is ambiguous */
  /* this is the case if control points meet each other */

  /* check v at umin and umax */
  for(k=2; k<4; k++)
  {
    if(k==2) i=0; else i=nurbs[nr].u_npnt-1;
#if TEST2
      printf("check at umin and umax\n");
      j=0;
      printf("x:%f ", nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+0]);
      printf("y:%f ", nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+1]);
      printf("z:%f ", nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+2]);
      printf("w:%f ", nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+3]);
      printf(" i:%d j:%d\n", i, j);
#endif

    /* check the first and last only */
    /* if they are equal its an ambiguous edge */ 
    /* difference to the first point */
    ambig[k]=1;
    j=nurbs[nr].v_npnt-1;
#if TEST2
      printf("x:%f ", nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+0]);
      printf("y:%f ", nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+1]);
      printf("z:%f ", nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+2]);
      printf("w:%f ", nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+3]);
      printf(" i:%d j:%d\n", i, j);
#endif
    for(n=0; n<4; n++)
    {
      diff1=( nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+n]-nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+n]);
      if( diff1*diff1 >ITER_TOL*ITER_TOL ) ambig[k]=0;
#if TEST2
      printf("diff1:%f ambig:%d\n",diff1,ambig[k]);
#endif
    }

    /* check if one ambiguous edge is collapsed (a point!) */
    /* this is the case if all control points on that edge are on the same spot */
    if(ambig[k])
    {
      collapsed[k]=1;
      for(j=1; j<nurbs[nr].v_npnt; j++)
      {
#if TEST2
        printf("x:%f ", nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+0]);
        printf("y:%f ", nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+1]);
        printf("z:%f ", nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+2]);
        printf("w:%f ", nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+3]);
        printf("\n");
#endif

        /* difference to the first point */
        for(n=0; n<4; n++)
        {
          diff1=( nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+n]-nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+n]);
          if( diff1*diff1 >ITER_TOL*ITER_TOL ) collapsed[k]=0;
        }
      }
    }
  }

  /* check u at vmin and vmax */
  for(k=0; k<2; k++)
  {
    if(k==0) j=0; else j=nurbs[nr].v_npnt-1;

#if TEST2
      printf("check at vmin and vmax\n");
    i=0;
      printf("x:%f ", nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+0]);
      printf("y:%f ", nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+1]);
      printf("z:%f ", nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+2]);
      printf("w:%f ", nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+3]);
      printf(" i:%d j:%d\n", i, j);
#endif
    /* check the first and last only */
    /* if they are equal its an ambiguous edge */ 
    /* difference to the first point */
    ambig[k]=1;
    i=nurbs[nr].u_npnt-1;
#if TEST2
      printf("x:%f ", nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+0]);
      printf("y:%f ", nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+1]);
      printf("z:%f ", nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+2]);
      printf("w:%f ", nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+3]);
      printf(" i:%d j:%d\n", i, j);
#endif
    for(n=0; n<4; n++)
    {
      diff1=(nurbs[nr].ctlarray[j*(nurbs[nr].v_stride)+n]-nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+n]);
      if( diff1*diff1 >ITER_TOL*ITER_TOL ) ambig[k]=0;
#if TEST2
      printf("diff1:%f ambig:%d\n",diff1,ambig[k]);
#endif
    }

    /* check if one ambiguous edge is collapsed (a point!) */
    /* this is the case if all control points on that edge are on the same spot */
    if(ambig[k])
    {
      collapsed[k]=1;
      for(i=1; i<nurbs[nr].u_npnt; i++)
      {
#if TEST2
        printf("x:%f ", nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+0]);
        printf("y:%f ", nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+1]);
        printf("z:%f ", nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+2]);
        printf("w:%f ", nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+3]);
        printf(" i:%d j:%d\n", i, j);
#endif

        /* difference to the first point */
        for(n=0; n<4; n++)
        {
          diff1=(nurbs[nr].ctlarray[j*(nurbs[nr].v_stride)+n]-nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+n]);
          if( diff1*diff1 >ITER_TOL*ITER_TOL ) collapsed[k]=0;
        }
      }
    }
  }


  if((ambig[2]==1)||(ambig[3]==1)) *axis=1;
  if((ambig[0]==1)||(ambig[1]==1)) *axis+=2;

  /* determine the case */
  /* -1error, 0plate, 1cyl, 2torus, 3ball, 4half-ball-bot, 5half-ball-top */
  /* ( half-ball-top: closed at top)  */
  i=0; for(k=0; k<4; k++)
  {
#if TEST2
    printf("collapsed=%d (1:yes)\n", collapsed[k]);
#endif
    i+=collapsed[k];
  }
  if(i==0)
  {
    if((*axis==1)||(*axis==2)) return(1);
    else if(*axis==3) return(2);
    else return(0);
  }
  else
  {
    if((collapsed[0]==1)&&(collapsed[1]==1)) { *axis=2; return(3); } 
    else if((collapsed[2]==1)&&(collapsed[3]==1)) { *axis=1; return(3); } 
    else if((collapsed[0]==1)&&(collapsed[1]==0)) { *axis=2; return(4); } 
    else if((collapsed[2]==1)&&(collapsed[3]==0)) { *axis=1; return(4); } 
    else if((collapsed[0]==0)&&(collapsed[1]==1)) { *axis=2; return(5); } 
    else if((collapsed[2]==0)&&(collapsed[3]==1)) { *axis=1; return(5); }
    else return(-1);
  }
}



void trimNurbs( Nurbs *nurbs, int nr, int patch, double TOL_AMBIG)
{
  int i,j,k,n,p;
  int curve, sum_p=0, sum_inverted=0, nurbsType=0,index=0, axis, buf=0, returnedPntsPerPnt  ;

  double convergTol,  normTol;
  int maxPass; 
  int numberOfControlPoints =nurbs[nr].u_npnt*nurbs[nr].v_npnt;
  snlCtrlPoint* controlPoints = new snlCtrlPoint [ numberOfControlPoints ];

  double v[3],maxv[3], minv[3];

  int *nppc;
  int disti=0, startFlag, pskip, nextUflag=0, nextVflag=0, lastUflag=-1, lastVflag=-1;
  double umin, umax,vmin,vmax, distp;
  double dist1, dist2, tol_ambig;

  PntProj *pntproj=NULL;
  snlSurfLocn* inverted=NULL;

  normTol=COS_TOL;
  convergTol=ITER_TOL;
  maxPass=MAX_PASS;

  if(nurbs[nr].ctlarray==NULL)
  {
    printf(" programmers error: nurbs ctlarray was not initialized\n");
    exit(0);
  }

  for(i=0; i<nurbs[nr].nc[patch]; i++) sum_p+=nurbs[nr].np[patch][i];

  // (1)      Create an array of control points from your cgx control points:
  for (i=0; i<nurbs[nr].u_npnt; i++)
  {
    for (j=0; j<nurbs[nr].v_npnt; j++)
    {
      controlPoints [ index++ ].components (
	  nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+0],
	  nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+1],
	  nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+2],
	  nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+3] );
    }
  }

  // (2)      Create an array of knots for each knot vector in your cgx surface:

  knot* knotsU = new knot [ nurbs[nr].u_nknt ];
  knot* knotsV = new knot [ nurbs[nr].v_nknt ];

  for(i=0; i<nurbs[nr].u_nknt; i++)  knotsU[i]=nurbs[nr].uknt[i];
  for(i=0; i<nurbs[nr].v_nknt; i++)  knotsV[i]=nurbs[nr].vknt[i];

  // (3)      Create an snlSurface:
  // You need to supply degreeU, degreeV, sizeU and sizeV.
  snlSurface* surface = new snlSurface ( nurbs[nr].u_exp, nurbs[nr].v_exp, nurbs[nr].u_npnt, nurbs[nr].v_npnt, controlPoints, knotsU, knotsV );

  // (4) Create an array of snlPoints to project onto the surface.
  snlPoint* toProject = new snlPoint [ sum_p ]; // You supply numPointsToProject.

  // Fill the snlPoint array with data from cgx
  // and determine the patch-x,y,z range  by looking into the points (tollerance of ambig)
  for(n=0; n<3; n++) { dist1=maxv[n]=-MAXVALUE; minv[n]=MAXVALUE; }
  p=0;
  for(curve=0; curve<nurbs[nr].nc[patch]; curve++)
  {
    j=0;
    for(i=0; i<nurbs[nr].np[patch][curve]; i++)
    {
      for(n=0;n<3; n++)
      {
	v[n]=nurbs[nr].xyz[patch][curve][j++];
        if(maxv[n]<v[n]) maxv[n]=v[n];
        if(minv[n]>v[n]) minv[n]=v[n];
      }
      toProject[p++].components (v[0],v[1],v[2]);
#if TEST
      printf("PNT P%d %f %f %f\n",p,v[0],v[1],v[2]);
#endif
    }
  }
  for(n=0; n<3; n++)
  {
    maxv[n]-=minv[n];
    if(dist1<maxv[n]) dist1=maxv[n];
  }
#if PRINT_SNL
  surface ->print_cpp();
  exit(-1);
#endif

  // determine the nurbs-x,y,z range  by looking into the points (tollerance of ambig)
  for(n=0; n<3; n++) { dist2=maxv[n]=-MAXVALUE; minv[n]=MAXVALUE; }
  for(i=0; i<nurbs[nr].u_npnt; i++)
  {
    for(j=1; j<nurbs[nr].v_npnt; j++)
    {
      for(n=0; n<3; n++)
      {
        if(maxv[n]<nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+n]) maxv[n]=nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+n];
        if(minv[n]>nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+n]) minv[n]=nurbs[nr].ctlarray[i*(nurbs[nr].v_stride*nurbs[nr].v_npnt)+j*(nurbs[nr].v_stride)+n];
      }
    }
  }
  for(n=0; n<3; n++)
  {
    maxv[n]-=minv[n];
    if(dist2<maxv[n]) dist2=maxv[n];
  }
  tol_ambig=dist1/dist2*TOL_AMBIG;

  // determine the nurbs-u,v range  by looking into the knots (to detect ambiguous points) 
  umin=nurbs[nr].uknt[0];
  umax=nurbs[nr].uknt[nurbs[nr].u_nknt-1];
  vmin=nurbs[nr].vknt[0];
  vmax=nurbs[nr].vknt[nurbs[nr].v_nknt-1];
#if TEST
  printf(" u_nknt:%d v_nknt:%d\n",nurbs[nr].u_nknt,nurbs[nr].v_nknt ); 
  printf(" umax:%f umin:%f vmax:%f vmin:%f\n", umax,umin,vmax,vmin); 
#endif

  // temporary buffer for the projected points
  pntproj= new PntProj [sum_p];

  // additional number of points per curve object
  nppc= new int [nurbs[nr].nc[patch]+1];

  // (5) Project points to surface.
  if(printFlag) printf("inversion of %d points from nurbs:%s\n", sum_p, nurbs[nr].name);
  returnedPntsPerPnt =MAX_RETURNED_PNTS_PER_PNT/5;

 finerProjection:;
  for(i=0; i<sum_p; i++) pntproj[i].n=0;
  inverted =  surface -> fastProject( toProject, sum_p, &sum_inverted, convergTol, normTol, maxPass,1, returnedPntsPerPnt );

  // sort all projections according to the point indexes in "toProject" 
  for(j=0; j<sum_inverted; j++)
  {
    pntproj[inverted[j].origPtIndex].paramU[pntproj[inverted[j].origPtIndex].n]
     =inverted[j].paramU;
    pntproj[inverted[j].origPtIndex].paramV[pntproj[inverted[j].origPtIndex].n]
     =inverted[j].paramV;
    pntproj[inverted[j].origPtIndex].dist[pntproj[inverted[j].origPtIndex].n]
     =inverted[j].dist;
    pntproj[inverted[j].origPtIndex].n++;
  }

  // determine the type of the nurbs 
  // (0plate, 1cyl, 2torus, 3ball, 4half-ball-bot, 5half-ball-top) 
  // ( half-ball-top: closed at top)  
  // axis=1 v is ambiguous (point has either vmax or vmin)
  // axis=2 u is ambiguous (point has either umax or umin)
  if((nurbsType=getNurbsType( nurbs, nr, &axis ))<0)
  { printf("ERROR: type of nurbs:%d not known\n",nurbsType); exit(-1); }
  if(printFlag) printf("nurbsType:%d (0plate, 1cyl, 2torus, 3ball, 4half-ball-bot, 5half-ball-top) axis:%d (1:u v:2 uv:3) tol_ambig:%f\n",nurbsType,axis,tol_ambig);

  // plate
  if(nurbsType==0)
  {
    n=0;
    for(curve=0; curve<nurbs[nr].nc[patch]; curve++)
    {
      j=0;
      // over all points of that curve
      // scip the last point to avoid differences between 1st and last
      for(k=0; k<nurbs[nr].np[patch][curve]-1; k++)
      {
        distp=MAXVALUE;
        for(i=0; i<pntproj[n].n; i++)
        {
          if(pntproj[n].dist[i]<distp)
          { distp=pntproj[n].dist[i]; disti=i; }
        }
        // check if the distance is "close"
        if((distp>MIN_DIST)&&(returnedPntsPerPnt<MAX_RETURNED_PNTS_PER_PNT))
        { returnedPntsPerPnt=MAX_RETURNED_PNTS_PER_PNT; delete[] inverted; goto finerProjection; }

        nurbs[nr].uv[patch][curve][j++]=pntproj[n].paramU[disti];
        nurbs[nr].uv[patch][curve][j++]=pntproj[n].paramV[disti];
#if TEST
        printf("pnt:%d dist:%f disti:%d uv:%f %f  orig-xyz: %f %f %f \n",n, distp,disti, pntproj[n].paramU[disti],pntproj[n].paramV[disti], toProject[n].x(), toProject[n].y(), toProject[n].z() );
#endif
        n++;
      }
  
      /* last point == 1st point */
      nurbs[nr].uv[patch][curve][j++]=nurbs[nr].uv[patch][curve][0];
      nurbs[nr].uv[patch][curve][j++]=nurbs[nr].uv[patch][curve][1];
      n++;
    }
  }
  else
  {
   // go to the first not ambiguous point and store the skiped ones 
   // go over all points and chose the closest one if ambiguous
   // finally care about the skipped leading points
   n=0;
   for(curve=0; curve<nurbs[nr].nc[patch]; curve++)
   {
    nppc[curve]=0;
    startFlag=1;
    pskip=0; 
    j=0;

    // over all points of that curve
    // scip the last point to avoid differences between 1st and last
    for(k=0; k<nurbs[nr].np[patch][curve]-1; k++)
    {
      // search the closest one to the nurbs
      distp=MAXVALUE;
      for(i=0; i<pntproj[n].n; i++) if(pntproj[n].dist[i]<distp)
      { distp=pntproj[n].dist[i]; disti=i; }

      // check if the distance is "close"
      if((distp>MIN_DIST)&&(returnedPntsPerPnt<MAX_RETURNED_PNTS_PER_PNT))
      { returnedPntsPerPnt=MAX_RETURNED_PNTS_PER_PNT; delete[] inverted; goto finerProjection; }
#if TEST
        printf("pnt:%d dist:%f disti:%d uv:%f %f  orig-xyz: %f %f %f \n",n, distp,disti, pntproj[n].paramU[disti],pntproj[n].paramV[disti], toProject[n].x(), toProject[n].y(), toProject[n].z() );
#endif

      // look for leading ambiguous points and skip the leading ones
      if(startFlag)
      {
        if((axis==2)||(axis==3))     // in udir
	{
          if(((pntproj[n].paramU[disti]- tol_ambig )<umin)||((pntproj[n].paramU[disti]+ tol_ambig )>umax)) { pskip++; goto skip; } 
	}
        if((axis==1)||(axis==3))     // in vdir 
	{
          if(((pntproj[n].paramV[disti]- tol_ambig )<vmin)||((pntproj[n].paramV[disti]+ tol_ambig )>vmax)) { pskip++; goto skip; }
	}
      }
      startFlag=0;

      // if the point is ambig look which coord would give the shortest dist to the prev point
      if((axis==2)||(axis==3))     // in udir
      {
       if(((pntproj[n].paramU[disti]-tol_ambig)<umin)||((pntproj[n].paramU[disti]+tol_ambig)>umax))
       {
	dist1=nurbs[nr].uv[patch][curve][j-2]-umin;
	dist2=nurbs[nr].uv[patch][curve][j-2]-umax;
	if((dist1*dist1)<(dist2*dist2)) pntproj[n].paramU[disti]=umin;
        else  pntproj[n].paramU[disti]=umax;
#if TEST
        printf("dist1:%f dist2:%f corr.u:%f\n", dist1, dist2, pntproj[n].paramU[disti]);
#endif
       }
      }
      if((axis==1)||(axis==3))     // in vdir 
      {
       if(((pntproj[n].paramV[disti]-tol_ambig)<vmin)||((pntproj[n].paramV[disti]+tol_ambig)>vmax))
       {
	dist1=nurbs[nr].uv[patch][curve][j-1]-vmin;
	dist2=nurbs[nr].uv[patch][curve][j-1]-vmax;
	if((dist1*dist1)<(dist2*dist2)) pntproj[n].paramV[disti]=vmin;
        else  pntproj[n].paramV[disti]=vmax;
#if TEST
        printf("dist1:%f dist2:%f corr.v:%f\n", dist1, dist2, pntproj[n].paramV[disti]);
#endif
       }
      }
#if TEST
      printf("corr pnt:%d uv:%f %f  orig-xyz: %f %f %f \n",n, pntproj[n].paramU[disti],pntproj[n].paramV[disti], toProject[n].x(), toProject[n].y(), toProject[n].z() );
#endif

      // if the coordinates of this point have to be used for the previously duplicated one
      if(nextUflag>0)
      { nurbs[nr].uv[patch][curve][nextUflag]=pntproj[n].paramU[disti]; nextUflag=0; } 
      if(nextVflag>0)
      { nurbs[nr].uv[patch][curve][nextVflag]=pntproj[n].paramV[disti]; nextVflag=0; }

      // if collapsed edges exist add an additional point
      if(nurbsType>2)
      {
        if(axis==1)     // in udir
        {
	  if(nurbsType==3)
          { if(((pntproj[n].paramU[disti]-tol_ambig)<umin)||((pntproj[n].paramU[disti]+tol_ambig)>umax)) buf=1; else buf=0; }
	  if(nurbsType==4)
          { if((pntproj[n].paramU[disti]-tol_ambig)<umin) buf=1; else buf=0; }
          if(nurbsType==5)
          { if((pntproj[n].paramU[disti]+tol_ambig)>umax) buf=1; else buf=0; }
          if(buf)  // on the collapsed edge
          {
            // additional point required       
            // use the v value from the previous point
#if TEST
            if(j>2) printf("use v-val from prev point:%f\n", nurbs[nr].uv[patch][curve][j-2] );
            else    printf("use v-val from prev point: (wait till end)\n" );
#endif
            nurbs[nr].uv[patch][curve][j++]=pntproj[n].paramU[disti];
            if(j>2) { nurbs[nr].uv[patch][curve][j]=nurbs[nr].uv[patch][curve][j-2];j++; }
            else lastVflag=j++;

            // include a new point
            nppc[curve]++;
            if( (nurbs[nr].uv[patch][curve] = (GLfloat *)realloc( (GLfloat *)nurbs[nr].uv[patch][curve], ((nurbs[nr].np[patch][curve]+1+nppc[curve])*2)*sizeof(GLfloat) ))==NULL )
            { printf(" ERROR: realloc failure in trimNurbs(), nurbs:%s can not be shaped\n\n", nurbs[nr].name);  return; }

            // use the v value from the next point
            // the value is not really known now, save a flag for later correction
            nurbs[nr].uv[patch][curve][j++]=pntproj[n].paramU[disti];
            distp=MAXVALUE; for(i=0; i<pntproj[n+1].n; i++) if(pntproj[n+1].dist[i]<distp)
            { distp=pntproj[n+1].dist[i]; nurbs[nr].uv[patch][curve][j]=pntproj[n+1].paramV[i]; }
            // check if the distance is "close"
            if((distp>MIN_DIST)&&(returnedPntsPerPnt<MAX_RETURNED_PNTS_PER_PNT))
            { returnedPntsPerPnt=MAX_RETURNED_PNTS_PER_PNT; delete[] inverted; goto finerProjection; }
            nextVflag=j++;
#if TEST
            printf("collapsed at u:%f add pnt:%d\n", pntproj[n].paramU[disti], j/2 );
#endif
          }
          else // not on the collapsed edge
          {
            nurbs[nr].uv[patch][curve][j++]=pntproj[n].paramU[disti];
            nurbs[nr].uv[patch][curve][j++]=pntproj[n].paramV[disti];
#if TEST
            printf(" ori:%d uv:%f %f\n", n, nurbs[nr].uv[patch][curve][j-2], nurbs[nr].uv[patch][curve][j-1] );
#endif
          }
        }
        else     // in vdir
        {
	  if(nurbsType==3)
          { if(((pntproj[n].paramV[disti]-tol_ambig)<vmin)||((pntproj[n].paramV[disti]+tol_ambig)>vmax)) buf=1; else buf=0; }
	  if(nurbsType==4)
          { if((pntproj[n].paramV[disti]-tol_ambig)<vmin) buf=1; else buf=0; }
          if(nurbsType==5)
          { if((pntproj[n].paramV[disti]+tol_ambig)>vmax) buf=1; else buf=0; }
          if(buf)  // on the collapsed edge
          {
            // additional point required       
            // use the u value from the previous point
#if TEST
            if(j>2) printf("use u-val from prev point:%f\n", nurbs[nr].uv[patch][curve][j-2] );
            else    printf("use u-val from prev point: (wait till end)\n" );
#endif
            if(j>2) { nurbs[nr].uv[patch][curve][j]=nurbs[nr].uv[patch][curve][j-2];j++; }
            else lastUflag=j++;
            nurbs[nr].uv[patch][curve][j++]=pntproj[n].paramV[disti];

            // include a new point
            nppc[curve]++;
            if( (nurbs[nr].uv[patch][curve] = (GLfloat *)realloc( (GLfloat *)nurbs[nr].uv[patch][curve], ((nurbs[nr].np[patch][curve]+1+nppc[curve])*2)*sizeof(GLfloat) ))==NULL )
            { printf(" ERROR: realloc failure in trimNurbs(), nurbs:%s can not be shaped\n\n", nurbs[nr].name);  return; }

            // use the u value from the next point
            // the value is not really known now, save a flag for later correction
            distp=MAXVALUE; for(i=0; i<pntproj[n+1].n; i++) if(pntproj[n+1].dist[i]<distp)
            { distp=pntproj[n+1].dist[i]; nurbs[nr].uv[patch][curve][j]=pntproj[n+1].paramU[i]; }

            // check if the distance is "close"
            if((distp>MIN_DIST)&&(returnedPntsPerPnt<MAX_RETURNED_PNTS_PER_PNT))
            { returnedPntsPerPnt=MAX_RETURNED_PNTS_PER_PNT; delete[] inverted; goto finerProjection; }
            nextUflag=j++;
            nurbs[nr].uv[patch][curve][j++]=pntproj[n].paramV[disti];
#if TEST
            printf("collapsed at v:%f add pnt:%d\n", pntproj[n].paramV[disti], j/2 );
#endif
          }
          else // not on the collapsed edge
          {
            nurbs[nr].uv[patch][curve][j++]=pntproj[n].paramU[disti];
            nurbs[nr].uv[patch][curve][j++]=pntproj[n].paramV[disti]; 
#if TEST
            printf(" ori:%d uv:%f %f\n", n, nurbs[nr].uv[patch][curve][j-2], nurbs[nr].uv[patch][curve][j-1] );
#endif
          }
        }
      }
      // no collapsed edges
      else
      {
        nurbs[nr].uv[patch][curve][j++]=pntproj[n].paramU[disti];
        nurbs[nr].uv[patch][curve][j++]=pntproj[n].paramV[disti];
      }
      
    skip:;
      n++;
    }

#if TEST
    printf("\nLEADING ambiguous points:%d\n\n",pskip);
#endif

    // store leading ambiguous points 
    for(k=n-nurbs[nr].np[patch][curve]+1; k<n-nurbs[nr].np[patch][curve]+1+pskip; k++)
    {
      // search the closest one to the nurbs
      distp=MAXVALUE;
      for(i=0; i<pntproj[k].n; i++) if(pntproj[k].dist[i]<distp)
      { distp=pntproj[k].dist[i]; disti=i; }

      // check if the distance is "close"
      if((distp>MIN_DIST)&&(returnedPntsPerPnt<MAX_RETURNED_PNTS_PER_PNT))
      { returnedPntsPerPnt=MAX_RETURNED_PNTS_PER_PNT; delete[] inverted; goto finerProjection; }
#if TEST
        printf("pnt:%d dist:%f disti:%d uv:%f %f  orig-xyz: %f %f %f \n", k, distp,disti, pntproj[k].paramU[disti],pntproj[k].paramV[disti], toProject[k].x(), toProject[k].y(), toProject[k].z() );
#endif

      //  look which coord would give the shortest dist to the last-1 point of the curve
      if((axis==2)||(axis==3))     // in udir
      {
       if(((pntproj[k].paramU[disti]-tol_ambig)<umin)||((pntproj[k].paramU[disti]+tol_ambig)>umax))
       {
	dist1=nurbs[nr].uv[patch][curve][j-2]-umin;
	dist2=nurbs[nr].uv[patch][curve][j-2]-umax;
	if((dist1*dist1)<(dist2*dist2)) pntproj[k].paramU[disti]=umin;
        else  pntproj[k].paramU[disti]=umax;
#if TEST
        printf("dist1:%f dist2:%f corr.u:%f\n", dist1, dist2, pntproj[k].paramU[disti]);
#endif
       }
      }
      if((axis==1)||(axis==3))     // in vdir 
      {
       if(((pntproj[k].paramV[disti]-tol_ambig)<vmin)||((pntproj[k].paramV[disti]+tol_ambig)>vmax))
       {
	dist1=nurbs[nr].uv[patch][curve][j-1]-vmin;
	dist2=nurbs[nr].uv[patch][curve][j-1]-vmax;
	if((dist1*dist1)<(dist2*dist2)) pntproj[k].paramV[disti]=vmin;
        else  pntproj[k].paramV[disti]=vmax;
#if TEST
        printf("dist1:%f dist2:%f corr.v:%f\n", dist1, dist2, pntproj[k].paramV[disti]);
#endif
       }
      }

      // if the coordinates of this point have to be used for the previously duplicated one
      if(nextUflag>0)
      { nurbs[nr].uv[patch][curve][nextUflag]=pntproj[k].paramU[disti]; nextUflag=0; } 
      if(nextVflag>0)
      { nurbs[nr].uv[patch][curve][nextVflag]=pntproj[k].paramV[disti]; nextVflag=0; }

      // if collapsed edges exist add an additional point
      if(nurbsType>2)
      {
        if(axis==1)     // in udir
        {
	  if(nurbsType==3)
          { if(((pntproj[k].paramU[disti]-tol_ambig)<umin)||((pntproj[k].paramU[disti]+tol_ambig)>umax)) buf=1; else buf=0; }
	  if(nurbsType==4)
          { if((pntproj[k].paramU[disti]-tol_ambig)<umin) buf=1; else buf=0; }
          if(nurbsType==5)
          { if((pntproj[k].paramU[disti]+tol_ambig)>umax) buf=1; else buf=0; }
          if(buf)  // on the collapsed edge
          {
            // additional point required       
            // use the v value from the previous point
#if TEST
            if(j>2) printf("use v-val from prev point:%f\n", nurbs[nr].uv[patch][curve][j-2] );
            else    printf("use v-val from prev point: (wait till end)\n" );
#endif
            nurbs[nr].uv[patch][curve][j]=nurbs[nr].uv[patch][curve][j-2]; j++;
            if(j>2) { nurbs[nr].uv[patch][curve][j]=nurbs[nr].uv[patch][curve][j-2];j++; }
            else lastVflag=j++;

            // include a new point
            nppc[curve]++;
            if( (nurbs[nr].uv[patch][curve] = (GLfloat *)realloc( (GLfloat *)nurbs[nr].uv[patch][curve], ((nurbs[nr].np[patch][curve]+1+nppc[curve])*2)*sizeof(GLfloat) ))==NULL )
            { printf(" ERROR: realloc failure in trimNurbs(), nurbs:%s can not be shaped\n\n", nurbs[nr].name);  return; }

            // use the v value from the next point
            // the value is not really known now, save a flag for later correction
            nurbs[nr].uv[patch][curve][j++]=pntproj[k].paramU[disti];
            distp=MAXVALUE; for(i=0; i<pntproj[k+1].n; i++) if(pntproj[k+1].dist[i]<distp)
            { distp=pntproj[k+1].dist[i]; nurbs[nr].uv[patch][curve][j]=pntproj[k+1].paramV[i]; }
            // check if the distance is "close"
            if((distp>MIN_DIST)&&(returnedPntsPerPnt<MAX_RETURNED_PNTS_PER_PNT))
            { returnedPntsPerPnt=MAX_RETURNED_PNTS_PER_PNT; delete[] inverted; goto finerProjection; }
            nextVflag=j++;
#if TEST
            printf("collapsed at u:%f add pnt:%d\n", pntproj[k].paramU[disti], j/2 );
#endif
          }
          else // not on the collapsed edge
          {
            nurbs[nr].uv[patch][curve][j++]=pntproj[k].paramU[disti];
            nurbs[nr].uv[patch][curve][j++]=pntproj[k].paramV[disti];
#if TEST
            printf(" ori:%d uv:%f %f\n", k, nurbs[nr].uv[patch][curve][j-2], nurbs[nr].uv[patch][curve][j-1] );
#endif
          }
        }
        else     // in vdir
        {
	  if(nurbsType==3)
          { if(((pntproj[k].paramU[disti]-tol_ambig)<umin)||((pntproj[k].paramU[disti]+tol_ambig)>umax)) buf=1; else buf=0; }
	  if(nurbsType==4)
          { if((pntproj[k].paramV[disti]-tol_ambig)<vmin) buf=1; else buf=0; }
          if(nurbsType==5)
          { if((pntproj[k].paramV[disti]+tol_ambig)>vmax) buf=1; else buf=0; }
          if(buf)  // on the collapsed edge
          {
            // additional point required       
            // use the u value from the previous point
#if TEST
            if(j>2) printf("use u-val from prev point:%f\n", nurbs[nr].uv[patch][curve][j-2] );
            else    printf("use u-val from prev point: (wait till end)\n" );
#endif
            if(j>2) { nurbs[nr].uv[patch][curve][j]=nurbs[nr].uv[patch][curve][j-2];j++; }
            else lastUflag=j++;
            nurbs[nr].uv[patch][curve][j++]=pntproj[k].paramV[disti];

            // include a new point
            nppc[curve]++;
            if( (nurbs[nr].uv[patch][curve] = (GLfloat *)realloc( (GLfloat *)nurbs[nr].uv[patch][curve], ((nurbs[nr].np[patch][curve]+1+nppc[curve])*2)*sizeof(GLfloat) ))==NULL )
            { printf(" ERROR: realloc failure in trimNurbs(), nurbs:%s can not be shaped\n\n", nurbs[nr].name);  return; }

            // use the u value from the next point
            // the value is not really known now, save a flag for later correction
            distp=MAXVALUE; for(i=0; i<pntproj[k+1].n; i++) if(pntproj[k+1].dist[i]<distp)
            { distp=pntproj[k+1].dist[i]; nurbs[nr].uv[patch][curve][j]=pntproj[k+1].paramU[i]; }

            // check if the distance is "close"
            if((distp>MIN_DIST)&&(returnedPntsPerPnt<MAX_RETURNED_PNTS_PER_PNT))
            { returnedPntsPerPnt=MAX_RETURNED_PNTS_PER_PNT; delete[] inverted; goto finerProjection; }
            nextUflag=j++;
            nurbs[nr].uv[patch][curve][j++]=pntproj[k].paramV[disti];
#if TEST
            printf("collapsed at v:%f add pnt:%d\n", pntproj[k].paramV[disti], j/2 );
#endif
          }
          else // not on the collapsed edge
          {
            nurbs[nr].uv[patch][curve][j++]=pntproj[k].paramU[disti];
            nurbs[nr].uv[patch][curve][j++]=pntproj[k].paramV[disti]; 
#if TEST
            printf(" ori:%d uv:%f %f\n", k, nurbs[nr].uv[patch][curve][j-2], nurbs[nr].uv[patch][curve][j-1] );
#endif
          }
        }
      }
      else
      {
        nurbs[nr].uv[patch][curve][j++]=pntproj[k].paramU[disti];
        nurbs[nr].uv[patch][curve][j++]=pntproj[k].paramV[disti];
      }
    }

    // if the coordinates of this point have to be used
    if(lastUflag>-1)
    { nurbs[nr].uv[patch][curve][lastUflag]=nurbs[nr].uv[patch][curve][j-2]; lastUflag=-1; } 
    if(lastVflag>-1)
    { nurbs[nr].uv[patch][curve][lastVflag]=nurbs[nr].uv[patch][curve][j-1]; lastVflag=-1; }

    // last point == 1st point
    nurbs[nr].uv[patch][curve][j++]=nurbs[nr].uv[patch][curve][0];
    nurbs[nr].uv[patch][curve][j++]=nurbs[nr].uv[patch][curve][1];

    // if the coordinates of this point have to be used for the previously duplicated one
    if(nextUflag>0)
    { nurbs[nr].uv[patch][curve][nextUflag]=nurbs[nr].uv[patch][curve][j-2]; } 
    if(nextVflag>0)
    { nurbs[nr].uv[patch][curve][nextVflag]=nurbs[nr].uv[patch][curve][j-1]; }

#if TEST
    // check if we have too many points now
    if(((nurbs[nr].np[patch][curve]+nppc[curve])*2)!=j)
    {
      printf(" ERROR too many points:%d/2 (!=%d) in curve)\n", j,nurbs[nr].np[patch][curve]);
      exit(-1);
    }
#endif
    n++;
   }
   // adapt the nurbs-parameter
   for(curve=0; curve<nurbs[nr].nc[patch]; curve++) nurbs[nr].np[patch][curve]+=nppc[curve];
  }

  delete[] pntproj;
  delete nppc;         // Release number of points per curve object
  delete surface;      // Release surface object.
  delete[] inverted;   // Release snlVertex array returned from projection function.
  delete[] toProject;  // Release points that were inverted onto surface.

  // The rest of the objects are deleted by libSNL.
}
