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

#include <cgx.h>

#define MAX_GSUR_PARAMERTER 7
extern SpecialSet specialset[1];
extern Eqal eqal;
extern char entity_k[SET_COLS];
extern Psets     *pset;

int writefbd(char *setname, Summen *anz, SumGeo *anzGeo, Points *pnt, Lines *line, Lcmb *lcmb, Gsur *surf, Gbod *body, Nurbs *nurbs, Sets *set, Shapes *shape )
{
  FILE *handle1;
  int  i, j, k, n, l, se, nr;
  char string[MAX_LINE_LENGTH], datout[MAX_LINE_LENGTH], buffer[MAX_LINE_LENGTH];
  int length, params, setNr;
  int bias_fbd;
  char attr;

  setNr=getSetNr(setname);
  if (setNr<0)
  {
    printf (" ERROR: set:%s does not exist\n", setname);
    return (-1);
  }

  length=sword( setname, string);
  string[length]='.';
  string[length+1]='f';
  string[length+2]='b';
  string[length+3]='d';
  string[length+4]='\0';
  strcpy( datout, string);
  /* Open the files and check to see that it was opened correctly */
  handle1 = fopen (datout, "w+b");
  if (handle1==NULL) { printf ("\nThe output file \"%s\" could not be opened.\n\n", datout); return(-1);}
  else  printf ("\n%s opened\n",datout);

  printf ("\n write fbd  \n");
  delSet(specialset->tmp);

  if (anzGeo->p)
  {
    if( compare( setname, "all", 3)==3)
      for (i=0; i<anzGeo->p; i++)
      {
        if( pnt[i].name != (char *)NULL )
          fprintf (handle1, " PNT %-5s %10.5f %14.5f %14.5f \n", pnt[i].name, pnt[i].px,pnt[i].py,pnt[i].pz);
      }
    else
    {
      for (i=0; i<set[setNr].anz_p; i++)
      {
        j=set[setNr].pnt[i];
        fprintf (handle1, " PNT %-5s %10.5f %14.5f %14.5f \n", pnt[j].name, pnt[j].px,pnt[j].py,pnt[j].pz);
      }
    }
  }

  if (anz->sets)
  {
    if( compare( setname, "all", 3)==3)
    {
      for (i=1; i<anz->sets; i++)
      {
        if(set[i].type) /* SEQ */
        {
          if ( ((compare( set[i].name, setname, length) == length)
          && (length == strlen(setname) ))
          || (( compare( setname, "all", 3) == 3) 
	      &&(set[i].name != (char *)NULL )) )
          {
            params=MAX_PARAM_PER_RECORD/2;
            if(set[i].anz_p)
            {
              strcpy(string,"   ");
              fprintf (handle1, " SEQA %s %s pnt ", set[i].name, string );
              for (j=0; j<set[i].anz_p; j++)
              {
                if(j==params)
        	   {
        	     params+=MAX_PARAM_PER_RECORD/2;
                  strcpy(string,"END");
                  fprintf (handle1, "\n SEQA %s %s pnt ", set[i].name, string );
                }
                fprintf (handle1, " %s", pnt[set[i].pnt[j]].name);
              }
              fprintf (handle1, "\n");
            }
          }
        }
      }
    }
    else
    {
      if ( set[setNr].name != (char *)NULL )
      {
        if(set[setNr].type==1)
	{
          if(set[setNr].anz_n)
          {
            params=MAX_PARAM_PER_RECORD/2;
            strcpy(string,"   ");
            fprintf (handle1, " SEQA %s %s nod ", set[setNr].name, string );
            for (j=0; j<set[setNr].anz_n; j++)
            {
              if(j==params)
              {
          	params+=MAX_PARAM_PER_RECORD/2;
                strcpy(string,"END");
                fprintf (handle1, "\n SEQA %s %s nod ", set[setNr].name, string );
              }
              fprintf(handle1, "%d ",  set[setNr].node[j]);
            }
            fprintf (handle1, "\n");
          }
	}

        for (i=0; i<set[setNr].anz_l; i++)
        {
          l=set[setNr].line[i];

          if (line[l].typ=='s')     
          {
            se=line[l].trk;
            if ( set[se].name != (char *)NULL )
            {
              params=MAX_PARAM_PER_RECORD/2;
              if(set[se].anz_p)
              {
                strcpy(string,"   ");
                fprintf (handle1, " SEQA %s %s pnt ", set[se].name, string );
                for (j=0; j<set[se].anz_p; j++)
                {
                  if(j==params)
          	  {
          	     params+=MAX_PARAM_PER_RECORD/2;
                    strcpy(string,"END");
                    fprintf (handle1, "\n SEQA %s %s pnt ", set[se].name, string );
                  }
                  fprintf (handle1, " %s", pnt[set[se].pnt[j]].name);
                }
                fprintf (handle1, "\n");
              }
            }
          }
        }
      }
    }
  }

  if (anzGeo->l)
  {
    if( compare( setname, "all", 3)==3)
    {
      for (i=0; i<anzGeo->l; i++)
      {
        /* load bias and division in buffer */
        bias_fbd=getBias_fbd(i,line);
        if(line[i].div<10) sprintf(buffer, "%d0%d", bias_fbd, line[i].div);
        else if(line[i].div<100) sprintf(buffer, "%d%d", bias_fbd, line[i].div);

        if (line[i].typ=='a') /* its an arc-line */
        {
          if( line[i].name != (char *)NULL )
            fprintf (handle1, " LINE %s %s %s %s %s\n",
            line[i].name, pnt[line[i].p1].name,
            pnt[line[i].p2].name, pnt[line[i].trk].name, buffer );
        }
        else if (line[i].typ=='s') /* its a seq-line */
        {
          /* check if the seq is correct defined */
          /* - the endpoints must be included in the seq */
          if( ((line[i].p1!=set[line[i].trk].pnt[0])&&(line[i].p2!=set[line[i].trk].pnt[0]))
            || ((line[i].p1!=set[line[i].trk].pnt[set[line[i].trk].anz_p-1])&&(line[i].p2!=set[line[i].trk].pnt[set[line[i].trk].anz_p-1])) )
	  {
            /* closer check */
            k=0;
            for(j=0; j<set[line[i].trk].anz_p; j++)
	    {
              if(set[line[i].trk].pnt[j]==line[i].p1) k++;
              if(set[line[i].trk].pnt[j]==line[i].p2) k++;
	    }
            if(k!=2)
	    {
              /* add line to the error-set */
              pre_seta( specialset->tmp, "l", line[i].name);
	    }
	  }

          if( line[i].name != (char *)NULL )
            fprintf (handle1, " LINE %s %s %s %s %s\n",
            line[i].name, pnt[line[i].p1].name,
            pnt[line[i].p2].name, set[line[i].trk].name, buffer );
        }
        else
        {
          if( line[i].name != (char *)NULL )
            fprintf (handle1, " LINE %s %s %s %s\n",
            line[i].name, pnt[line[i].p1].name, pnt[line[i].p2].name, buffer );
        }
      }
    }
    else
    {
      for (i=0; i<set[setNr].anz_l; i++)
      {
        j=set[setNr].line[i];

        /* load bias and division in buffer */
        bias_fbd=getBias_fbd(j,line);
        if(line[j].div<10) sprintf(buffer, "%d0%d", bias_fbd, line[j].div);
        else if(line[j].div<100) sprintf(buffer, "%d%d", bias_fbd, line[j].div);

        if (line[j].typ=='a') /* its an arc-line */
        {
          fprintf (handle1, " LINE %s %s %s %s %s\n",
            line[j].name, pnt[line[j].p1].name,
            pnt[line[j].p2].name, pnt[line[j].trk].name, buffer );
        }
        else if (line[j].typ=='s') /* its an seq-line */
        {
          /* check if the seq is correct defined */
          /* - the endpoints must be included in the seq */
          if( ((line[j].p1!=set[line[j].trk].pnt[0])&&(line[j].p2!=set[line[j].trk].pnt[0]))
            || ((line[j].p1!=set[line[j].trk].pnt[set[line[j].trk].anz_p-1])&&(line[j].p2!=set[line[j].trk].pnt[set[line[j].trk].anz_p-1])) )
	  {
            /* closer check */
            k=0;
            for(n=0; n<set[line[j].trk].anz_p; n++)
	    {
              if(set[line[j].trk].pnt[n]==line[j].p1) k++;
              if(set[line[j].trk].pnt[n]==line[j].p2) k++;
	    }
            if(k!=2)
	    {
              /* add line to the error-set */
              pre_seta( specialset->tmp, "l", line[j].name);
	    }
	  }

          fprintf (handle1, " LINE %s %s %s %s %s\n",
            line[j].name, pnt[line[j].p1].name, 
            pnt[line[j].p2].name, set[line[j].trk].name, buffer );
        }
        else
        {
          fprintf (handle1, " LINE %s %s %s %s\n",
            line[j].name, pnt[line[j].p1].name, pnt[line[j].p2].name, buffer );
        }
      }
    }
  }

  if (anzGeo->c)
  {
    if( compare( setname, "all", 3)==3)
    {
      for (i=0; i<anzGeo->c; i++)
      {
        if( lcmb[i].name != (char *)NULL )
        {
          fprintf (handle1, " LCMB %s ", lcmb[i].name );
          k=0;
          for (j=0; j<lcmb[i].nl; j++)
	  {
            if (k<8) 
	    {
              fprintf (handle1, " %1c %s", lcmb[i].o[j], line[lcmb[i].l[j]].name );
              k++;
            }
            else
	    {
              k=0; fprintf (handle1, "\n LCMB %s ADD ", lcmb[i].name );
              fprintf (handle1, " %1c %s", lcmb[i].o[j], line[lcmb[i].l[j]].name );
              k++;
            }
          }
          fprintf (handle1, " \n");
        }
      }
    }
    else
    {
      for (i=0; i<set[setNr].anz_c; i++)
      {
          n=set[setNr].lcmb[i];
          fprintf (handle1, " LCMB %s ", lcmb[n].name );
          k=0;
          for (j=0; j<lcmb[n].nl; j++)
	  {
            if (k<8) 
	    {
              fprintf (handle1, " %1c %s", lcmb[n].o[j], line[lcmb[n].l[j]].name );
              k++;
            }
            else
	    {
              k=0; fprintf (handle1, "\n LCMB %s ADD ", lcmb[n].name );
              fprintf (handle1, " %1c %s", lcmb[n].o[j], line[lcmb[n].l[j]].name );
              k++;
            }
          }
          fprintf (handle1, " \n");
      }
    }
  }

  if (anzGeo->sh)
  {
    if( compare( setname, "all", 3)==3)
    {
      for (i=0; i<anzGeo->sh; i++)
      {
        if( shape[i].name != (char *)NULL )
        {
          if(shape[i].type==0) fprintf (handle1, " SHPE %s PLN %s %s %s \n", shape[i].name, pnt[shape[i].p[0]].name, pnt[shape[i].p[1]].name, pnt[shape[i].p[2]].name );
	}
      }
    }
    else
    {
      for (i=0; i<set[setNr].anz_sh; i++)
      {
        n=set[setNr].shp[i];
        if(shape[n].type==0) fprintf (handle1, " SHPE %s PLN %s %s %s \n", shape[n].name, pnt[shape[n].p[0]].name, pnt[shape[n].p[1]].name, pnt[shape[n].p[2]].name );
      }
    }
  }

  if (anzGeo->nurs)
  {
    if( compare( setname, "all", 3)==3)
      for (nr=0; nr<anzGeo->nurs; nr++)
      {
        if( nurbs[nr].name != (char *)NULL )
        {
          fprintf(handle1," NURS %s DEFINE %4d %4d %4d %4d %4d %4d\n", nurbs[nr].name, nurbs[nr].u_exp, nurbs[nr].v_exp, nurbs[nr].u_npnt, nurbs[nr].v_npnt, nurbs[nr].u_nknt, nurbs[nr].v_nknt);

          for (j=0; j<nurbs[nr].u_npnt; j++)
            for (k=0; k<nurbs[nr].v_npnt; k++)
              fprintf (handle1, " NURS %s CONTROL %4d %4d %s %f\n", nurbs[nr].name, j+1, k+1, pnt[nurbs[nr].ctlpnt[j][k]].name, nurbs[nr].weight[j][k] );

          for (j=0; j<nurbs[nr].u_nknt; j++)
            fprintf(handle1, " NURS %s KNOT  U %4d %f\n", nurbs[nr].name, j+1, nurbs[nr].uknt[j]);
          for (j=0; j<nurbs[nr].v_nknt; j++)
            fprintf(handle1, " NURS %s KNOT  V %4d %f\n", nurbs[nr].name, j+1, nurbs[nr].vknt[j]);

          fprintf(handle1, " NURS %s END\n", nurbs[nr].name);
        }
      }
    else
    {
      for (i=0; i<set[setNr].anz_nurs; i++)
      {
        nr=set[setNr].nurs[i];

          fprintf(handle1," NURS %s DEFINE %4d %4d %4d %4d %4d %4d\n", nurbs[nr].name, nurbs[nr].u_exp, nurbs[nr].v_exp, nurbs[nr].u_npnt, nurbs[nr].v_npnt, nurbs[nr].u_nknt, nurbs[nr].v_nknt);

          for (j=0; j<nurbs[nr].u_npnt; j++)
            for (k=0; k<nurbs[nr].v_npnt; k++)
              fprintf (handle1, " NURS %s CONTROL %4d %4d %s %f\n", nurbs[nr].name, j+1, k+1, pnt[nurbs[nr].ctlpnt[j][k]].name, nurbs[nr].weight[j][k]);

          for (j=0; j<nurbs[nr].u_nknt; j++)
            fprintf(handle1, " NURS %s KNOT  U %4d %f\n", nurbs[nr].name, j+1, nurbs[nr].uknt[j]);
          for (j=0; j<nurbs[nr].v_nknt; j++)
            fprintf(handle1, " NURS %s KNOT  V %4d %f\n", nurbs[nr].name, j+1, nurbs[nr].vknt[j]);

          fprintf(handle1, " NURS %s END\n", nurbs[nr].name);

      }
    }
  }

  /* must follow the nurs, shapes */
  if (anzGeo->s)
  {
    if( compare( setname, "all", 3)==3)
    {
      for (i=0; i<anzGeo->s; i++)
      {
        if( surf[i].name != (char *)NULL )
        {
          if(surf[i].sh==-1) fprintf (handle1, " GSUR %s %1c BLEND", surf[i].name, surf[i].ori );
          else if(surf[i].sh>-1) fprintf (handle1, " GSUR %s %1c %s", surf[i].name, surf[i].ori, shape[surf[i].sh].name );
          params=0; 
          for (j=0; j<surf[i].nl; j++)
          {
            if (params==MAX_GSUR_PARAMERTER)
            { params=0; fprintf (handle1, "\n GSUR %s ADD", surf[i].name); }
            else  params++; 
            if (surf[i].typ[j]=='l')
            {
              fprintf (handle1, " %1c %s", surf[i].o[j], line[surf[i].l[j]].name );
            }
            else
            {
              fprintf (handle1, " %1c %s", surf[i].o[j], lcmb[surf[i].l[j]].name );
            }
          }
          fprintf (handle1, " \n");
        }
      }
    }
    else
    {
      for (i=0; i<set[setNr].anz_s; i++)
      {
          n=set[setNr].surf[i];
          if(surf[n].sh==-1) fprintf (handle1, " GSUR %s %1c BLEND", surf[n].name, surf[n].ori );
          else if(surf[n].sh>-1) fprintf (handle1, " GSUR %s %1c %s", surf[n].name, surf[n].ori, shape[surf[n].sh].name );
          params=0; 
          for (j=0; j<surf[n].nl; j++)
          {
            if (params==MAX_GSUR_PARAMERTER)
            { params=0; fprintf (handle1, "\n GSUR %s ADD", surf[n].name); }
            else  params++; 
            if (surf[n].typ[j]=='l')
            {
              fprintf (handle1, " %1c %s", surf[n].o[j], line[surf[n].l[j]].name );
            }
            else
            {
              fprintf (handle1, " %1c %s", surf[n].o[j], lcmb[surf[n].l[j]].name );
            }
        }
        fprintf (handle1, " \n");
      }
    }
  }

  if (anzGeo->b)
  {
    if( compare( setname, "all", 3)==3)
      for (i=0; i<anzGeo->b; i++)
      {
        if( body[i].name != (char *)NULL )
        {
          fprintf (handle1, " GBOD %s NORM", body[i].name );
          for (j=0; j<body[i].ns; j++)
            fprintf (handle1, " %1c %s", body[i].o[j], surf[body[i].s[j]].name );
          fprintf (handle1, " \n");
        }
      }
    else
    {
      for (i=0; i<set[setNr].anz_b; i++)
      {
        n=set[setNr].body[i];
        fprintf (handle1, " GBOD %s NORM", body[n].name );
        for (j=0; j<body[n].ns; j++)
          fprintf (handle1, " %1c %s", body[n].o[j], surf[body[n].s[j]].name );
        fprintf (handle1, " \n");
      }
    }
  }

  /* write the sets */
  for (i=1; i<anz->sets; i++)
  {
    if((!set[i].type)&&( set[i].name!=(char *)NULL)&&( set[i].name[0]!='_')) /* ordinary SET */
    {
      /*  FAM vertraegt nur einen wert hinter der spec.  */
      for (j=0; j<set[i].anz_p; j++) if( getIndex(&set[setNr].pnt,set[setNr].anz_p,set[i].pnt[j]) >-1)
      {
        if( pnt[set[i].pnt[j]].name!=(char *)NULL) fprintf (handle1, " SETA %s p %s \n", set[i].name, pnt[set[i].pnt[j]].name);
      }
      for (j=0; j<set[i].anz_l; j++) if( getIndex(&set[setNr].line,set[setNr].anz_l,set[i].line[j]) >-1)
      {
        if( line[set[i].line[j]].name!=(char *)NULL) fprintf (handle1, " SETA %s l %s \n", set[i].name, line[set[i].line[j]].name);
      }
      for (j=0; j<set[i].anz_c; j++) if( getIndex(&set[setNr].lcmb,set[setNr].anz_c,set[i].lcmb[j]) >-1)
      {
        if( lcmb[set[i].lcmb[j]].name!=(char *)NULL) fprintf (handle1, " SETA %s c %s \n", set[i].name, lcmb[set[i].lcmb[j]].name);
      }
      for (j=0; j<set[i].anz_s; j++) if( getIndex(&set[setNr].surf,set[setNr].anz_s,set[i].surf[j]) >-1)
      {
        if( surf[set[i].surf[j]].name!=(char *)NULL) fprintf (handle1, " SETA %s s %s \n", set[i].name, surf[set[i].surf[j]].name);
      }
      for (j=0; j<set[i].anz_b; j++) if( getIndex(&set[setNr].body,set[setNr].anz_b,set[i].body[j]) >-1)
      {
        if( body[set[i].body[j]].name!=(char *)NULL) fprintf (handle1, " SETA %s b %s \n", set[i].name, body[set[i].body[j]].name);
      }
      for (j=0; j<set[i].anz_nurs; j++) if( getIndex(&set[setNr].nurs,set[setNr].anz_nurs,set[i].nurs[j]) >-1)
      {
        if( nurbs[set[i].nurs[j]].name!=(char *)NULL) fprintf (handle1, " SETA %s S %s \n", set[i].name, nurbs[set[i].nurs[j]].name);
      }
      for (j=0; j<set[i].anz_sh; j++) if( getIndex(&set[setNr].shp,set[setNr].anz_sh,set[i].shp[j]) >-1)
      {
        if( shape[set[i].shp[j]].name!=(char *)NULL) fprintf (handle1, " SETA %s sh %s \n", set[i].name, shape[set[i].shp[j]].name);
      }
    }
  }
  if( compare( setname, "all", 3)==3)
  {
    for (i=0; i<anz->sets; i++)
    {
      if(set[i].name!=(char *)NULL)
      {
        if(set[i].eattr==1)  attr='R';
        else if(set[i].eattr==2)  attr='I';
        else if(set[i].eattr==-1)  attr='U';
        else if(set[i].eattr==-2)  attr='G';
        else attr=' ';
        if(set[i].etyp==1) fprintf (handle1, " ELTY %s HE8%c \n", set[i].name, attr);
        if(set[i].etyp==4) fprintf (handle1, " ELTY %s HE20%c \n", set[i].name, attr);
        if(set[i].etyp==6) fprintf (handle1, " ELTY %s TE10 \n", set[i].name);
        if(set[i].etyp==7) fprintf (handle1, " ELTY %s TR3%c \n", set[i].name, attr);
        if(set[i].etyp==8) fprintf (handle1, " ELTY %s TR6 \n", set[i].name);
        if(set[i].etyp==9) fprintf (handle1, " ELTY %s QU4 \n", set[i].name);
        if(set[i].etyp==10) fprintf (handle1, " ELTY %s QU8 \n", set[i].name);
        if(set[i].etyp==11) fprintf (handle1, " ELTY %s BE2 \n", set[i].name);
        if(set[i].etyp==12) fprintf (handle1, " ELTY %s BE3 \n", set[i].name);
      }
    }
  }
  else
  { 
    if(set[setNr].eattr==1)  attr='R';
    else if(set[setNr].eattr==2)  attr='I';
    else if(set[setNr].eattr==-1)  attr='U';
    else if(set[setNr].eattr==-2)  attr='G';
    else attr=' ';
    if(set[setNr].etyp==1) fprintf (handle1, " ELTY %s HE8%c \n", set[setNr].name, attr);
    if(set[setNr].etyp==4) fprintf (handle1, " ELTY %s HE20%c \n", set[setNr].name, attr);
    if(set[setNr].etyp==6) fprintf (handle1, " ELTY %s TE10 \n", set[setNr].name);
    if(set[setNr].etyp==7) fprintf (handle1, " ELTY %s TR3%c \n", set[setNr].name, attr);
    if(set[setNr].etyp==8) fprintf (handle1, " ELTY %s TR6 \n", set[setNr].name);
    if(set[setNr].etyp==9) fprintf (handle1, " ELTY %s QU4 \n", set[setNr].name);
    if(set[setNr].etyp==10) fprintf (handle1, " ELTY %s QU8 \n", set[setNr].name);
    if(set[setNr].etyp==11) fprintf (handle1, " ELTY %s BE2 \n", set[setNr].name);
    if(set[setNr].etyp==12) fprintf (handle1, " ELTY %s BE3 \n", set[setNr].name);
  }

  /* write the quality thresholds */
  if(eqal.aspr) fprintf (handle1, " EQAL aspr %lf\n", eqal.aspr);
  if(eqal.jbir) fprintf (handle1, " EQAL jbir %lf\n", eqal.jbir);
  if(eqal.mca)  fprintf (handle1, " EQAL mca %lf\n", eqal.mca);

  /* check if specialset->tmp was created */
  i= getSetNr(specialset->tmp);
  if(i>-1)
  { 
    printf (" WARNING: lines in set:%s use seqences which do not include the corner-points. Try to merg the points\n", specialset->tmp);
  }

  /* write the plot,plus commands */
  for (j=0; j<anzGeo->psets; j++ )
  {
    if((pset[j].type[0]!='n')&&(pset[j].type[0]!='e')&&(pset[j].type[0]!='f'))
    {
      if(j==0) fprintf (handle1, " plot %s %s %c\n", pset[j].type, set[pset[j].nr].name, entity_k[pset[j].col]);
      else fprintf (handle1, " plus %s %s %c\n", pset[j].type, set[pset[j].nr].name, entity_k[pset[j].col]);
    }
  }

  fprintf (handle1, "\n");
  fclose(handle1);
  return (1);
}





