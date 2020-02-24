/***********************************************************
 *  Copyright Univ. of Texas M.D. Anderson Cancer Center
 *  1992.
 *
 *	Some routines modified from Numerical Recipes in C,
 *	including error report, array or matrix declaration
 *	and releasing.
 ****/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>

/***********************************************************
 *	Report error message to stderr, then exit the program
 *	with signal 1.
 ****/
void nrerror(char error_text[])

{
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}

/***********************************************************
 *	Allocate an array with index from nl to nh inclusive.
 *
 *	Original matrix and vector from Numerical Recipes in C
 *	don't initialize the elements to zero. This will
 *	be accomplished by the following functions.
 ****/
double *AllocVector(short nl, short nh)
{
  double *v;
  short i;

  v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
  if (!v) nrerror("allocation failure in vector()");

  v -= nl;
  for(i=nl;i<=nh;i++) v[i] = 0.0;	/* init. */
  return v;
}

/***********************************************************
 *	Allocate a matrix with row index from nrl to nrh
 *	inclusive, and column index from ncl to nch
 *	inclusive.
 ****/
double **AllocMatrix(short nrl,short nrh,
					 short ncl,short nch)
{
  short i,j;
  double **m;

  m=(double **) malloc((unsigned) (nrh-nrl+1)
						*sizeof(double*));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m -= nrl;

  for(i=nrl;i<=nrh;i++) {
    m[i]=(double *) malloc((unsigned) (nch-ncl+1)
						*sizeof(double));
    if (!m[i]) nrerror("allocation failure 2 in matrix()");
    m[i] -= ncl;
  }

  for(i=nrl;i<=nrh;i++)
    for(j=ncl;j<=nch;j++) m[i][j] = 0.0;
  return m;
}

/***********************************************************
 *	Allocate a 3D array with x index from nxl to nxh
 *	inclusive, y index from nyl to nyh and z index from nzl to nzh
 *	inclusive.
 ****/
double ***Alloc3D(short nxl,short nxh,
					 short nyl,short nyh,
					     short nzl, short nzh)
{
  double ***t;
  short i,j,k;


  t=(double ***) malloc((unsigned) (nxh-nxl+1)*sizeof(double **));
  if (!t) nrerror("allocation failure 1 in matrix()");
  t -= nxl;

  for(i=nxl;i<=nxh;i++) {
    t[i]=(double **) malloc((unsigned) (nyh-nyl+1)*sizeof(double *));
    if (!t[i]) nrerror("allocation failure 2 in matrix()");
    t[i] -= nyl;
    for(j=nyl;j<=nyh;j++) {
         t[i][j]=(double *) malloc((unsigned) (nzh-nzl+1)*sizeof(double));
         if (!t[i][j]) nrerror("allocation failure 3 in matrix()");
         t[i][j] -= nzl;}

  }

  for(i=nxl;i<=nxh;i++)
    for(j=nyl;j<=nyh;j++)
       for(k=nzl; k<=nzh;k++) t[i][j][k] = 0.0;
  return t;
}
/***********************************************************
 *Index to 3D array.
 ****/
long Getindex(int x, int y, int z, int Size)
{
     return (Size*Size*x + Size*y + z);
}
/***********************************************************
 *	Release the memory.
 ****/
void FreeVector(double *v,short nl,short nh)
{
  free((char*) (v+nl));
}

/***********************************************************
 *	Release the memory.
 ****/
void FreeMatrix(double **m,short nrl,short nrh,
			    short ncl,short nch)
{
  short i;

  for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
  free((char*) (m+nrl));
}


/***********************************************************
 *	Release the memory.
 ****/
void Free3D(double ***t,short nxl,short nxh,
			    short nyl,short nyh, short nzl, short nzh)
{
  short i,j;

  for(i=nxh;i>=nxl;i--)
   {for(j=nyl;j>=nyl;j--) free((char*) (t[i][j]+nzl));
   }
  for(i=nxh;i>=nxl;i--) free((char*) (t[i]+nyl));
  free((char*) (t+nxl));
}
/***********************************************************
 *	function to calculate dot product
 ****/

double Dot(double dot1[], double dot2[], short l)
  {  if(l==2)
        return (dot1[0]*dot2[0] + dot1[1]*dot2[1]);
     else
        return (dot1[0]*dot2[0] + dot1[1]*dot2[1] + dot1[2]*dot2[2]);
  }
