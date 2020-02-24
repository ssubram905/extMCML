/***********************************************************
 *  Copyright Univ. of Texas M.D. Anderson Cancer Center
 *  1992.
 *
 *	Launch, move, and record photon weight.
 ****/

#include "mcml.h"

#define STANDARDTEST 0
  /* testing program using fixed rnd seed. */

#define PARTIALREFLECTION 0
  /* 1=split photon, 0=statistical reflection. */

#define COSZERO (1.0-1.0E-12)
  /* cosine of about 1e-6 rad. */

#define COS90D  1.0E-6
  /* cosine of about 1.57 - 1e-6 rad. */

/***********************************************************
 *	A random number generator from Numerical Recipes in C.
 ****/
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC 1.0E-9

void Display(PhotonStruct * Photon_Ptr, InputStruct * In_Ptr);
void Transmission(InputStruct  *	In_Ptr,
				  PhotonStruct *	Photon_Ptr, double ni, double nt, double *ux1, double *uy1, double *uz1);
double RFresnel(double n1,	/* incident refractive index.*/
				double n2,	/* transmit refractive index.*/
				double ca1,	/* cosine of the incident */
							/* angle. 0<a1<90 degrees. */
				double * ca2_Ptr);

float ran3(int *idum)
{
  static int inext,inextp;
  static long ma[56];
  static int iff=0;
  long mj,mk;
  int i,ii,k;

  if (*idum < 0 || iff == 0) {
    iff=1;
    mj=MSEED-(*idum < 0 ? -*idum : *idum);
    mj %= MBIG;
    ma[55]=mj;
    mk=1;
    for (i=1;i<=54;i++) {
      ii=(21*i) % 55;
      ma[ii]=mk;
      mk=mj-mk;
      if (mk < MZ) mk += MBIG;
      mj=ma[ii];
    }
    for (k=1;k<=4;k++)
      for (i=1;i<=55;i++) {
	ma[i] -= ma[1+(i+30) % 55];
	if (ma[i] < MZ) ma[i] += MBIG;
      }
    inext=0;
    inextp=31;
    *idum=1;
  }
  if (++inext == 56) inext=1;
  if (++inextp == 56) inextp=1;
  mj=ma[inext]-ma[inextp];
  if (mj < MZ) mj += MBIG;
  ma[inext]=mj;
  return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC


/***********************************************************
 *	Generate a random number between 0 and 1.  Take a
 *	number as seed the first time entering the function.
 *	The seed is limited to 1<<15.
 *	We found that when idum is too large, ran3 may return
 *	numbers beyond 0 and 1.
 ****/
double RandomNum(void)
{
  static Boolean first_time=1;
  static int idum;	/* seed for ran3. */

  if(first_time) {
#if STANDARDTEST /* Use fixed seed to test the program. */
    idum = - 1;
#else
    idum = -(int)time(NULL)%(1<<15);
	  /* use 16-bit integer as the seed. */
#endif
    ran3(&idum);
    first_time = 0;
    idum = 1;
  }

  return( (double)ran3(&idum) );
}

/***********************************************************
 *	Compute the specular reflection.
 *
 *	If the first layer is a turbid medium, use the Fresnel
 *	reflection from the boundary of the first layer as the
 *	specular reflectance.
 *
 *	If the first layer is glass, multiple reflections in
 *	the first layer is considered to get the specular
 *	reflectance.
 *
 *	The subroutine assumes the Layerspecs array is correctly
 *	initialized.
 ****/
double Rspecular(InputStruct * In_Ptr, LayerStruct * Layerspecs_Ptr, PhotonStruct * Photon_Ptr)
{
  double ni = In_Ptr->layerspecs[0].n;
  double nt = In_Ptr->layerspecs[1].n;
  //printf("Ni and Nt: %lf, %lf",ni,nt);
  double ca1 = Photon_Ptr->ca1;
  double ux1,uz1,uy1,ca2;
  double r;
  r = RFresnel(ni, nt, ca1, &ca2);
  Photon_Ptr->bd == -1;
  Transmission(In_Ptr, Photon_Ptr, ni, nt, &ux1, &uy1, &uz1);
  Photon_Ptr->ux = ux1;
  Photon_Ptr->uy = uy1;
  Photon_Ptr->uz = uz1;
  return (r);
}

/*****************************************
Cosine Intensity Distribution with phase
******************************************/
void cosine(InputStruct * In_Ptr, PhotonStruct *Photon_Ptr)
{
  double cost, sint;	/* cosine and sine of the */
						/* polar deflection angle theta. */
  double cosp, sinp;	/* cosine and sine of the */
						/* azimuthal angle psi. */
  double ux = Photon_Ptr->ux;
  double uy = Photon_Ptr->uy;
  double uz = Photon_Ptr->uz;
  double psi;
  double r;
  double h = In_Ptr->h;
  double theta = asin(In_Ptr->r/(In_Ptr->r+h));

  cost = asin(RandomNum());

  if(cost > 1) cost = 1;
  sint = sqrt(1.0 - cost*cost);  	/* sqrt() is faster than sin(). */

  r = h*sint/cost;
  if(cost - cos(theta) < COS90D )
  {
      Photon_Ptr->dead = 1;
      return;
  }
  psi = 2.0*PI*RandomNum(); /* spin psi 0-2pi. */
  cosp = cos(psi);
  if(psi<PI)
    sinp = sqrt(1.0 - cosp*cosp);
	  /* sqrt() is faster than sin(). */
  else
    sinp = - sqrt(1.0 - cosp*cosp);


  Photon_Ptr->ux = sint*cosp;
  Photon_Ptr->uy = sint*sinp;
  Photon_Ptr->uz = cost;

  Photon_Ptr->x+= r*cosp;
  Photon_Ptr->y+= r*sinp;
  Photon_Ptr->z+= In_Ptr->r - sqrt(In_Ptr->r*In_Ptr->r - Photon_Ptr->y*Photon_Ptr->y);
 // printf("\n\nPhoton position: %lf, %lf, %lf", Photon_Ptr->x, Photon_Ptr->y, Photon_Ptr->z);

}


/***********************************************************
 *	Initialize a photon packet.
 ****/
void LaunchPhoton(InputStruct  * In_Ptr,
				  LayerStruct  * Layerspecs_Ptr,
				  PhotonStruct * Photon_Ptr)
{
  Photon_Ptr->dead 	= 0;
  Photon_Ptr->layer = 1;
  Photon_Ptr->s	= 0;
  Photon_Ptr->sleft= 0;

  Photon_Ptr->x 	= 0.0;
  Photon_Ptr->y	 	= 0.0;
  Photon_Ptr->z	 	= 0.0;
  Photon_Ptr->ux	= 0.0;
  Photon_Ptr->uy	= 0.0;
  Photon_Ptr->uz	= 1.0;
  Photon_Ptr->bd	= 0.0;
  Photon_Ptr->ca1	= 1.0;
  In_Ptr->Rsp_onerun = 0.0;

  switch (In_Ptr->source_type)
  {
      case 1: Photon_Ptr->ca1 = 1.0;
              break;

      case 2: cosine(In_Ptr, Photon_Ptr);
              break;

  }
  if(Photon_Ptr->dead == 1)
  {
      return;
  }
  In_Ptr->Rsp_onerun = Rspecular(In_Ptr, Layerspecs_Ptr, Photon_Ptr);
  Photon_Ptr->w	 	= 1.0 - In_Ptr->Rsp_onerun;
  /*if((Layerspecs_Ptr[1].mua == 0.0)&& (Layerspecs_Ptr[1].mus == 0.0))
  { //glass layer.
    Photon_Ptr->layer 	= 2;
    Photon_Ptr->z	= Layerspecs_Ptr[2].z0;
  }*/
}
/***********************************************************
 *	Choose (sample) a new theta angle for photon propagation
 *	according to the anisotropy.
 *
 *	If anisotropy g is 0, then
 *		cos(theta) = 2*rand-1.
 *	otherwise
 *		sample according to the Gegenbauer-Kernel function.
 *
 *	Returns the cosine of the polar deflection angle theta.
 ****/
double SpinTheta(double g, double alpha)
{
  double cost;

  if(g == 0.0)
    cost = 2*RandomNum() -1;
  else
  {
    if(alpha == 0.5)
    {
        double temp = (1-g*g)/(1-g+2*g*RandomNum());
        cost = (1+g*g - temp*temp)/(2*g);
    }
    else
    {
        double temp1 = pow((1-g),2*alpha);
        double temp2 = pow((1+g),2*alpha);
        double temp3 = 1.0/temp2 +RandomNum()*(temp2-temp1)/(temp1*temp2);
        double temp = pow(temp3,(-1.0/alpha));
        cost = (1+g*g - temp)/(2*g);
    }
  }
    if(cost < -1) cost = -1;
	else if(cost > 1) cost = 1;

  return(cost);
}


/***********************************************************
 *	Choose a new direction for photon propagation by
 *	sampling the polar deflection angle theta and the
 *	azimuthal angle psi.
 *
 *	Note:
 *  	theta: 0 - pi so sin(theta) is always positive
 *  	feel free to use sqrt() for cos(theta).
 *
 *  	psi:   0 - 2pi
 *  	for 0-pi  sin(psi) is +
 *  	for pi-2pi sin(psi) is -
 ****/
void Spin(double g, double alpha,
		  PhotonStruct * Photon_Ptr)
{
  double cost, sint;	/* cosine and sine of the */
						/* polar deflection angle theta. */
  double cosp, sinp;	/* cosine and sine of the */
						/* azimuthal angle psi. */
  double ux = Photon_Ptr->ux;
  double uy = Photon_Ptr->uy;
  double uz = Photon_Ptr->uz;
  double psi;

  cost = SpinTheta(g, alpha);
  sint = sqrt(1.0 - cost*cost);
	/* sqrt() is faster than sin(). */

  psi = 2.0*PI*RandomNum(); /* spin psi 0-2pi. */
  cosp = cos(psi);
  if(psi<PI)
    sinp = sqrt(1.0 - cosp*cosp);
	  /* sqrt() is faster than sin(). */
  else
    sinp = - sqrt(1.0 - cosp*cosp);

  if(fabs(uz) > COSZERO)  { 	/* normal incident. */
    Photon_Ptr->ux = sint*cosp;
    Photon_Ptr->uy = sint*sinp;
    Photon_Ptr->uz = cost*SIGN(uz);
	  /* SIGN() is faster than division. */
  }
  else  {		/* regular incident. */
    double temp = sqrt(1.0 - uz*uz);
    Photon_Ptr->ux = sint*(ux*uz*cosp - uy*sinp)
					/temp + ux*cost;
    Photon_Ptr->uy = sint*(uy*uz*cosp + ux*sinp)
					/temp + uy*cost;
    Photon_Ptr->uz = -sint*cosp*temp + uz*cost;
  }
}

/***********************************************************
 *	Move the photon s away in the current layer of medium.
 ****/
void Hop(PhotonStruct *	Photon_Ptr)
{
  double s = Photon_Ptr->s;
  Photon_Ptr->prevx = Photon_Ptr->x;
  Photon_Ptr->prevy = Photon_Ptr->y;
  Photon_Ptr->prevz = Photon_Ptr->z;

  Photon_Ptr->x += s*Photon_Ptr->ux;
  Photon_Ptr->y += s*Photon_Ptr->uy;
  Photon_Ptr->z += s*Photon_Ptr->uz;
}

/***********************************************************
 *	If uz != 0, return the photon step size in glass,
 *	Otherwise, return 0.
 *
 *	The step size is the distance between the current
 *	position and the boundary in the photon direction.
 *
 *	Make sure uz !=0 before calling this function.
 ****/
void StepSizeInGlass(PhotonStruct *  Photon_Ptr,
					 InputStruct  *  In_Ptr)
{
    Photon_Ptr->s = In_Ptr->dz*10;
}

/***********************************************************
 *	Pick a step size for a photon packet when it is in
 *	tissue.
 *	If the member sleft is zero, make a new step size
 *	with: -log(rnd)/(mua+mus).
 *	Otherwise, pick up the leftover in sleft.
 *
 *	Layer is the index to layer.
 *	In_Ptr is the input parameters.
 ****/
void StepSizeInTissue(PhotonStruct * Photon_Ptr,
					  InputStruct  * In_Ptr)
{
  short  layer = Photon_Ptr->layer;
  double mua = In_Ptr->layerspecs[layer].mua;
  double mus = In_Ptr->layerspecs[layer].mus;

  if(Photon_Ptr->sleft == 0.0) {  /* make a new step. */
    double rnd;

    do rnd = RandomNum();
      while( rnd <= 0.0 );    /* avoid zero. */
	Photon_Ptr->s = -log(rnd)/(mua+mus);
  }
  else {	/* take the leftover. */
	Photon_Ptr->s = Photon_Ptr->sleft/(mua+mus);
	Photon_Ptr->sleft = 0.0;
  }
}

/***********************************************************
 *	Check if the step will hit the boundary.
 *	Return 1 if hit boundary.
 *	Return 0 otherwise.
 *
 * 	If the projected step hits the boundary, the members
 *	s and sleft of Photon_Ptr are updated.
 ****/
Boolean HitBoundary(PhotonStruct *  Photon_Ptr,
					InputStruct  *  In_Ptr)
{
  double dx2,dx1;  /* dx1-distance between new and old positions .
                      dx2-distance between center and new position.  */
  long double dl_b;
  double temp;
  short  layer = Photon_Ptr->layer;
  double rv[3], inv[3];

  double s = Photon_Ptr->s;
  double ux = Photon_Ptr->ux;
  double uz = Photon_Ptr->uz;
  double uy = Photon_Ptr->uy;
  Boolean hit;
  double r = In_Ptr->r;
 // printf("Radius: %lf", r);
  double y2,z2,y1,z1;
  double y_int,z_int;
  double r0 = In_Ptr->layerspecs[layer].r0;  // radius of top layer
  double r1 = In_Ptr->layerspecs[layer].r1;  // radius of bottom layer

  z1 = Photon_Ptr->z; // old photon position
  y1 = Photon_Ptr->y;

  /* Distance to the boundary. */
  z2= z1 + s*Photon_Ptr->uz;   // new photon position
  y2= y1 + s*Photon_Ptr->uy;
  dx2= sqrt((z2-r)*(z2-r)+(y2)*(y2));   // Distance of new point from center

  double d[2],f[2], e[2];
  d[0] = y2 - y1;   // Vector in the direction of the photon propagation
  d[1] = z2 - z1;
  f[0] = y1 - 0.0;  // vector from center of circle to photon position
  f[1] = z1 - r;
  e[0] = y1;        // previous photon position
  e[1] = z1;
  /*Solving the parameteric equation of a line intersecting
    a circle in the y-z plane with parameter t. */

  if (dx2>r0) // distance greater than top layer radius, photon is crossing up.
  {  double a = Dot(d ,d,2) ;
     double b = 2*Dot(f, d,2 ) ;
     double c = Dot( f, f,2 ) - r0*r0 ;
     double discriminant = b*b-4*a*c;
     if( discriminant < 0 ){
      hit =0;
      return hit;                  // no intersection
     }
     else
     {  /* ray didn't totally miss sphere, so there is a solution to
           the equation.*/
       discriminant = sqrt( discriminant );
       double t1 = (-b + discriminant)/(2*a);
       double t2 = (-b - discriminant)/(2*a);
       //printf("\n t1 and t2 (distance greater than r0): %lf, %lf:", t1 , t2);
       if( t1 >= 0 && t1 <= 1 )
       { dl_b=t1;
         hit =1;
       }
       else  if( t2 >= 0 && t2 <= 1 )
       { dl_b = t2;
        hit = 1;
       }
       else {
       hit = 0;
       return hit;
      }
     }
     Photon_Ptr->bd = 1;
  }
  else if (dx2<r1) // distance less than top layer radius, photon is crossing down
  {  double a = Dot( d , d ,2) ;
     double b = 2*Dot(f, d ,2) ;
     double c = Dot( f, f ,2) - r1*r1 ;
     double discriminant = b*b-4*a*c;
     //printf("\n a,b,c and disc (distance less than r1): %lf ,%lf ,%lf, %lf:", a,b,c,discriminant);
     if( discriminant < 0 ){
      hit =0;
      return hit;                  // no intersection
     }
     else
     { /* ray didn't totally miss sphere, so there is a solution to
           the equation.*/
       discriminant = sqrt( discriminant );
       double t1 = (-b + discriminant)/(2*a);
       double t2 = (-b - discriminant)/(2*a);
       //printf("\n t1 and t2 (distance less than r1): %lf, %lf:", t1 , t2);
       if( t1 >= 0 && t1 <= 1 )
       { dl_b = t1;
        hit = 1;
       }
       else if( t2 >= 0 && t2 <= 1 )
      { dl_b = t2;
        hit = 1;
       }
      else {
      hit = 0;
      return hit;
      }
     }
     Photon_Ptr->bd = -1;
  }
  else
  {
     hit =0;
     return hit;
  }

  y_int = y1 + dl_b*d[0];
  z_int = z1 + dl_b*d[1];
  if((z_int-z1==0)&&(y_int-y1==0)) // To avoid atan2() exception.
  {
      hit = 0;
      return hit;
  }
  dl_b = sqrt((y_int-y1)*(y_int-y1)+(z_int-z1)*(z_int-z1));
  double mut = In_Ptr->layerspecs[layer].mua+ In_Ptr->layerspecs[layer].mus;
  Photon_Ptr->sleft = (Photon_Ptr->s - dl_b)*mut;
  Photon_Ptr->s    = dl_b;
  hit = 1;


  //Update photon positions.
  Photon_Ptr->z = Photon_Ptr->z + Photon_Ptr->s*Photon_Ptr->uz;
  Photon_Ptr->y = Photon_Ptr->y + Photon_Ptr->s*Photon_Ptr->uy;
  Photon_Ptr->x = Photon_Ptr->x + Photon_Ptr->s*Photon_Ptr->ux;

  rv[0] = 0.0;
  rv[1] = Photon_Ptr->y - 0.0;                 // unit vector in the direction of the radius
  rv[2] = Photon_Ptr->z-In_Ptr->r;
  temp = sqrt(rv[1]*rv[1]+rv[2]*rv[2]);
  rv[1] /= temp;
  rv[2] /= temp;

  inv[0] =  Photon_Ptr->ux;  // input vector
  inv[1] =  Photon_Ptr->uy;
  inv[2] =  Photon_Ptr->uz;

  Photon_Ptr->ca1 = fabs(Dot(rv,inv,3));

  return hit;
}

/***********************************************************
 *	Drop photon weight inside the tissue (not glass).
 *
 *  The photon is assumed not dead.
 *
 *	The weight drop is dw = w*mua/(mua+mus).
 *
 *	The dropped weight is assigned to the absorption array
 *	elements.
 ****/
void Drop(InputStruct  *	In_Ptr,
		  PhotonStruct *	Photon_Ptr,
		  OutStruct *		Out_Ptr)
{
  double dwa;		/* absorbed weight.*/
  double x = Photon_Ptr->x;
  double y = Photon_Ptr->y;
  double z = Photon_Ptr->z;
  short izd, ixd, iyd;	/* LW 5/20/98. To avoid out of short range.*/
  short  iz, ix, iy;	/* index to z & r. */
  short  layer = Photon_Ptr->layer;
  double mua, mus;
  long ind;
  short nr = In_Ptr->nr;
  short nz = In_Ptr->nz;
  short  r = In_Ptr->nz/2;
  /* update photon weight. */
  mua = In_Ptr->layerspecs[layer].mua;
  mus = In_Ptr->layerspecs[layer].mus;
  dwa = Photon_Ptr->w * mua/(mua+mus);
  Photon_Ptr->w -= dwa;
  /* compute array indices. */
  izd = fabs(Photon_Ptr->z)/In_Ptr->dz;
  if(izd>In_Ptr->nz-1) iz=In_Ptr->nz-1;
  else iz = izd;
  ixd = fabs(x)/In_Ptr->dr;
  if(ixd>In_Ptr->nr/2-1) 
    {ix=In_Ptr->nr/2-1;
     Photon_Ptr->dead = 1;
    }
  else   ix= ixd;
  iyd = fabs(y)/In_Ptr->dz;
  if(iyd>sqrt(2*r*izd-izd*izd))
  { iyd=sqrt(2*r*izd-izd*izd);
    dwa=0;
  }
  if(iyd>In_Ptr->nz/2-1) iy=In_Ptr->nz/2-1;
  else iy = iyd;
  
  //if((ix < 0 || iy < 0 || iz <0) || (ix > nr/2-1 || iy > nz/2-1 || iz > nz -1))  return;

   /* assign dwa to the absorption array element. */
  if(x >= 0 && y>= 0)  Out_Ptr->A_xyz[nr/2+ix][nz/2+iy][iz]	+= dwa/mua;
  else if(x < 0 && y>= 0)  Out_Ptr->A_xyz[nr/2-1-ix][nz/2+iy][iz]	+= dwa/mua;
  else if(x < 0 && y< 0)  Out_Ptr->A_xyz[nr/2-1-ix][nz/2-1-iy][iz]	+= dwa/mua;
  else if(x >= 0 && y< 0)  Out_Ptr->A_xyz[nr/2+ix][nz/2-1-iy][iz]	+= dwa/mua;
}

/***********************************************************
 *	The photon weight is small, and the photon packet tries
 *	to survive a roulette.
 ****/
void Roulette(PhotonStruct * Photon_Ptr)
{
  if(Photon_Ptr->w == 0.0)
    Photon_Ptr->dead = 1;
  else if(RandomNum() < CHANCE) /* survived the roulette.*/
    Photon_Ptr->w /= CHANCE;
  else
    Photon_Ptr->dead = 1;
}

/***********************************************************
 *	Compute the Fresnel reflectance.
 *
 *	Make sure that the cosine of the incident angle a1
 *	is positive, and the case when the angle is greater
 *	than the critical angle is ruled out.
 *
 * 	Avoid trigonometric function operations as much as
 *	possible, because they are computation-intensive.
 ****/
double RFresnel(double n1,	/* incident refractive index.*/
				double n2,	/* transmit refractive index.*/
				double ca1,	/* cosine of the incident */
							/* angle. 0<a1<90 degrees. */
				double * ca2_Ptr)  /* pointer to the */
							/* cosine of the transmission */
							/* angle. a2>0. */
{
  double r;

  if(n1==n2) {			  	/** matched boundary. **/
    *ca2_Ptr = ca1;
    r = 0.0;
  }
  else if(ca1>COSZERO) {	/** normal incident. **/
    *ca2_Ptr = ca1;
    r = (n2-n1)/(n2+n1);
    r *= r;
  }
  else if(ca1<COS90D)  {	/** very slant. **/
    *ca2_Ptr = 0.0;
    r = 1.0;
  }
  else  {			  		/** general. **/
    double sa1, sa2;
	  /* sine of the incident and transmission angles. */
    double ca2;

    sa1 = sqrt(1-ca1*ca1);
    sa2 = n1*sa1/n2;
    if(sa2>=1.0) {
	  /* double check for total internal reflection. */
      *ca2_Ptr = 0.0;
      r = 1.0;
    }
    else  {
      double cap, cam;	/* cosines of the sum ap or */
						/* difference am of the two */
						/* angles. ap = a1+a2 */
						/* am = a1 - a2. */
      double sap, sam;	/* sines. */

      *ca2_Ptr = ca2 = sqrt(1-sa2*sa2);

      cap = ca1*ca2 - sa1*sa2; /* c+ = cc - ss. */
      cam = ca1*ca2 + sa1*sa2; /* c- = cc + ss. */
      sap = sa1*ca2 + ca1*sa2; /* s+ = sc + cs. */
      sam = sa1*ca2 - ca1*sa2; /* s- = sc - cs. */
      r = 0.5*sam*sam*(cam*cam+cap*cap)/(sap*sap*cam*cam);
		/* rearranged for speed. */
    }
  }
  //printf("Reflection: %lf", r);
  return(r);
}

/***********************************************************
 *	Record the photon weight exiting the first layer(uz<0),
 *	no matter whether the layer is glass or not, to the
 *	reflection array.
 *
 *	Update the photon weight as well.
 ****/
void RecordR(double			Refl,	/* reflectance. */
			 InputStruct  *	In_Ptr,
			 PhotonStruct *	Photon_Ptr,
			 OutStruct *	Out_Ptr)
{
  double x = Photon_Ptr->x;
  double y = Photon_Ptr->y;
  double z = Photon_Ptr->z;
  double r = In_Ptr->r;
  short nr = In_Ptr->nr;
  short nz = In_Ptr->nz;

  short  ix, iy;	/* index to r & angle. */
  double ixd, iyd, izd;	/* LW 5/20/98. To avoid out of short range.*/
  if (y > fabs(2*r*z) - z*z) y =  fabs(2*r*z) - z*z;
  izd = fabs(Photon_Ptr->z)/In_Ptr->dz;

  ixd = fabs(x)/In_Ptr->dr;
  if(ixd>In_Ptr->nr/2-1) ix=In_Ptr->nr/2-1;
  else   ix= ixd;

  iyd = fabs(y)/In_Ptr->dz;
  if(iyd>sqrt(2*r*(izd)-izd*izd)) iyd=sqrt(2*r*(izd)-izd*izd);
  if(iyd>In_Ptr->nz/2-1) iy=In_Ptr->nz/2-1;
  else iy = iyd;

  /* assign photon to the reflection array element. */
  if(x >= 0 && y>= 0)  Out_Ptr->Rd_xy[nr/2+ix][nz/2+iy]	+= Photon_Ptr->w*(1.0-Refl);
  else if(x < 0 && y>= 0)  Out_Ptr->Rd_xy[nr/2-1-ix][nz/2+iy]	+= Photon_Ptr->w*(1.0-Refl);
  else if(x < 0 && y< 0)  Out_Ptr->Rd_xy[nr/2-1-ix][nz/2-1-iy]	+= Photon_Ptr->w*(1.0-Refl);
  else if(x >= 0 && y< 0)  Out_Ptr->Rd_xy[nr/2+ix][nz/2-1-iy]	+= Photon_Ptr->w*(1.0-Refl);

  Photon_Ptr->w *= Refl;
}

/***********************************************************
 *	Record the photon weight exiting the last layer(uz>0),
 *	no matter whether the layer is glass or not, to the
 *	transmittance array.
 *
 *	Update the photon weight as well.
 ****/
void RecordT(double 		Refl,
			 InputStruct  *	In_Ptr,
			 PhotonStruct *	Photon_Ptr,
			 OutStruct *	Out_Ptr)
{
  double x = Photon_Ptr->x;
  double y = Photon_Ptr->y;
  double z = Photon_Ptr->z;
  double r = In_Ptr->r;
  short nr = In_Ptr->nr;
  short nz = In_Ptr->nz;

  short  ix, iy;	/* index to r & angle. */
  double ixd, iyd, izd;	/* LW 5/20/98. To avoid out of short range.*/
  if (y > fabs(2*r*z) - z*z) y =  fabs(2*r*z) - z*z;
  izd = fabs(Photon_Ptr->z)/In_Ptr->dz;

  ixd = fabs(x)/In_Ptr->dr;
  if(ixd>In_Ptr->nr/2-1) ix=In_Ptr->nr/2-1;
  else   ix= ixd;

  iyd = fabs(y)/In_Ptr->dz;
  if(iyd>sqrt(2*r*(izd)-izd*izd)) iyd=sqrt(2*r*(izd)-izd*izd);
  if(iyd>In_Ptr->nz/2-1) iy=In_Ptr->nz/2-1;
  else iy = iyd;

  /* assign photon to the transmittance array element. */
  if(x >= 0 && y>= 0)  Out_Ptr->Tt_xy[nr/2+ix][nz/2+iy]	+= Photon_Ptr->w*(1.0-Refl);
  else if(x < 0 && y>= 0)  Out_Ptr->Tt_xy[nr/2-1-ix][nz/2+iy]	+= Photon_Ptr->w*(1.0-Refl);
  else if(x < 0 && y< 0)  Out_Ptr->Tt_xy[nr/2-1-ix][nz/2-1-iy]	+= Photon_Ptr->w*(1.0-Refl);
  else if(x >= 0 && y< 0)  Out_Ptr->Tt_xy[nr/2+ix][nz/2-1-iy]	+= Photon_Ptr->w*(1.0-Refl);


  Photon_Ptr->w *= Refl;
}
/***********************************************************
 *	Decide the reflected uz and uy with respect to original
 *  co-ordinates. Reflection is found by mirroring the input
 *  photon position about the radius of artery layer.
 ****/
void Reflection(InputStruct  *	In_Ptr,
				  PhotonStruct *	Photon_Ptr, double *ux1, double *uy1, double *uz1)
{
        double r[3], inv[3];
        double temp;

        r[0] = 0.0;
        r[1] = Photon_Ptr->y - 0.0;                 // unit vector in the direction of the radius
        r[2] = Photon_Ptr->z-In_Ptr->r;
        temp = sqrt(r[1]*r[1]+r[2]*r[2]);
        r[1] /= temp;
        r[2] /= temp;
        //printf("\n\nDirection cosine, Radius: %lf, %lf, %lf", r[0], r[1], r[2]);

        inv[0] =  Photon_Ptr->ux;  // input vector
        inv[1] =  Photon_Ptr->uy;
        inv[2] =  Photon_Ptr->uz;

        //printf("\nDirection cosine, Input: %lf, %lf, %lf", Photon_Ptr->ux, Photon_Ptr->uy, Photon_Ptr-> uz);

        *ux1 = -2*Dot(r,inv,3)*r[0] + inv[0];
        *uy1 = -2*Dot(r,inv,3)*r[1] + inv[1];
        *uz1 = -2*Dot(r,inv,3)*r[2] + inv[2];

        //printf("\nDirection cosines, refl: %lf , %lf , %lf ", *ux1, *uy1, *uz1);
        return;

}
/***********************************************************
 *	Decide the transmitted uz and uy with respect to original
 *  co-ordinates. Transmitted directions are found by finding
 *  an unit vector with angle ca2 w.r.t the radius of the circle
 *  at the point of crossing.
 ****/
void Transmission(InputStruct  *	In_Ptr,
				  PhotonStruct *	Photon_Ptr, double n1, double n2, double *ux1, double *uy1, double *uz1)
{
        double r[3], inv[3];
        double temp, ca1;

        r[0] = 0.0;
        r[1] = Photon_Ptr->y - 0.0;                 // unit vector in the direction of the radius
        r[2] = Photon_Ptr->z-In_Ptr->r;
        temp = sqrt(r[1]*r[1]+r[2]*r[2]);
        r[1] /= temp;
        r[2] /= temp;
        if(Photon_Ptr->bd == 1)
        { r[1] *= -1;
          r[2] *= -1;
        }
        /*else if(Photon_Ptr->bd == 1)
        { r[1] *= 1;
          r[2] *= 1;
        }
        else
        { *ux1 = Photon_Ptr->ux;
          *uy1 = Photon_Ptr->uy;
          *uz1 = Photon_Ptr->uz;
          return;
        }*/

        //printf("\n\nDirection cosines, Radius: %lf , %lf , %lf ", r[0], r[1], r[2]);

        inv[0] =  Photon_Ptr->ux;  // input vector
        inv[1] =  Photon_Ptr->uy;
        inv[2] =  Photon_Ptr->uz;

        //printf("\nDirection cosines, Input: %lf , %lf , %lf ", inv[0], inv[1], inv[2]);
        if(n1==n2)
        {
            *ux1 = Photon_Ptr->ux;
            *uy1 = Photon_Ptr->uy;
            *uz1 = Photon_Ptr->uz;
            //printf("\nDirection cosines, trans: %lf , %lf , %lf ", *ux1, *uy1, *uz1);
            return;
        }
        else
        {
            ca1 = -Dot(r,inv,3);
            //printf("\n CA1, trans: %lf ", ca1);
            *ux1 = n1/n2*inv[0]+ r[0]*(n1/n2*ca1 - sqrt(1- (n1*n1)/(n2*n2)*(1 - ca1*ca1)));
            *uy1 = n1/n2*inv[1]+ r[1]*(n1/n2*ca1 - sqrt(1- (n1*n1)/(n2*n2)*(1 - ca1*ca1)));  //direction cosines of transmitted vector
            *uz1 = n1/n2*inv[2]+ r[2]*(n1/n2*ca1 - sqrt(1- (n1*n1)/(n2*n2)*(1 - ca1*ca1)));
            //printf("\nDirection cosines, trans: %lf , %lf , %lf ", *ux1, *uy1, *uz1);
            return;
        }
}
/***********************************************************
 *	Decide whether the photon will be transmitted or
 *	reflected on the upper boundary (uz<0) of the current
 *	layer.
 *
 *	If "layer" is the first layer, the photon packet will
 *	be partially transmitted and partially reflected if
 *	PARTIALREFLECTION is set to 1,
 *	or the photon packet will be either transmitted or
 *	reflected determined statistically if PARTIALREFLECTION
 *	is set to 0.
 *
 *	Record the transmitted photon weight as reflection.
 *
 *	If the "layer" is not the first layer and the photon
 *	packet is transmitted, move the photon to "layer-1".
 *
 *	Update the photon parmameters.
 ****/
void CrossUpOrNot(InputStruct  *	In_Ptr,
				  PhotonStruct *	Photon_Ptr,
				  OutStruct *		Out_Ptr)
{
  double ca1 = Photon_Ptr->ca1; /* z directional cosine. */
  double ca2;	/* cosines of transmission alpha. always */
				/* positive. */
  double r=0.0;	/* reflectance */
  short  layer = Photon_Ptr->layer;
  double ni = In_Ptr->layerspecs[layer].n;
  double nt = In_Ptr->layerspecs[layer-1].n;
  double ux1, uz1, uy1;
  double theta = asin(0.05/In_Ptr->r);

  /* Get r. */
  if( ca1 <= In_Ptr->layerspecs[layer].cos_crit0)
    r=1.0;		      /* total internal reflection. */
  else r = RFresnel(ni, nt, ca1, &ca2);


#if PARTIALREFLECTION
  if(layer == 1 && r<1.0) {	/* partially transmitted. */
    Transmission(In_Ptr, Photon_Ptr, ni, nt, &ux1, &uy1, &uz1);
    Photon_Ptr->ux = ux1;
    Photon_Ptr->uy = uy1;
    Photon_Ptr->uz = uz1;
                           /* transmitted photon. */
    if (Photon_Ptr->z <= In_Ptr->r*cos(theta))  RecordR(r, In_Ptr, Photon_Ptr, Out_Ptr);
    else if (Photon_Ptr->z >= In_Ptr->r + In_Ptr->r*(cos(theta)))
    RecordT(r, In_Ptr, Photon_Ptr, Out_Ptr);

    Reflection(In_Ptr, Photon_Ptr, &ux1, &uy1, &uz1);
    Photon_Ptr->ux = ux1;
    Photon_Ptr->uy = uy1;
    Photon_Ptr->uz = uz1;

    if (Photon_Ptr->z <= In_Ptr->r*(1-cos(theta)))  RecordR(r, In_Ptr, Photon_Ptr, Out_Ptr);
    else if (Photon_Ptr->z >= In_Ptr->r + In_Ptr->r*(cos(theta)))
    RecordT(r, In_Ptr, Photon_Ptr, Out_Ptr);

  }
  else if(RandomNum() > r) {/* transmitted to layer-1. */
    Photon_Ptr->layer--;
    Transmission(In_Ptr, Photon_Ptr, ni, nt, &ux1, &uy1, &uz1);
    Photon_Ptr->ux = ux1;
    Photon_Ptr->uy = uy1;
    Photon_Ptr->uz = uz1;
      }
  else			      		/* reflected. */
   {    Reflection(In_Ptr, Photon_Ptr, &ux1, &uy1, &uz1);
        Photon_Ptr->ux = ux1;
        Photon_Ptr->uy = uy1;
        Photon_Ptr->uz = uz1;


   }
#else
  if(RandomNum() > r) {		/* transmitted to layer-1. */
    if(layer==1)  {
      Transmission(In_Ptr, Photon_Ptr, ni, nt, &ux1, &uy1, &uz1);
      Photon_Ptr->ux = ux1;
      Photon_Ptr->uy = uy1;
      Photon_Ptr->uz = uz1;

      if (Photon_Ptr->z <= In_Ptr->r*(1-cos(theta))) {
      RecordR(0.0, In_Ptr, Photon_Ptr, Out_Ptr);
      }
      else if (Photon_Ptr->z >= In_Ptr->r + In_Ptr->r*(cos(theta))){
      RecordT(0.0, In_Ptr, Photon_Ptr, Out_Ptr);
      }
      Photon_Ptr->dead = 1;
    }
    else {
      //printf("Photon Transmitted");
      Photon_Ptr->layer--;
      Transmission(In_Ptr, Photon_Ptr, ni, nt, &ux1, &uy1, &uz1);
      Photon_Ptr->ux = ux1;
      Photon_Ptr->uy = uy1;
      Photon_Ptr->uz = uz1;
       }
  }
  else			      		/* reflected. */
   {    //printf("Photon reflected:");
        Reflection(In_Ptr, Photon_Ptr, &ux1, &uy1, &uz1);
        Photon_Ptr->ux = ux1;
        Photon_Ptr->uy = uy1;
        Photon_Ptr->uz = uz1;


   }
#endif
}

/***********************************************************
 *	Decide whether the photon will be transmitted  or be
 *	reflected on the bottom boundary (uz>0) of the current
 *	layer.
 *
 *	If the photon is transmitted, move the photon to
 *	"layer+1". If "layer" is the last layer, record the
 *	transmitted weight as transmittance. See comments for
 *	CrossUpOrNot.
 *
 *	Update the photon parmameters.
 ****/
void CrossDnOrNot(InputStruct  *	In_Ptr,
				  PhotonStruct *	Photon_Ptr,
				  OutStruct *		Out_Ptr)
{
  double ca1 = Photon_Ptr->ca1; /* z directional cosine. */
  double ca2;	/* cosines of transmission alpha. */
  double r=0.0;	/* reflectance */
  short  layer = Photon_Ptr->layer;
  double ni = In_Ptr->layerspecs[layer].n;
  double nt = In_Ptr->layerspecs[layer+1].n;
  double ux1, uz1, uy1;

  /* Get r. */
  if( ca1 <= In_Ptr->layerspecs[layer].cos_crit1)
    r=1.0;		/* total internal reflection. */
  else r = RFresnel(ni, nt, ca1, &ca2);

#if PARTIALREFLECTION
  if(layer == In_Ptr->num_layers && r<1.0) {
    Transmission(In_Ptr, Photon_Ptr, ni, nt, &ux1, &uy1, &uz1);
    Photon_Ptr->ux = ux1;
    Photon_Ptr->uy = uy1;
    Photon_Ptr->uz = uz1;
  }
  else if(RandomNum() > r) {/* transmitted to layer+1. */
    if(layer == In_Ptr->num_layers) Photon_Ptr->layer--;
    else Photon_Ptr->layer++;
    Transmission(In_Ptr, Photon_Ptr, ni, nt, &ux1, &uy1, &uz1);
    Photon_Ptr->ux = ux1;
    Photon_Ptr->uy = uy1;
    Photon_Ptr->uz = uz1;
 }
   else			      		/* reflected. */
   {    Reflection(In_Ptr, Photon_Ptr, &ux1, &uy1, &uz1);
        Photon_Ptr->ux = ux1;
        Photon_Ptr->uy = uy1;
        Photon_Ptr->uz = uz1;


   }
#else
  if(RandomNum() > r)
  {		/* transmitted to layer+1. */
     if(layer == In_Ptr->num_layers) Photon_Ptr->layer--;
      else Photon_Ptr->layer++;
      Transmission(In_Ptr, Photon_Ptr, ni, nt, &ux1, &uy1, &uz1);
      Photon_Ptr->ux = ux1;
      Photon_Ptr->uy = uy1;
      Photon_Ptr->uz = uz1;
  }
   else			      		/* reflected. */
   {    //printf("\nPhoton reflected:");
        Reflection(In_Ptr, Photon_Ptr, &ux1, &uy1, &uz1);
        Photon_Ptr->ux = ux1;
        Photon_Ptr->uy = uy1;
        Photon_Ptr->uz = uz1;
   }
#endif
}

/***********************************************************
 ****/
void CrossOrNot(InputStruct  *	In_Ptr,
				PhotonStruct *	Photon_Ptr,
				OutStruct    *	Out_Ptr)
{
  if(Photon_Ptr->bd == 1)
    CrossUpOrNot(In_Ptr, Photon_Ptr, Out_Ptr);
  else if (Photon_Ptr->bd == -1)
    CrossDnOrNot(In_Ptr, Photon_Ptr, Out_Ptr);
}

/***********************************************************
 *	Move the photon packet in glass layer.
 *	Horizontal photons are killed because they will
 *	never interact with tissue again.
 ****/
void HopInGlass(InputStruct  * In_Ptr,
				PhotonStruct * Photon_Ptr,
				OutStruct    * Out_Ptr)
{
  double dl;     /* step size. 1/cm */

  if(Photon_Ptr->uz == 0.0) {
	/* horizontal photon in glass is killed. */
    Photon_Ptr->dead = 1;
  }
  else {
    StepSizeInGlass(Photon_Ptr, In_Ptr);
    while(!(HitBoundary(Photon_Ptr, In_Ptr))) Hop(Photon_Ptr);
    CrossOrNot(In_Ptr, Photon_Ptr, Out_Ptr);
    Photon_Ptr->sleft = 0.0;
  }
}

/***********************************************************
 *	Set a step size, move the photon, drop some weight,
 *	choose a new photon direction for propagation.
 *
 *	When a step size is long enough for the photon to
 *	hit an interface, this step is divided into two steps.
 *	First, move the photon to the boundary free of
 *	absorption or scattering, then decide whether the
 *	photon is reflected or transmitted.
 *	Then move the photon in the current or transmission
 *	medium with the unfinished stepsize to interaction
 *	site.  If the unfinished stepsize is still too long,
 *	repeat the above process.
 ****/
void HopDropSpinInTissue(InputStruct  *  In_Ptr,
						 PhotonStruct *  Photon_Ptr,
						 OutStruct    *  Out_Ptr)
{
  StepSizeInTissue(Photon_Ptr, In_Ptr);

  if(HitBoundary(Photon_Ptr, In_Ptr)) {
      //printf("\nHit boundary Display");
      //Display(Photon_Ptr , In_Ptr);
     CrossOrNot(In_Ptr, Photon_Ptr, Out_Ptr);

  }
  else {
    //printf("\nNo hit boundary Display");
    //Display(Photon_Ptr , In_Ptr);
    Hop(Photon_Ptr);
    Drop(In_Ptr, Photon_Ptr, Out_Ptr);
    Spin(In_Ptr->layerspecs[Photon_Ptr->layer].g, In_Ptr->layerspecs[Photon_Ptr->layer].alpha,
		Photon_Ptr);
    Photon_Ptr->ca1 = 0.0;


  }
}

/***********************************************************
 ****/
void HopDropSpin(InputStruct  *  In_Ptr,
				 PhotonStruct *  Photon_Ptr,
				 OutStruct    *  Out_Ptr)
{
  short layer = Photon_Ptr->layer;
  if((In_Ptr->layerspecs[layer].mua == 0.0)
  && (In_Ptr->layerspecs[layer].mus == 0.0))
	/* glass layer. */
    HopInGlass(In_Ptr, Photon_Ptr, Out_Ptr);
  else
   {       //printf("Press ENTER to continue: ");
           //while (1)
           //{ int c=getchar();
           //if (c=='\n' || c==EOF) break; }
       HopDropSpinInTissue(In_Ptr, Photon_Ptr, Out_Ptr);
   }


  if( Photon_Ptr->w < In_Ptr->Wth && !Photon_Ptr->dead)
    Roulette(Photon_Ptr);
}

void Display(PhotonStruct * Photon_Ptr, InputStruct * In_Ptr)
{   short layer = Photon_Ptr->layer;

    printf("\nPhoton position(x,y,z) : %lf , %lf, %lf  ", Photon_Ptr->x, Photon_Ptr->y, Photon_Ptr->z);
    printf("\nPhoton direction(ux,uy,uz): %lf, %lf, %lf ", Photon_Ptr->ux, Photon_Ptr->uy, Photon_Ptr->uz);
    printf("\nIncident angle (w.r.t radius): %lf ", Photon_Ptr->ca1);
    printf("\nPhoton Weight: %lf ", Photon_Ptr->w);
    printf("\nPhoton layer: %hd ", Photon_Ptr->layer);
    printf("\nPhoton step size: %lf ", Photon_Ptr->s);
    printf("\nPhoton sleft: %lf ", Photon_Ptr->sleft);
    printf("\nBoundary Hit? %hd", Photon_Ptr->bd);
    printf("\n\n Layer Properties ");
    printf("\nLayer mua, mus, g : %lf, %lf, %lf", In_Ptr->layerspecs[layer].mua, In_Ptr->layerspecs[layer].mus, In_Ptr->layerspecs[layer].g);
    printf("\nLayer r0, r1 : %lf, %lf", In_Ptr->layerspecs[layer].r0, In_Ptr->layerspecs[layer].r1);

}
