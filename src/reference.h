#include <malloc.h>
#include <stdlib.h>

#ifdef _WIN32
#define alloca _alloca
#endif

#define dbl_prec 53
#define binSize 45

///////////////////////////////// copied from CAMPARY package/////////////////////////////
inline double FPadd_rn(const double x, const double y){
  return x + y;
}
inline double FPmul_rn(const double x, const double y){
  return x * y;
}
inline double fma_d_rn_cpu(const double x, const double y, double xy) {
  double C,px,qx,hx,py,qy,hy,lx,ly,fma;
  /*@ assert C == 0x1p27+1; */
	C = 0x1p27+1;

  px = x*C;
  qx = x-px;
  hx = px+qx;
  lx = x-hx;

  py = y*C;
  qy = y-py;
  hy = py+qy;
  ly = y-hy;

  fma = -x*y+hx*hy;
  fma += hx*ly;
  fma += hy*lx;
  fma += lx*ly;
  return fma;
}

inline double FPfma_rn(const double x, const double y, const double z){
	#ifdef FP_FAST_FMA
  //#warning cpu has fma
	  return fma(x, y, z);
	#else
   	return fma_d_rn_cpu(x, y, z);
	#endif
}
/* Computes fl(a*b) and err(a*b). */
inline double two_prod(const double a, const double b, double *err){
	const double p = FPmul_rn(a, b);
 	*err = FPfma_rn(a, b, -p);
 	return p;
}
inline double two_sum(const double a, const double b, double *err){
  double s = FPadd_rn(a, b);
  double aa = FPadd_rn(s, -b);
  double bb = FPadd_rn(s, -aa);
  double da = FPadd_rn(a, -aa);
  double db = FPadd_rn(b, -bb);
  *err = FPadd_rn(da, db);
  return s;
}
/* Computes fl(a+b) and err(a+b). Assumes |a| >= |b| */
inline double fast_two_sum(const double a, const double b, double *err){
  double s = FPadd_rn(a, b);
  double z = FPadd_rn(s, -a);
  *err = FPadd_rn(b, -z);
  return s;
}

inline void fast_VecSumErrBranchmul(const double *x, double *r,int sX, int sR){
  int ptr = 0, i = 1;
  double e = x[0];

  while(i<sX && ptr<sR){
    r[ptr] = fast_two_sum(e, x[i], &e); i++;
    if(e == 0.) e = r[ptr]; else ptr++;
  }
  if(ptr<sR && e!=0.){ r[ptr] = e; ptr++; }
  for(i=ptr; i<sR; i++) r[i] = 0.;
}

inline void fast_VecSumErrBranch(double *x, int sX, int sR){
	int ptr = 0, i = 1;
	double e = x[0];

  	while(i<sX && ptr<sR){
		x[ptr] = fast_two_sum(e, x[i], &e); i++;
		if(e == 0.) e = x[ptr]; else ptr++;
  	}
  	if(ptr<sR && e!=0.){ x[ptr] = e; ptr++; }
		for(i=ptr; i<sR; i++) x[i] = 0.;
}

inline void fast_VecSumErr(double *x, int sX){
	double e;
	x[0] = fast_two_sum(x[0], x[1], &e);
	for(int i=2; i<sX-1; i++) x[i-1] = fast_two_sum(e, x[i], &e);
	x[sX-2] = fast_two_sum(e, x[sX-1], &x[sX-1]);
}

inline void merge(double const *x, double const *y, double *z,int K,int L){
  int i=0, j=0;
  for(int n=0; n<K+L; n++){
    if(i==K || (j<L && fabs(y[j])>fabs(x[i]))){ z[n] = y[j]; j++; }
    else{ z[n] = x[i]; i++; }
  }
}

inline void renorm_rand2L(int sX, int sR, double x[]){
	double pr;
	int i, j, ptr = 0;

  x[0] = two_sum(x[0], x[1], &x[1]);
  for(i=2; i<sX; i++){
		pr = two_sum(x[i-1], x[i], &x[i]);
    for(j=i-2; j>0; j--) pr = two_sum(x[j], pr, &x[j+1]);
		x[0] = fast_two_sum(x[0], pr, &x[1]);
	}

	i = 2;
  if(x[1] == 0.) pr = x[0];
  else { pr = x[1]; ptr++;}
  while(ptr<sR && i<sX){
    x[ptr] = fast_two_sum(pr, x[i], &pr); i++;
    if(pr == 0.) pr = x[ptr]; else ptr++;
  }
  if(ptr<sR && pr!=0.){ x[ptr] = pr; ptr++; }
  for(i=ptr; i<sX; i++) x[i] = 0.;
}

// Camapary mul
/** Multiplies x and y and returns the result in r, as an ulp-nonoverlapping expansion.
    K - size of x, L - size of y, R - size of r.
    Computes all partial products based on scalar products and then accumulates
    them in a special designed structure that has a fixed-point flavour.
    double-precision = 53, bin size = 45; **/
// Multiplication_accurate with relative error <= 2^{-(p-1)R}( 1 + (r+1)2^{-p} )
// template <int K, int L, int R>
inline void certifiedMul(int K, int L, int R, const double *x, const double *y, double *r)
{
    int const LBN = R * dbl_prec / binSize + 2;
    double *B = (double *)alloca((LBN + 2) * sizeof(double));
    double *lim = (double *)alloca((LBN + 2) * sizeof(double));
    int i;

    int *exp_x = (int *)alloca((K) * sizeof(int));
    int *exp_y = (int *)alloca((L) * sizeof(int));
    for (i = 0; i < K; i++)
        frexp(x[i], &exp_x[i]);
    for (i = 0; i < L; i++)
        frexp(y[i], &exp_y[i]);

    double factor = ldexp(1.0, -binSize); /* 2^(-45) */
    int exp_start = exp_x[0] + exp_y[0];
    lim[0] = ldexp(1.5, exp_start - binSize + dbl_prec - 1);
    B[0] = lim[0];
    for (i = 1; i < LBN + 2; i++)
    {
        lim[i] = FPmul_rn(lim[i - 1], factor);
        B[i] = lim[i];
    }

#ifdef UNROLL_MUL_MAIN
    unrollPartialProds<L, R, 0, K, 1>::unrollX_all(exp_start, LBN, x, exp_x, y, exp_y, B);
#else
    double p, e;
    int j, sh, l;
    for (i = 0; i < K; i++)
    {
        for (j = 0; j < L; j++)
        {
            l = exp_start - (exp_x[i] + exp_y[j]);
            sh = (int)(l / binSize);
            l = l - sh * binSize;
            if (sh < LBN - 1)
            {
                p = two_prod(x[i], y[j], &e);
                if (l < 30)
                { // binSize - 2*(dbl_prec-binSize-1) - 1){
                    B[sh] = fast_two_sum(B[sh], p, &p);
                    B[sh + 1] = FPadd_rn(B[sh + 1], p);

                    B[sh + 1] = fast_two_sum(B[sh + 1], e, &e);
                    B[sh + 2] = FPadd_rn(B[sh + 2], e);
                }
                else if (l < 37)
                { // binSize - (dbl_prec-binSize-1) - 1){
                    B[sh] = fast_two_sum(B[sh], p, &p);
                    B[sh + 1] = FPadd_rn(B[sh + 1], p);

                    B[sh + 1] = fast_two_sum(B[sh + 1], e, &e);
                    B[sh + 2] = fast_two_sum(B[sh + 2], e, &e);
                    B[sh + 3] = FPadd_rn(B[sh + 3], e);
                }
                else
                {
                    B[sh] = fast_two_sum(B[sh], p, &p);
                    B[sh + 1] = fast_two_sum(B[sh + 1], p, &p);
                    B[sh + 2] = FPadd_rn(B[sh + 2], p);

                    B[sh + 2] = fast_two_sum(B[sh + 2], e, &e);
                    B[sh + 3] = FPadd_rn(B[sh + 3], e);
                }
            }
        }
    }
#endif

    /* unbias the B's */
    for (i = 0; i < LBN; i++)
        B[i] = FPadd_rn(B[i], -lim[i]);
    // recheck if right implementation used if result is not correct
    fast_VecSumErrBranchmul(B, r, LBN, R);
}
