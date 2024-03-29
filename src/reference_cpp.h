///////////////////////////////// copied from CAMPARY package/////////////////////////////

#define TWO_SUM_FLOPS 6
#define FAST_TWO_SUM_FLOPS 3
#define TWO_MULT_FMA_FLOPS 3

// #define UNROLL_RENORM

/* ========== TwoSum ========== */

/* Computes fl(a+b) and err(a+b). */
static inline double two_sum(const double a, const double b, double &err){
  double s = FPadd_rn(a, b);
  double aa = FPadd_rn(s, -b);
  double bb = FPadd_rn(s, -aa);
  double da = FPadd_rn(a, -aa);
  double db = FPadd_rn(b, -bb);
  err = FPadd_rn(da, db);
  return s;
}

/* Computes fl(a+b) and err(a+b). Assumes |a| >= |b| */
static inline double fast_two_sum(const double a, const double b, double &err){
  double s = FPadd_rn(a, b);
  double z = FPadd_rn(s, -a);
  err = FPadd_rn(b, -z);
  return s;
}

/* ========== VecSum ========== */

template<int from, int to, int step=1>
struct unrollVecSum{
  static inline void fast_VecSum(const double *x, double s, double *r){
    s = fast_two_sum(x[from], s, &r[from+1]);
    unrollVecSum<from+step,to,step>::fast_VecSum(x, s, r);
  }
  static inline void VecSum(const double *x, double s, double *r){
    s = two_sum(x[from], s, r[from+1]);
    unrollVecSum<from+step,to,step>::VecSum(x, s, r);
  }
};

// Terminal case
template<int from, int step>
struct unrollVecSum<from,from,step>{
  static inline void fast_VecSum(const double *x, double s, double *r){ r[0] = s; }
  static inline void VecSum(const double *x, double s, double *r){ r[0] = s; }
};

template<int sX>
static inline void fast_VecSum(double *x, double *r, int _sX){
	double s = fast_two_sum(x[sX-2], x[sX-1], &r[sX-1]);
	#ifdef UNROLL_RENORM
		unrollVecSum<sX-3,-1,-1>::fast_VecSum(x, s, r);
	#else
		for(int i=sX-3; i>0; i--) s = fast_two_sum(x[i], s, r[i+1]);
  	r[0] = fast_two_sum(x[0], s, r[1]);
	#endif
}

/* ========== Renormalization ========== */

/** Algorithm for normalizing the array x. The result satisfies |x_i|<ulp(x_{i-1}); P-nonoverlapping expansion.
    sX is the size of array x and sR is the required size for the normalized expansion
    Requirement: sR <= sX **/
template<int sX, int sR>
static inline void fast_renorm3L(double *x, int _sX, double *r, int _sR){
  double e, f[sX];
  int ptr = 0, i = 2;

  fast_VecSum<sX>(x,f,sX);

  if(f[1] == 0.) e = f[0]; else { e = f[1]; ptr++; }
  while(ptr<((sR<sX)?sR+1:sR) && i<sX){
    f[ptr] = fast_two_sum(e, f[i], &e); i++;
    if(e == 0.) e = f[ptr]; else ptr++;
  }
  if(i==sX && ptr<((sR<sX)?sR+1:sR) && e!=0.){ f[ptr] = e; ptr++; }
  for(i=ptr; i<((sR<sX)?sR+1:sR); i++) f[i] = 0.;

  for(int j=0; j<((sR<sX)?sR-1:sX-2); j++){
    r[j] = fast_two_sum(f[j], f[j+1], &e);
    for(int i=j+1; i<((sR<sX)?sR-1:sX-2); i++) f[i] = fast_two_sum(e, f[i+1], &e);
    f[(sR<sX)?sR-1:sX-2] = fast_two_sum(e, f[(sR<sX)?sR:sX-1], &f[(sR<sX)?sR:sX-1]);
  }
  r[(sR<sX)?sR-1:sX-2] = f[(sR<sX)?sR-1:sX-2];
  r[(sR<sX)?sR-1:sX-1] = f[(sR<sX)?sR-1:sX-1];

  // for(i=sR; i<sX; i++) r[i] = 0.;
}

/* ========== Addition ========== */

/** Merges the arrays x of size K and y of size L, and puts the result in the array z
    All 3 arrays x, y and z are sorted in decreasing order of magnitude
		K - size of x, L - size of y **/
template <int K, int L>
static inline void merge(double const *x, double const *y, double *z){
  int i=0, j=0;
  for(int n=0; n<K+L; n++){
    if(i==K || (j<L && fabs(y[j])>fabs(x[i]))){ z[n] = y[j]; j++; }
    else{ z[n] = x[i]; i++; }
  }
}

template<int sX>
static inline void VecSum(const double *x, double *r){
  double s = fast_two_sum(x[sX-2], x[sX-1], r[sX-1]);
	#ifdef UNROLL_RENORM
  	unrollVecSum<sX-3,-1,-1>::VecSum(x, s, r);
	#else
	  for(int i=sX-3; i>0; i--) s = two_sum(x[i], s, r[i+1]);
  	r[0] = two_sum(x[0], s, r[1]);
	#endif
}

/** VecSum for addition of a FP expansion with a FP number **/
template<int sX>
static inline void VecSum_4Add1(const double *x, const double y, double *r){
  double s = two_sum(x[sX-1], y, r[sX]);
	#ifdef UNROLL_RENORM
  	unrollVecSum<sX-2,-1,-1>::VecSum(x, s, r);
	#else
	  for(int i=sX-2; i>0; i--) s = two_sum(x[i], s, r[i+1]);
  	r[0] = two_sum(x[0], s, r[1]);
	#endif
}


/** Algorithm for normalizing the array x. The result satisfies |x_i|<=ulp(x_{i-1}); ulp-nonoverlapping expansion.
    sX is the size of array x and sR is the required size for the normalized expansion
    Requirement: sR <= sX **/
// Renormalize
template<int sX, int sR>
static inline void renorm2L(const double *x, double *r){
  double e, f[sX];
  int ptr = 0, i = 2;

  VecSum<sX>(x,f);

  if(f[1] == 0.) e = f[0];
  else { r[0] = f[0]; e = f[1]; ptr++; }
  while(ptr<sR && i<sX){
    r[ptr] = fast_two_sum(e, f[i], e); i++;
    if(e == 0.) e = r[ptr]; else ptr++;
  }
  if(ptr<sR && e!=0.){ r[ptr] = e; ptr++; }
  for(i=ptr; i<sR; i++) r[i] = 0.;
}

/** Addition of a FP expansion with a FP number **/
// Renormalize - special case for adding an expansion with a FP number
template<int sX, int sR>
static inline void renorm2L_4Add1(const double *x, const double y, double *r){
  double e, f[sX+1];
  int ptr = 0, i = 2;

  VecSum_4Add1<sX>(x,y,f);

  if(f[1] == 0.) e = f[0];
  else { r[0] = f[0]; e = f[1]; ptr++; }
  while(ptr<sR && i<sX+1){
    r[ptr] = fast_two_sum(e, f[i], e); i++;
    if(e == 0.) e = r[ptr]; else ptr++;
  }
  if(ptr<sR && e!=0.){ r[ptr] = e; ptr++; }
  for(i=ptr; i<sR; i++) r[i] = 0.;
}

/** Adds x and y and returns the result in r, as an ulp-nonoverlapping expansion. 
		K - size of x, L - size of y, R - size of r **/
// Addition_accurate with relative error < 4.5 * 2^{-(prec-1)R}
template <int K, int L, int R>
static inline void certifiedAdd(double *x, double *y, double *r, int _K, int _L, int _R){
  if(K == 1) renorm2L_4Add1<L,R>(y, x[0], r);
  else if(L == 1) renorm2L_4Add1<K,R>(x, y[0], r);
  else{
    double f[K+L];
    merge<K,L>( x, y, f );
    // renorm2L<K+L,R>( f, r );
    fast_renorm3L<K+L,R>(f, K+L, r, R);
  }
}

/* ========== Multiplication ========== */

/* Computes fl(a*b) and err(a*b). */
static inline double two_prod(const double a, const double b, double &err){
	const double p = FPmul_rn(a, b);
 	err = FPfma_rn(a, b, -p);
 	return p;
}

/** Algorithm for normalizing the array x that contains random numbers. 
    After the first level	the result satisfies |x_i|<uls(x_{i-1}); S-nonoverlapping expansion.
		In the end, the result satisfies |x_i|<=ulp(x_{i-1}); ulp-nonoverlapping expansion.
    sX is the size of array x and sR is the required size for the normalized expansion
    Requirement: sR <= sX **/
// Renormalize_random
template<int sX, int sR>
static inline void renorm_rand2L(const double *x, double *r){
  double pr, f[sX];
  int i, j, ptr = 0;

  f[0] = two_sum(x[0], x[1], f[1]);
  for(i=2; i<sX; i++){
    pr = two_sum(f[i-1], x[i], f[i]);
    for(j=i-2; j>0; j--) pr = two_sum(f[j], pr, f[j+1]);
    f[0] = fast_two_sum(f[0], pr, f[1]);
  }

  i = 2;
  if(f[1] == 0.) pr = f[0];
  else { r[0] = f[0]; pr = f[1]; ptr++; }
  while(ptr<sR && i<sX){
    r[ptr] = fast_two_sum(pr, f[i], pr); i++;
    if(pr == 0.) pr = r[ptr]; else ptr++;
  }
  if(ptr<sR && pr!=0.){ r[ptr] = pr; ptr++; }
  for(i=ptr; i<sR; i++) r[i] = 0.;
}

/** Multiplies x and y and returns the result in r
    K - size of x, L - size of y, R - size of z. Constraints: K>=L.
    Uses a generalization of Bailey's multiplication algorithm that simulates the paper-and-pencil method
		Constraints on the algorithm: K>=L; R<=2*K*L (if R>=2*K*L the extra terms will be 0) **/
// Multiplication_quick with relative error <= 2^{-(p-1)(R+1)} (128/127 (K+L)- 129/254 R - 385/254 + 2^{p-1} + 2^{-p-r}(R^2+R)((R+1)!)^2 )
template <int K, int L, int R>
static inline void baileyMul_renorm(double *x, double *y, double *z, int _K, int _L, int _R){
  double p, e, f[R+1];
  int i, j, n = R;
  for(i=K+L-1; i<R; i++) f[i] = 0.0;

  //computes the last term of the result f[R]
  //computes the necessary products (if any) using just simple multiplication
  const int nn = R;
  if (nn<L){
    f[nn] = FPmul_rn(x[0], y[nn]);
    for(i=1; i<=nn; i++){
      #if defined __CUDA_ARCH__ || FP_FAST_FMA
        f[nn] = FPfma_rn(x[i], y[nn-i], f[nn]);
      #else
        f[nn] = FPadd_rn(f[nn], FPmul_rn(x[i], y[nn-i]));
      #endif
    }
  }else if (nn<K){
    f[nn] = FPmul_rn(x[nn], y[0]);
    for(j=1; j<L; j++){
      #if defined __CUDA_ARCH__ || FP_FAST_FMA
        f[nn] = FPfma_rn(x[nn-j], y[j], f[nn]);
      #else
        f[nn] = FPadd_rn(f[nn], FPmul_rn(x[nn-j], y[j]));
      #endif
    }
  }else if (nn<K+L-1){
    f[nn] = FPmul_rn(x[nn-L+1], y[L-1]);
    for(i=nn-L+2; i<K; i++){
      #if defined __CUDA_ARCH__ || FP_FAST_FMA
        f[nn] = FPfma_rn(x[i], y[nn-i], f[nn]);
      #else
        f[nn] = FPadd_rn(f[nn], FPmul_rn(x[i], y[nn-i]));
      #endif
    }
  }else
    f[nn] = 0.0;

  // computes the last R-K elements of the result
  // we will have (K+L-1 - n) products to compute and sum
  for (n=(R+1>K+L-1)?K+L-2:R-1; n>=K; n--){
    f[n] = two_prod(x[n-L+1], y[L-1], e);
    for(j=n+1; j<R; j++) f[j] = two_sum(f[j], e, e);
    f[R] = FPadd_rn(f[R], e);

    for(i=n-L+2; i<K; i++){
      p = two_prod(x[i],y[n-i],e);
      for(j=n+1; j<R; j++) f[j] = two_sum(f[j], e, e);
      f[R] = FPadd_rn(f[R], e);

      f[n] = two_sum(f[n], p, e);
      for(j=n+1; j<R; j++) f[j] = two_sum(f[j], e, e);
      f[R] = FPadd_rn(f[R], e);
    }
  }
  // computes the elements with the same number of products inside
  // we will have L products to compute and sum
  for (n=(R+1>K)?K-1:R-1; n>=L; n--){
    f[n] = two_prod(x[n], y[0], e);
    for(i=n+1; i<R; i++) f[i] = two_sum(f[i], e, e);
    f[R] = FPadd_rn(f[R], e);

    for(j=1; j<L; j++){
      p = two_prod(x[n-j], y[j], e);
      for(i=n+1; i<R; i++) f[i] = two_sum(f[i], e, e);
      f[R] = FPadd_rn(f[R], e);

      f[n] = two_sum(f[n], p, e);
      for(i=n+1; i<R; i++) f[i] = two_sum(f[i], e, e);
      f[R] = FPadd_rn(f[R], e);
    }
  }
  // computes the first L doubles of the result
  // we will have (n+1) prducts to compute and sum
  for (n=(R+1>L)?L-1:R-1; n>=0; n--){
    f[n] = two_prod(x[0], y[n], e);
    for(j=n+1; j<R; j++) f[j] = two_sum(f[j], e, e);
    f[R] = FPadd_rn(f[R], e);

    for(i=1; i<=n; i++){
      p = two_prod(x[i], y[n-i], e);
      for(j=n+1; j<R; j++) f[j] = two_sum(f[j], e, e);
      f[R] = FPadd_rn(f[R], e);

      f[n] = two_sum(f[n], p, e);
      for(j=n+1; j<R; j++) f[j] = two_sum(f[j], e, e);
      f[R] = FPadd_rn(f[R], e);
    }
  }
  // renorm_rand2L<R+1,R>(f, z);	
  fast_renorm3L<R+1,R>(f, R+1, z, R);
}

// Calculate flop count for multiplication
template <int K, int L, int R>
unsigned int calculate_multiplication_flops(){
  int i, j, n = R;

  int flops = 0;

  //computes the last term of the result f[R]
  //computes the necessary products (if any) using just simple multiplication
  const int nn = R;
  if (nn<L){
    flops += 1;
    for(i=1; i<=nn; i++){
      #if defined __CUDA_ARCH__ || FP_FAST_FMA
        flops += 2;
      #else
        flops += 2;
      #endif
    }
  }else if (nn<K){
    flops += 1;
    for(j=1; j<L; j++){
      #if defined __CUDA_ARCH__ || FP_FAST_FMA
        flops += 2;
      #else
        flops += 2;
      #endif
    }
  }else if (nn<K+L-1){
    flops += 1;
    for(i=nn-L+2; i<K; i++){
      #if defined __CUDA_ARCH__ || FP_FAST_FMA
        flops += 2;
      #else
        flops += 2;
      #endif
    }
  }

  // computes the last R-K elements of the result
  // we will have (K+L-1 - n) products to compute and sum
  for (n=(R+1>K+L-1)?K+L-2:R-1; n>=K; n--){
    flops += TWO_MULT_FMA_FLOPS;
    for(j=n+1; j<R; j++) 
    {
      flops += TWO_SUM_FLOPS;
    }
    flops += 1;

    for(i=n-L+2; i<K; i++){
      flops += TWO_MULT_FMA_FLOPS;
      for(j=n+1; j<R; j++) 
      {
        flops += TWO_SUM_FLOPS;
      }
      flops += 1;

      flops += TWO_SUM_FLOPS;
      for(j=n+1; j<R; j++) 
      {
        flops += TWO_SUM_FLOPS;
      }
      flops += 1;
    }
  }
  // computes the elements with the same number of products inside
  // we will have L products to compute and sum
  for (n=(R+1>K)?K-1:R-1; n>=L; n--){
    flops += TWO_MULT_FMA_FLOPS;
    for(i=n+1; i<R; i++) 
    {
      flops += TWO_SUM_FLOPS;
    }
    flops += 1;

    for(j=1; j<L; j++){
      flops += TWO_MULT_FMA_FLOPS;
      for(i=n+1; i<R; i++) 
      {
        flops += TWO_SUM_FLOPS;
      }
      flops += 1;

      flops += TWO_SUM_FLOPS;
      for(i=n+1; i<R; i++) 
      {
        flops += TWO_SUM_FLOPS;
      }
      flops += 1;
    }
  }
  // computes the first L doubles of the result
  // we will have (n+1) prducts to compute and sum
  for (n=(R+1>L)?L-1:R-1; n>=0; n--){
    flops += TWO_MULT_FMA_FLOPS;
    for(j=n+1; j<R; j++) 
    {
      flops += TWO_SUM_FLOPS;
    }
    flops += 1;

    for(i=1; i<=n; i++){
      flops += TWO_MULT_FMA_FLOPS;
      for(j=n+1; j<R; j++) 
      {
        flops += TWO_SUM_FLOPS;
      }
      flops += 1;

      flops += TWO_SUM_FLOPS;
      for(j=n+1; j<R; j++) 
      {
        flops += TWO_SUM_FLOPS;
      }
      flops += 1;
    }
  }

  return flops;
}

template<int sX, int sR>
static inline void fast_VecSumErrBranch(const double *x, double *r){
  int ptr = 0, i = 1;
  double e = x[0];

  while(i<sX && ptr<sR){
    r[ptr] = fast_two_sum(e, x[i], e); i++;
    if(e == 0.) e = r[ptr]; else ptr++;
  }
  if(ptr<sR && e!=0.){ r[ptr] = e; ptr++; }
  for(i=ptr; i<sR; i++) r[i] = 0.;
}

/** Multiplies x and y and returns the result in r, as an ulp-nonoverlapping expansion. 
    K - size of x, L - size of y, R - size of z. Constraints: K>=L.
    The algorithm computes the partial products in a paper-and-pencil fashion and then 
    accumulates them in a special designed structure that has a fixed-point flavour.
		double-precision = 53, bin size = 45;
    For operations using double-double, triple-double and quad-double we provide specialized 
    versions that use a generalization of the QD library's multiplication algorithm. **/
// Multiplication_accurate with relative error <= 2^{-(p-1)R}( 1 + (r+1)2^{-p} + 2^{-p+2}(K+L-R-2) )
template <int K, int L, int R>
static inline void truncatedMul(double *x, double *y, double *r, int _K, int _L, int _R){
	int const LBN = R*dbl_prec/binSize + 2;
  double B[LBN+2], lim[LBN+2];
  int i;

	int exp_x[(R+1<K)?R+1:K], exp_y[(R+1<L)?R+1:L];
	for(i=0; i<((R+1<K)?R+1:K); i++) frexp(x[i], &exp_x[i]);
	for(i=0; i<((R+1<L)?R+1:L); i++) frexp(y[i], &exp_y[i]);

	double factor = ldexp(1.0, -binSize); // 2^(-45)
	int exp_start = exp_x[0] + exp_y[0];
	lim[0] = ldexp(1.5, exp_start - binSize + dbl_prec-1); B[0] = lim[0];
	for(i=1; i<LBN+2; i++){ lim[i] = FPmul_rn(lim[i-1], factor); B[i] = lim[i]; }

	#ifdef UNROLL_MUL_MAIN
	  unrollPartialProds< L, R, 0, (R<K?R:K), 1 > :: unrollX_trunc(exp_start, LBN, x, exp_x, y, exp_y, B);
	#else
		int j, l, sh;
		double p, e;
		for(i=0; i<(R<K?R:K); i++){
			for(j=0; j<(R-i<L?R-i:L); j++){
				l  = exp_start - (exp_x[i]+exp_y[j]);
				sh = (int)(l/binSize);
		  	l  = l - sh*binSize;
			  if(sh < LBN-1){
					p = two_prod(x[i], y[j], e);
			    if(l < 30){ // binSize - 2*(dbl_prec-binSize-1) - 1){
						B[sh] = fast_two_sum(B[sh], p, p);
						B[sh+1] = FPadd_rn(B[sh+1], p);

						B[sh+1] = fast_two_sum(B[sh+1], e, e);
						B[sh+2] = FPadd_rn(B[sh+2], e);
					}else if(l < 37){ // binSize - (dbl_prec-binSize-1) - 1){
						B[sh] = fast_two_sum(B[sh], p, p);
		        B[sh+1] = FPadd_rn(B[sh+1], p);

						B[sh+1] = fast_two_sum(B[sh+1], e, e);
						B[sh+2] = fast_two_sum(B[sh+2], e, e);
		        B[sh+3] = FPadd_rn(B[sh+3], e);
					}else{
						B[sh] = fast_two_sum(B[sh], p, p);
						B[sh+1] = fast_two_sum(B[sh+1], p, p);
		        B[sh+2] = FPadd_rn(B[sh+2], p);

						B[sh+2] = fast_two_sum(B[sh+2], e, e);
		        B[sh+3] = FPadd_rn(B[sh+3], e);
		}}}}
	#endif
	
	//computation of the error correction terms; using just simple multiplication
  if (R < L){
  	#ifdef UNROLL_MUL_ERROR
  		unrollPartialProds< L, R, 0, R+1, 1 > :: unroll_simple1(exp_start, LBN, x, exp_x, y, exp_y, B);
  	#else
  		double p;
  		int sh;
			for(i=0; i<=R; i++){
		  	sh = (int)((exp_start - (exp_x[i]+exp_y[R-i])) / binSize);
		    if(sh < LBN){
					p = FPmul_rn(x[i], y[R-i]);
					B[sh] = fast_two_sum(B[sh], p, p);
		      B[sh+1] = fast_two_sum(B[sh+1], p, p);
		      B[sh+2] = FPadd_rn(B[sh+2], p);
		  }}
    #endif
	}else if(R < K){
  	#ifdef UNROLL_MUL_ERROR
  		unrollPartialProds< L, R, 0, L, 1 > :: unroll_simple2(exp_start, LBN, x, exp_x, y, exp_y, B);
  	#else
  		double p;
  		int sh;
		  for(i=0; i<L; i++){
				sh = (int)((exp_start - (exp_x[R-i]+exp_y[i])) / binSize);
		    if(sh < LBN){
		      p = FPmul_rn(x[R-i], y[i]);
					B[sh] = fast_two_sum(B[sh], p, p);
		      B[sh+1] = fast_two_sum(B[sh+1], p, p);
		      B[sh+2] = FPadd_rn(B[sh+2], p);
		  }}
		#endif
  }else if(R < K+L-1){
  	#ifdef UNROLL_MUL_ERROR
			unrollPartialProds< L, R, R-L+1, K, 1 > :: unroll_simple1(exp_start, LBN, x, exp_x, y, exp_y, B);
  	#else
  		double p;
  		int sh;
		  for(i=R-L+1; i<K; i++){
				sh = (int)((exp_start - (exp_x[i]+exp_y[R-i])) / binSize);
		    if(sh < LBN){
		      p = FPmul_rn(x[i], y[R-i]);
					B[sh] = fast_two_sum(B[sh], p, p);
		      B[sh+1] = fast_two_sum(B[sh+1], p, p);
		      B[sh+2] = FPadd_rn(B[sh+2], p);
			}}
		#endif
	}

	/* unbias the B's */
	for (i=0; i<LBN; i++) B[i] = FPadd_rn(B[i], -lim[i]);
  fast_VecSumErrBranch<LBN,R>(B, r);
}
