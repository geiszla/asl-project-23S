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
