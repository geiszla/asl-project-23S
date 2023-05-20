#include <math.h>

#ifdef _WIN32
#define alloca _alloca
#endif


#define dbl_prec 53
#define binSize 45

const double trennung = 0.000000000000000001;// 10^-18
// compile with g++ -std=c++11 ./main.cpp
// error free transforms 
/*
*/


///////////////////////////////// copied from CAMPARY package/////////////////////////////
double FPadd_rn(const double x, const double y){
  return x + y;
}
static inline double FPmul_rn(const double x, const double y){
  return x * y;
}
static inline double fma_d_rn_cpu(const double x, const double y, double xy) {
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

static inline double FPfma_rn(const double x, const double y, const double z){
	#ifdef FP_FAST_FMA
  //#warning cpu has fma
	  return fma(x, y, z);
	#else
   	return fma_d_rn_cpu(x, y, z);
	#endif
}
/* Computes fl(a*b) and err(a*b). */
static inline double two_prod(const double a, const double b, double &err){
	const double p = FPmul_rn(a, b);
 	err = FPfma_rn(a, b, -p);
 	return p;
}
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

// Sorry a function can't start with a 2 therefore I took the nearest solution of it
void twoSum(const double a,const  double b, double *s_res, double *e_res){
    double s = a+b;
    double t = s-b; 
    double e = (a-t) + (b-(s-t));
    // write ouput
    *s_res = s; *e_res = e;
    return;
}

void twoMultFMA(const double a,const double b,double *pi_res, double *e_res){
    double pi = a*b;
    double e = fma(a,b,-pi);
    *pi_res = pi; *e_res = e;
    return;

}


// call it with the array x and a array of the same size 
void vecSum(double *x, double *e_res, int in_out_size){
    int length = in_out_size;
    double* s = (double *)alloca(length*sizeof(double));
    s[length-1] = x[length-1];
    for(int i = length-2; i>=0; i--){
        double s_tmp,e_tmp; 
        twoSum(x[i],s[i+1],&s_tmp, &e_tmp);
        s[i] = s_tmp; e_res[i] = e_tmp;
    }
    e_res[0]=s[0];


    return;
}

/**
Implementation of VecSumErrBranch algorithm (Algorithm 7)

Input: e vector size n (S-nonoverlapping), output vector size m
Output: f vector size m

**/
void vecSumErrBranch(double* e, int n, int m, double *f){
   double* err = (double *)alloca(n*sizeof(double));

   int j = 0;
   err[0] = e[0];
   for (int i = 0; i <= n-2; i++) {
      twoSum(err[i], e[i+1], &f[j], &err[i+1]);
      if (err[i+1] != 0) {
         if (j >= m - 1){ // enough output terms
            return;
         }
         j++;
      } else {
         err[i+1] = f[j];
      }
   }
   if (err[n-1] != 0 && j < m) {
      f[j] = err[n-1];
   }

   
}

/**
Implementation of VecSumErr algorithm (Algorithm 8)

Input: f vector size n
Output: g vector size n

Only correct if n>2!!

Info: Probably faster to do it inplace in f then creating new g, but I tried to do it as similar as the algorithm.
**/
void vecSumErr(double* f, int n, double* g){
	int m = n-1;
	double* err = (double *)alloca(n*sizeof(double));

	err[0] = f[0];
	
	for (int i=0; i<= m-1; i++){
		twoSum(err[i],f[i+1],&g[i],&err[i+1]);
	}
	g[m] = err[m];

}



/** implementation of Algorithm 6 renormalization
 **/

void renormalizationalgorithm(double x[],int size_of_x , double f[], int m){
    double* err = (double *)alloca((size_of_x)*sizeof(double));
    double* f_tmp = (double *)alloca((m+1)*sizeof(double));
    for (int i = 0; i<=m; i++){
        f_tmp[i] = 0;
    }
    vecSum(x,err,size_of_x);

    vecSumErrBranch(err,size_of_x,m+1,f_tmp);
   
    for (int i = 0; i<=(m-2); i++){
        
        vecSumErr(&(f_tmp[i]), m-i+1, &(f_tmp[i]));
      
       f[i] = f_tmp[i];
    }
    f[m-1] = f_tmp[m-1];


}

// camapry merge only for testing
static inline void merge(double const *x, double const *y, double *z,int K,int L){
  int i=0, j=0;
  for(int n=0; n<K+L; n++){
    if(i==K || (j<L && fabs(y[j])>fabs(x[i]))){ z[n] = y[j]; j++; }
    else{ z[n] = x[i]; i++; }
  }
}


static inline void renorm_rand2L(int sX, int sR, double x[]){
	double pr;
	int i, j, ptr = 0;

  x[0] = two_sum(x[0], x[1], x[1]);
  for(i=2; i<sX; i++){
		pr = two_sum(x[i-1], x[i], x[i]);
    for(j=i-2; j>0; j--) pr = two_sum(x[j], pr, x[j+1]);
		x[0] = fast_two_sum(x[0], pr, x[1]);
	}

	i = 2;
  if(x[1] == 0.) pr = x[0];
  else { pr = x[1]; ptr++;}
  while(ptr<sR && i<sX){
    x[ptr] = fast_two_sum(pr, x[i], pr); i++;
    if(pr == 0.) pr = x[ptr]; else ptr++;
  }
  if(ptr<sR && pr!=0.){ x[ptr] = pr; ptr++; }
  for(i=ptr; i<sX; i++) x[i] = 0.;
}

/**Implementation of FP exansionaddition with k terms 
 * Input a and b of length k 
 * Output r of length k
*/
void addition(double *a, double *b, double *s, int length_a,int length_b, int length_result){
     double*  tmp = (double *)alloca((length_a+length_b)*sizeof(double));
     merge(a,b,tmp,length_a,length_b);
     renormalizationalgorithm(tmp,length_a+length_b,s,length_result);
    
    return;
}



/**Implementation of FP exansion multiplication with k terms 
 * Input a and b of length k 
 * Output r of length k
*/
void multiplication(double *a, double *b, double *r,int sizea, int sizeb, int sizer){
    int k = sizea;
    double* err = (double *)alloca((sizea*sizea -1)*sizeof(double));
  
    double* r_ext = (double *)alloca(sizea*sizeof(double));
    twoMultFMA(a[0],b[0],&(r_ext[0]), &(err[0]));
    for(int n=1; n<sizea; n++){
        double* e_tmp = (double *)alloca(n*sizeof(double));
        double* p = (double *)alloca(n*sizeof(double));
        for(int i = 0; i<=n; i++){
            twoMultFMA(a[i],b[n-i],&(p[i]), &(e_tmp[i]));
        }
        double* tmp = (double *)alloca((n*n+n)*sizeof(double));
        double* tmp2 = (double *)alloca((n*n+n+1)*sizeof(double));
        for(int b =0; b<=n; b++){
            tmp2[b]=p[b];
        }
        for(int b =n+1; b<n*n+n; b++){
            tmp2[b]=err[b-n];
        }
        // write result into tmp
        vecSum(tmp2, tmp,1);
        // r_N = tmp[0]
        r_ext[n] = 0;
        // now write e 0:n^2 -1
        for(int b =0; b< n*n-1; b++){
            err[b]= tmp[b+1];
        }

        // now compute e[0:(n+1)^2 -1] <- e[0:n^2 + n -1],e[0:n]
        int count = 0;
        for(int b =0; b<=((n^2) +n-1); b++){
            err[count] = err[b]; count++;
        }
        for(int b =0; b<=n; b++){
            err[count] = err[b]; count++;
        }
        // here compute line 9-11
        for ( int i =1; i<k; i++){
            r_ext[k] = r_ext[k] + a[i]*b[k-i];
        }
        // compute line 12 -14
         for ( int i =0; i<k*k; i++){
            r_ext[k] = r_ext[k] + err[i];
        }

    }

    for(int i = 1; i<=sizea*sizea-1; i++){
        r_ext[k]= r_ext[k]+ err[i];
    }
    renormalizationalgorithm(r_ext,k+1,r,sizea);
    // return
}
// helper
int exponent(double d)
{
  int result;
  frexp(d,&result);
  return result;
}

void deposit(double *p, double *b1, double *b2) {
	twoSum(*b1, *p, b1, p);
	*b2 = *b2 + *p;
}
/**
Implementation of Accumulate algorithm (Algorithm 2, Paper 2)

Input: 	double p, e
		vector<double> b
		int sh, l
Output: vector<double> b

**/

void accumulate(double p, double e, double *b, int sh, int l) {
	int c = dbl_prec - binSize - 1;
	if (l < binSize - 2*c - 1){
		deposit(&p, &b[sh], &b[sh+1]);
		deposit(&e, &b[sh+1], &b[sh+2]);
	} else if (l < binSize - c){
		deposit(&p, &b[sh], &b[sh+1]);
		twoSum(b[sh+1], e, &b[sh+1], &e);
		deposit(&e, &b[sh+2], &b[sh+3]);
	} else {
		twoSum(b[sh], p, &b[sh], &p);
		deposit(&p, &b[sh+1], &b[sh+2]);
		deposit(&e, &b[sh+2], &b[sh+3]);
	}
}

/**
Implementation of Renormalize algorithm (Algorithm 3, Paper 2)

Input: 	vector<double> x (size n)
		int n, k
Output: vector<double> r (ulp-nonoverlapping) (size k)
**/

void renormalize(double *x,double* r, int n, int k) {
	double eps = x[0];
	
	int j = 0;
	int i = 1;
	while (i < n && j < k) {
		twoSum(eps, x[i], &r[j], &eps);
		if (eps == 0) { // no overflow
			eps = r[j];
		} else {
			j++;
		}
		i++;
	}
	if (eps != 0 && j < k) {
		r[j] = eps;
	}
}


/**
Implementation of Multiplication algorithm (Algorithm 1, Paper 2)

Input: x vector size n (ulp-nonoverlapping)
		y vector size m (ulp-nonoverlapping)
Output: pi vector size r (ulp-nonoverlapping) = x*y

Constraint: n >=m
**/

void mult2(double* x, double* y, double*pi, int n, int m, int r){
    
	double* B = (double *)alloca((r*dbl_prec/binSize+2)*sizeof(double));
	// get sum of first exponents
	int e = exponent(x[0]) + exponent(y[0]);
	
	// initialize each Bin with starting value
	for (int i=0; i<r*dbl_prec/binSize+2;i++){
		B[i] = ldexp(1.5,e-(i+1)*binSize+dbl_prec-1); // 1.5*2^(e-(i+1)b+p-1)
	}
	int j,l,sh;
	double p, err;
	for(int i=0; i<= fmin(n-1, r); i++){
		for(j=0;j<= fmin(m-1,r-1-i); j++){
			//p = two_prod(x[i], y[j], &err); this leads to more similar result as CAMPARY
			twoMultFMA(x[i], y[j],&p,&err);
			l = e- exponent(x[i]) - exponent(y[j]); 
			sh = floor(l/binSize); // bin of the first pair
			l = l- sh*binSize; // number of leading bits
			accumulate(p,err,B,sh,l); // add to correct bins
		}
		j-=1;
		if (j < m-1){ //  I don't get what this part does
			p = x[i]*y[j];
			l = e- exponent(x[i]) - exponent(y[j]); 
			sh = floor(l/binSize); 
			l = l- sh*binSize; 
			accumulate(p,0.,B,sh,l);
		}
	}
	for (int i=0;  i<r*dbl_prec/binSize+2;i++){
		B[i] = B[i]-ldexp(1.5,e-(i+1)*binSize+dbl_prec-1); // B_i - 1.5*2^(e-(i+1)b+p-1)
	}

	
	renormalize(B,pi,r*dbl_prec/binSize+2,r);
}


