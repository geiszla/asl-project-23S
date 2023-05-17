#include <math.h>
#define dbl_prec 53
#define binSize 45
// the following functions are defined in basefunctions.cpp
void twoMultFMA(const double a,const double b,double *pi_res, double *e_res);
int exponent(double d);
void accumulate(double p, double e, double *b, int sh, int l);
void renormalize(double *x,double* r, int n, int k);



/*
Optimization 0: Precompute 
*/
void mult2_0(double* x, double* y,double*pi, int n, int m, int r){
    int bins = r*dbl_prec/binSize+2;
		double *B = (double *)alloca(bins*sizeof(double));
    // get exponents
    int *exp_x = (int *)alloca(n*sizeof(int));
	int *exp_y =(int *)alloca(m*sizeof(int));
    for(int i=0;i<n;i++){
        exp_x[i] = exponent(x[i]);
    }
        for(int j=0;j<m;j++){
        exp_y[j] = exponent(y[j]);
    }
	// get sum of first exponents
	int e = exp_x[0] + exp_x[0];
	

    double *B_start = (double *)alloca(bins*sizeof(double));
	// initialize each Bin with starting value
    B_start[0] = ldexp(1.5,e+dbl_prec-1-binSize); // since 1.5*2^(e-(i+1)b+p-1) == 1.5*2^(e-b+p-1)*2^(-b*i)
    B[0] = B_start[0];
    double const C1 = ldexp(1.0, -binSize); // multiply in each iteration last value by 2^-b instead of 2^(-b*i)
	for (int i=1; i<bins;i++){
        B_start[i] = B_start[i-1]*C1; // 1.5*2^(e-(i+1)b+p-1)
		B[i] = B_start[i];
	}
	int j,l,sh;
	double p, err;
    int max_terms = fmin(n-1, r);
	for(int i=0; i<= max_terms; i++){
		for(j=0;j<= fmin(m-1,r-1-i); j++){
			//p = two_prod(x[i], y[j], &err); this leads to more similar result as CAMPARY
			twoMultFMA(x[i], y[j],&p,&err);
			l = e- exp_x[i] - exp_y[j]; 
			sh = floor(l/binSize); // bin of the first pair
			l = l- sh*binSize; // number of leading bits
			accumulate(p,err,B,sh,l); // add to correct bins
		}
		j-=1;
		if (j < m-1){ 
			p = x[i]*y[j];
			l = e- exp_x[i] - exp_y[j]; 
			sh = floor(l/binSize); 
			l = l- sh*binSize; 
			accumulate(p,0.,B,sh,l);
		}
	}
	for (int i=0;  i<bins;i++){
		B[i] = B[i]-B_start[i]; // B_i - 1.5*2^(e-(i+1)b+p-1)
	}

	
	renormalize(B,pi,bins,r);
}
/*
Optimization 0: Precompute 
Optimization 1: replace complex with simpler functions: fmin(a,b) -> a<b?a:b, floor -> (int),a/b -> a*b^-1

*/
void mult2_1(double* x, double* y,double*pi, int n, int m, int r){
    int bins = r*dbl_prec/binSize+2;
    int b_inv = 1/binSize;
	double *B = (double *)alloca(bins*sizeof(double));
    // get exponents
    int *exp_x = (int *)alloca(n*sizeof(int));
		int *exp_y = (int *)alloca(m*sizeof(int));
    for(int i=0;i<n;i++){
        exp_x[i] = exponent(x[i]);
    }
        for(int j=0;j<m;j++){
        exp_y[j] = exponent(y[j]);
    }
	// get sum of first exponents
	int e = exp_x[0] + exp_x[0];
	

    double *B_start = (double *)alloca(bins*sizeof(double));
	// initialize each Bin with starting value
    B_start[0] = ldexp(1.5,e+dbl_prec-1-binSize); // since 1.5*2^(e-(i+1)b+p-1) == 1.5*2^(e-b+p-1)*2^(-b*i)
    B[0] = B_start[0];
    double const C1 = ldexp(1.0, -binSize); // multiply in each iteration last value by 2^-b instead of 2^(-b*i)
	for (int i=1; i<bins;i++){
        B_start[i] = B_start[i-1]*C1; // 1.5*2^(e-(i+1)b+p-1)
		B[i] = B_start[i];
	}
	int j,l,sh;
	double p, err;
    int max_terms_x = (n-1)<r?(n-1):r;
	for(int i=0; i<= max_terms_x; i++){ // 
        int max_terms_y = ((m-1)<(r-1-i)?(m-1):(r-1-i));
		for(j=0;j<= max_terms_y; j++){ //(m-1)<(r-1-i)?(m-1):(r-1-i)
			//p = two_prod(x[i], y[j], &err); this leads to more similar result as CAMPARY
			twoMultFMA(x[i], y[j],&p,&err);
			l = e- exp_x[i] - exp_y[j]; 
			sh = (int)(l*b_inv); // bin of the first pair
			l = l- sh*binSize; // number of leading bits
			accumulate(p,err,B,sh,l); // add to correct bins
		}
		j-=1;
		if (j < m-1){ 
			p = x[i]*y[j];
			l = e- exp_x[i] - exp_y[j]; 
			sh = (int)(l*b_inv); 
			l = l- sh*binSize; 
			accumulate(p,0.,B,sh,l);
		}
	}
	for (int i=0;  i<bins;i++){
		B[i] = B[i]-B_start[i]; // B_i - 1.5*2^(e-(i+1)b+p-1)
	}

	
	renormalize(B,pi,bins,r);
}

/*
Optimization 0: Precompute 
Optimization 1: replace complex with simpler functions: fmin(a,b) -> a<b?a:b, floor -> (int),a/b -> a*b^-1
Optimization 2: remove function calls accumulate,deposit, exponent, twoMultFma, renormalize

*/
void mult2_2(double* x, double* y,double*pi, int n, int m, int r){
    int bins = r*dbl_prec/binSize+2;
    int b_inv = 1/binSize;
	double *B = (double *)alloca(bins*sizeof(double));
    // get exponents
    int *exp_x = (int *)alloca(n*sizeof(int));
	int *exp_y = (int *)alloca(m*sizeof(int));
    for(int i=0;i<n;i++){
		frexp(x[i],&exp_x[i]);
    }
    for(int j=0;j<m;j++){
		frexp(y[j],&exp_y[j]);
    }
	// get sum of first exponents
	int e = exp_x[0] + exp_x[0];
	

    double *B_start = (double *)alloca(bins*sizeof(double));
	// initialize each Bin with starting value
    B_start[0] = ldexp(1.5,e+dbl_prec-1-binSize); // since 1.5*2^(e-(i+1)b+p-1) == 1.5*2^(e-b+p-1)*2^(-b*i)
    B[0] = B_start[0];
    double const C1 = ldexp(1.0, -binSize); // multiply in each iteration last value by 2^-b instead of 2^(-b*i)
	for (int i=1; i<bins;i++){
        B_start[i] = B_start[i-1]*C1; // 1.5*2^(e-(i+1)b+p-1)
		B[i] = B_start[i];
	}
	int j,l,sh;
	double p, err;
    int max_terms_x = (n-1)<r?(n-1):r;
	for(int i=0; i<= max_terms_x; i++){ // 
        int max_terms_y = ((m-1)<(r-1-i)?(m-1):(r-1-i));
		for(j=0;j<= max_terms_y; j++){ //(m-1)<(r-1-i)?(m-1):(r-1-i)
			
			//twoMultFMA(x[i], y[j],&p,&err);
			double temp_p = x[i]*y[j];
			err = fma(x[i],y[j],-temp_p);
			p = temp_p;


			l = e- exp_x[i] - exp_y[j]; 
			sh = (int)(l*b_inv); // bin of the first pair
			l = l- sh*binSize; // number of leading bits

			//accumulate
			int c = dbl_prec - binSize - 1;
			double s,t;
			if (l < binSize - 2*c - 1){
				//deposit(&p, &B[sh], &B[sh+1]);
				//twoSum(B[sh], p, &B[sh], &p);
				s = B[sh]+p;
				t = s-p; 
				B[sh] = s; 
				p = (B[sh]-t) + (p-(s-t));


				B[sh+1] = B[sh+1] + p;

				//deposit(&err, &B[sh+1], &B[sh+2]);
				//twoSum(B[sh+1], err,&B[sh+1], &err);
				s = B[sh+1]+err;
				t = s-err;
				B[sh+1] = s;
				err = (B[sh+1]-t) + (err-(s-t));


				B[sh+2] = B[sh+2] + err;
			} else if (l < binSize - c){
				//deposit(&p, &B[sh], &B[sh+1]);
				//twoSum(B[sh], p, &B[sh], &p);
				s = B[sh]+p;
				t = s-p; 
				B[sh] = s; 
				p = (B[sh]-t) + (p-(s-t));

				B[sh+1] = B[sh+1] + p;

				//twoSum(B[sh+1], err, &B[sh+1], &err);
				s = B[sh+1]+err;
				t = s-err;
				B[sh+1] = s;
				err = (B[sh+1]-t) + (err-(s-t));

				//deposit(&err, &B[sh+2], &B[sh+3]);
				//twoSum(B[sh+2], err, &B[sh+2], &err);
				s = B[sh+2]+err;
				t = s-err;
				B[sh+2] = s;
				err = (B[sh+2]-t) + (err-(s-t));

				B[sh+3] = B[sh+3]+ err;
			} else {
				//twoSum(B[sh], p, &B[sh], &p);
				s = B[sh]+p;
				t = s-p; 
				B[sh] = s; 
				p = (B[sh]-t) + (p-(s-t));
				//deposit(&p, &B[sh+1], &B[sh+2]);
				//twoSum(B[sh+1],p,&B[sh+1], &p);
				s = B[sh+1]+p;
				t = s-p; 
				B[sh+1] = s; 
				p = (B[sh+1]-t) + (p-(s-t));

				B[sh+2] = B[sh+2]+p;
				//deposit(&err, &B[sh+2], &B[sh+3]);
				//twoSum(B[sh+2],err, &B[sh+2], &err);
				s = B[sh+2]+err;
				t = s-err;
				B[sh+2] = s;
				err = (B[sh+2]-t) + (err-(s-t));

				B[sh+3] = B[sh+3]+ err;
			}
			
		}
		j-=1;
		if (j < m-1){ 
			p = x[i]*y[j];
			l = e- exp_x[i] - exp_y[j]; 
			sh = (int)(l*b_inv); 
			l = l- sh*binSize;
			//accumulate with ignoring errors
			int c = dbl_prec - binSize - 1;
			double s,t;
			if (l < binSize - 2*c - 1){
				//deposit(&p, &B[sh], &B[sh+1]);
				//twoSum(B[sh], p, &B[sh], &p);
				s = B[sh]+p;
				t = s-p; 
				B[sh] = s; 
				p = (B[sh]-t) + (p-(s-t));


				B[sh+1] = B[sh+1] + p;

			} else if (l < binSize - c){
				//deposit(&p, &B[sh], &B[sh+1]);
				//twoSum(B[sh], p, &B[sh], &p);
				s = B[sh]+p;
				t = s-p; 
				B[sh] = s; 
				p = (B[sh]-t) + (p-(s-t));

				B[sh+1] = B[sh+1] + p;

			} else {
				//twoSum(B[sh], p, &B[sh], &p);
				s = B[sh]+p;
				t = s-p; 
				B[sh] = s; 
				p = (B[sh]-t) + (p-(s-t));
				//deposit(&p, &B[sh+1], &B[sh+2]);
				//twoSum(B[sh+1],p,&B[sh+1], &p);
				s = B[sh+1]+p;
				t = s-p; 
				B[sh+1] = s; 
				p = (B[sh+1]-t) + (p-(s-t));

				B[sh+2] = B[sh+2]+p;

			}
		}
	}
	for (int i=0;  i<bins;i++){
		B[i] = B[i]-B_start[i]; // B_i - 1.5*2^(e-(i+1)b+p-1)
	}

	
	//renormalize(B,pi,bins,r);
 	double eps = B[0];
	
	j = 0;
	int i = 1;
	while (i < bins && j < r) {
		//twoSum(eps, B[i], &pi[j], &eps);
		double s = eps+B[i];
		double t = s-B[i]; 
		pi[j] = s; 
		eps= (eps-t) + (B[i]-(s-t));
		
		if (eps == 0) { // no overflow
			eps = pi[j];
		} else {
			j++;
		}
		i++;
	}
	if (eps != 0 && j < r) {
		pi[j] = eps;
	} 
}


/*
Optimization 0: Precompute 
Optimization 1: replace complex with simpler functions: fmin(a,b) -> a<b?a:b, floor -> (int),a/b -> a*b^-1
Optimization 2: remove function calls accumulate,deposit, exponent, twoMultFma, renormalize
Optimization 3: unroll most inner loop twice
*/
void mult2_3(double* x, double* y,double*pi, int n, int m, int r){
    int bins = r*dbl_prec/binSize+2;
    int b_inv = 1/binSize;
	int c = dbl_prec - binSize - 1;
	const int COND1 = binSize - 2*c - 1;
	const int COND2 = binSize - c;

	double *B = (double *)alloca(bins*sizeof(double));
    // get exponents
    int *exp_x = (int *)alloca(n*sizeof(int));
	int *exp_y = (int *)alloca(m*sizeof(int));
    for(int i=0;i<n;i++){
		frexp(x[i],&exp_x[i]);
    }
    for(int j=0;j<m;j++){
		frexp(y[j],&exp_y[j]);
    }
	// get sum of first exponents
	int e = exp_x[0] + exp_x[0];
	

    double *B_start = (double *)alloca(bins*sizeof(double));
	// initialize each Bin with starting value
    B_start[0] = ldexp(1.5,e+dbl_prec-1-binSize); // since 1.5*2^(e-(i+1)b+p-1) == 1.5*2^(e-b+p-1)*2^(-b*i)
    B[0] = B_start[0];
    double const C1 = ldexp(1.0, -binSize); // multiply in each iteration last value by 2^-b instead of 2^(-b*i)
	for (int i=1; i<bins;i++){
        B_start[i] = B_start[i-1]*C1; // 1.5*2^(e-(i+1)b+p-1)
		B[i] = B_start[i];
	}
	int j;
	int l_0, sh_0, l_1, sh_1;
	double p_0, err_0, p_1, err_1;
    int max_terms_x = (n-1)<r?(n-1):r;
	int i = 0;
	for(i=0; i<= max_terms_x; i++){ //
		double x_0 = x[i];
		double exp_x_0 = exp_x[i];
        int max_terms_y = ((m-1)<(r-1-i)?(m-1):(r-1-i));
		for(j=0;j<= max_terms_y; j+=2){ 
			double y_0 = y[j];
			double y_1 = y[j+1];
			double exp_y_0 = exp_y[j];
			double exp_y_1 = exp_y[j+1];

			// unroll i + 0, j + 0

			double temp_p_0 = x_0*y_0;
			err_0 = fma(x_0,y_0,-temp_p_0);
			p_0 = temp_p_0;


			l_0 = e- exp_x_0 - exp_y_0; 
			sh_0 = (int)(l_0*b_inv); 
			l_0 = l_0- sh_0*binSize; 

			// read from memory
			double B0_0 = B[sh_0 + 0];
			double B1_0 = B[sh_0 + 1];
			double B2_0 = B[sh_0 + 2];
			double B3_0 = B[sh_0 + 3];

			//accumulate
			
			double s0_0,t0_0, s1_0,t1_0, s2_0,t2_0; // temporary vairables, after the _ decides the loop unrolling
			if (l_0 < COND1){

				s0_0 = B0_0+p_0;
				t0_0= s0_0-p_0; 
				B0_0 = s0_0; 
				p_0 = (B0_0-t0_0) + (p_0-(s0_0-t0_0));


				B1_0 = B1_0 + p_0;


				s1_0 = B1_0+err_0;
				t1_0 = s1_0-err_0;
				B1_0 = s1_0;
				err_0 = (B1_0-t1_0) + (err_0-(s1_0-t1_0));


				B2_0 = B2_0 + err_0;
			} else if (l_0 < COND2){

				s0_0 = B0_0+p_0;
				t0_0 = s0_0-p_0; 
				B0_0 = s0_0; 
				p_0 = (B0_0-t0_0) + (p_0-(s0_0-t0_0));

				B1_0 = B1_0 + p_0;


				s1_0 = B1_0+err_0;
				t1_0 = s1_0-err_0;
				B1_0 = s1_0;
				err_0 = (B1_0-t1_0) + (err_0-(s1_0-t1_0));


				s2_0 = B2_0+err_0;
				t2_0 = s2_0-err_0;
				B2_0 = s2_0;
				err_0 = (B2_0-t2_0) + (err_0-(s2_0-t2_0));

				B3_0 = B3_0+ err_0;
			} else {

				s0_0 = B0_0+p_0;
				t0_0 = s0_0-p_0; 
				B0_0 = s0_0; 
				p_0 = (B0_0-t0_0) + (p_0-(s0_0-t0_0));

				s1_0 = B1_0+p_0;
				t1_0 = s1_0-p_0; 
				B1_0 = s1_0; 
				p_0 = (B1_0-t1_0) + (p_0-(s1_0-t1_0));

				B2_0 = B2_0+p_0;

				s2_0 = B2_0+err_0;
				t2_0 = s2_0-err_0;
				B2_0 = s2_0;
				err_0 = (B2_0-t2_0) + (err_0-(s2_0-t2_0));

				B3_0 = B3_0+ err_0;
			}
			// write back to memory
			B[sh_0 + 0] = B0_0;
			B[sh_0 + 1] = B1_0;
			B[sh_0 + 2] = B2_0;
			B[sh_0 + 3] = B3_0;

			// unroll i + 0, j + 1

			double temp_p_1 = x_0*y_1;
			err_1 = fma(x_0,y_1,-temp_p_1);
			p_1 = temp_p_1;


			l_1 = e- exp_x_0 - exp_y_1; 
			sh_1 = (int)(l_1*b_inv);
			l_1 = l_1- sh_1*binSize; 

			if (abs(sh_0 -sh_1) > 4){
				std::cout << "unrolling!"<< sh_0<<" "<< sh_1 << std::endl;
			}

			// read from memory --> HERE is a RAW dependancy therefore loop unrolling is useless
			double B0_1 = B[sh_1 + 0];
			double B1_1 = B[sh_1 + 1];
			double B2_1 = B[sh_1 + 2];
			double B3_1 = B[sh_1 + 3];
			
			//accumulate 
			double s0_1,t0_1, s1_1,t1_1, s2_1,t2_1; // temporary vairables, after the _ decides the loop unrolling
			if (l_1 < COND1){

				s0_1 = B0_1+p_1;
				t0_1= s0_1-p_1; 
				B0_1 = s0_1; 
				p_1 = (B0_1-t0_1) + (p_1-(s0_1-t0_1));


				B1_1 = B1_1 + p_1;


				s1_1 = B1_1+err_1;
				t1_1 = s1_1-err_1;
				B1_1 = s1_1;
				err_1 = (B1_1-t1_1) + (err_1-(s1_1-t1_1));


				B2_1 = B2_1 + err_1;
			} else if (l_1 < COND2){

				s0_1 = B0_1+p_1;
				t0_1 = s0_1-p_1; 
				B0_1 = s0_1; 
				p_1 = (B0_1-t0_1) + (p_1-(s0_1-t0_1));

				B1_1 = B1_1 + p_1;

				
				s1_1 = B1_1+err_1;
				t1_1 = s1_0-err_1;
				B1_1 = s1_1;
				err_1 = (B1_1-t1_1) + (err_1-(s1_1-t1_1));


				s2_1 = B2_1+err_1;
				t2_1 = s2_1-err_1;
				B2_1 = s2_1;
				err_1 = (B2_1-t2_1) + (err_1-(s2_1-t2_1));

				B3_1 = B3_1+ err_1;
			} else {
				s0_1 = B0_1+p_1;
				t0_1 = s0_1-p_1; 
				B0_1 = s0_1; 
				p_1 = (B0_1-t0_1) + (p_1-(s0_1-t0_1));

				s1_1 = B1_1+p_1;
				t1_1 = s1_1-p_1; 
				B1_1 = s1_1; 
				p_1 = (B1_1-t1_1) + (p_1-(s1_1-t1_1));

				B2_1 = B2_1+p_1;

				s2_1 = B2_1+err_1;
				t2_1 = s2_1-err_1;
				B2_1 = s2_1;
				err_1 = (B2_1-t2_1) + (err_1-(s2_1-t2_1));

				B3_1 = B3_1+ err_1;
			}
			// write back to memory
			B[sh_1 + 0] = B0_1;
			B[sh_1 + 1] = B1_1;
			B[sh_1 + 2] = B2_1;
			B[sh_1 + 3] = B3_1;			

		}
		// finish remaining elements in j
		for(;j<= max_terms_y; j++){ 
			double y_0 = y[j];
			double exp_y_0 = exp_y[j];

			double temp_p_0 = x_0*y_0;
			err_0 = fma(x_0,y_0,-temp_p_0);
			p_0 = temp_p_0;


			l_0 = e- exp_x_0 - exp_y_0; 
			sh_0 = (int)(l_0*b_inv); 
			l_0 = l_0- sh_0*binSize;

			// read from memory
			double B0_0 = B[sh_0 + 0];
			double B1_0 = B[sh_0 + 1];
			double B2_0 = B[sh_0 + 2];
			double B3_0 = B[sh_0 + 3];
			//accumulate

			double s0_0,t0_0, s1_0,t1_0, s2_0,t2_0; // temporary vairables, after the _ decides the loop unrolling
			if (l_0 < COND1){

				s0_0 = B0_0+p_0;
				t0_0= s0_0-p_0; 
				B0_0 = s0_0; 
				p_0 = (B0_0-t0_0) + (p_0-(s0_0-t0_0));


				B1_0 = B1_0 + p_0;


				s1_0 = B1_0+err_0;
				t1_0 = s1_0-err_0;
				B1_0 = s1_0;
				err_0 = (B1_0-t1_0) + (err_0-(s1_0-t1_0));


				B2_0 = B2_0 + err_0;

			} else if (l_0 < COND2){

				s0_0 = B0_0+p_0;
				t0_0 = s0_0-p_0; 
				B0_0 = s0_0; 
				p_0 = (B0_0-t0_0) + (p_0-(s0_0-t0_0));

				B1_0 = B1_0 + p_0;

				s1_0 = B1_0+err_0;
				t1_0 = s1_0-err_0;
				B1_0 = s1_0;
				err_0 = (B1_0-t1_0) + (err_0-(s1_0-t1_0));

				s2_0 = B2_0+err_0;
				t2_0 = s2_0-err_0;
				B2_0 = s2_0;
				err_0 = (B2_0-t2_0) + (err_0-(s2_0-t2_0));

				B3_0 = B3_0+ err_0;
			} else {

				s0_0 = B0_0+p_0;
				t0_0 = s0_0-p_0; 
				B0_0 = s0_0; 
				p_0 = (B0_0-t0_0) + (p_0-(s0_0-t0_0));

				s1_0 = B1_0+p_0;
				t1_0 = s1_0-p_0; 
				B1_0 = s1_0; 
				p_0 = (B1_0-t1_0) + (p_0-(s1_0-t1_0));

				B2_0 = B2_0+p_0;

				s2_0 = B2_0+err_0;
				t2_0 = s2_0-err_0;
				B2_0 = s2_0;
				err_0 = (B2_0-t2_0) + (err_0-(s2_0-t2_0));

				B3_0 = B3_0+ err_0;
			}
			// write back to memory
			B[sh_0 + 0] = B0_0;
			B[sh_0 + 1] = B1_0;
			B[sh_0 + 2] = B2_0;
			B[sh_0 + 3] = B3_0;			

		}
		j-=1;
		if (j < m-1){
			double p;
			int l, sh;

			p = x_0*y[j];
			l = e- exp_x_0 - exp_y[j]; 
			sh = (int)(l*b_inv); 
			l = l- sh*binSize;

			// read from memory
			double B0 = B[sh+0];
			double B1 = B[sh+1];
			double B2 = B[sh+2];

			//accumulate with ignoring errors

			double s,t;
			if (l < binSize - 2*c - 1){
	
				s = B0+p;
				t = s-p; 
				B0 = s; 
				p = (B0-t) + (p-(s-t));


				B[sh+1] = B[sh+1] + p;

			} else if (l < binSize - c){

				s = B0+p;
				t = s-p; 
				B0 = s; 
				p = (B0-t) + (p-(s-t));

				B1 = B1 + p;

			} else {

				s = B0+p;
				t = s-p; 
				B0 = s; 
				p = (B0-t) + (p-(s-t));

				s = B1+p;
				t = s-p; 
				B1 = s; 
				p = (B1-t) + (p-(s-t));

				B2 = B2+p;

			}
			// write to memory
			B[sh+0] = B0;
			B[sh+1] = B1;
			B[sh+2] = B2;

			
		}
	}


	for (int i=0;  i<bins;i++){
		B[i] = B[i]-B_start[i]; // B_i - 1.5*2^(e-(i+1)b+p-1)
	}

	//renormalize
	double eps = B[0];
	
	j = 0;
	i = 1;
	while (i < bins && j < r) {
		//twoSum(eps, B[i], &pi[j], &eps);
		double s = eps+B[i];
		double t = s-B[i]; 
		pi[j] = s; 
		eps= (eps-t) + (B[i]-(s-t));
		
		if (eps == 0) { // no overflow
			eps = pi[j];
		} else {
			j++;
		}
		i++;
	}
	if (eps != 0 && j < r) {
		pi[j] = eps;
	} 
}