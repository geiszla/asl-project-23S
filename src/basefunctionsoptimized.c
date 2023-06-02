#include <math.h>
#include <immintrin.h>
const double trennung = 0.000000000000000001; // 10^-18
// compile with g++ -std=c++11 ./main.cpp
// error free transforms
/*
 */

// Sorry a function can't start with a 2 therefore I took the nearest solution of it
void twoSum(const double a, const double b, double *s_res, double *e_res)
{
    double s = a + b;
    double t = s - b;
    double e = (a - t) + (b - (s - t));
    // write ouput
    *s_res = s;
    *e_res = e;
    return;
}

void twoMultFMA(const double a, const double b, double *pi_res, double *e_res)
{
    double pi = a * b;
    double e = fma(a, b, -pi);
    *pi_res = pi;
    *e_res = e;
    return;
}

// call it with the array x and a array of the same size
void vecSum(double *x, double *e_res, int in_out_size)
{
    int length = in_out_size;
    double *s = (double *)alloca(length * sizeof(double));
    s[length - 1] = x[length - 1];
    for (int i = length - 2; i >= 0; i--)
    {
        double s_tmp, e_tmp;
        twoSum(x[i], s[i + 1], &s_tmp, &e_tmp);
        s[i] = s_tmp;
        e_res[i + 1] = e_tmp;
    }
    e_res[0] = s[0];

    return;
}

/**
Implementation of VecSumErrBranch algorithm (Algorithm 7)

Input: e vector size n (S-nonoverlapping), output vector size m
Output: f vector size m

**/
void vecSumErrBranch(double *e, int n, int m, double *f)
{
    double *err = (double *)alloca(n * sizeof(double));

    int j = 0;
    err[0] = e[0];
    for (int i = 0; i <= n - 2; i++)
    {
        twoSum(err[i], e[i + 1], &f[j], &err[i + 1]);
        if (err[i + 1] != 0)
        {
            if (j >= m - 1)
            { // enough output terms
                return;
            }
            j++;
        }
        else
        {
            err[i + 1] = f[j];
        }
    }
    if (err[n - 1] != 0 && j < m)
    {
        f[j] = err[n - 1];
    }
}

/**
Implementation of VecSumErr algorithm (Algorithm 8)

Input: f vector size n
Output: g vector size n

Only correct if n>2!!

Info: Probably faster to do it inplace in f then creating new g, but I tried to do it as similar as the algorithm.
**/
void vecSumErr(double *f, int n, double *g)
{
    int m = n - 1;
    double *err = (double *)alloca(n * sizeof(double));

    err[0] = f[0];

    for (int i = 0; i <= m - 1; i++)
    {
        twoSum(err[i], f[i + 1], &g[i], &err[i + 1]);
    }
    g[m] = err[m];
}

/** implementation of Algorithm 6 renormalization
 **/

void renormalizationalgorithm(double *x, int size_of_x, double *f, int m)
{

    double *err = (double *)alloca((size_of_x) * sizeof(double));
    double *f_tmp = (double *)alloca((m + 1) * sizeof(double));
    for (int i = 0; i <= m; i++)
    {
        f_tmp[i] = 0;
    }
    vecSum(x, err, size_of_x);

    vecSumErrBranch(err, size_of_x, m + 1, f_tmp);

    for (int i = 0; i <= (m - 2); i++)
    {

        vecSumErr(&(f_tmp[i]), m - i + 1, &(f_tmp[i]));

        f[i] = f_tmp[i];
    }
    f[m - 1] = f_tmp[m - 1];

    return;
}

/**Implementation of FP exansionaddition with k terms
 * Input a and b of length k
 * Output r of length k
 */
void addition(double *a, double *b, double *s, int length_a, int length_b, int length_result)
{
    double *tmp = (double *)alloca((length_a + length_b) * sizeof(double));
    merge(a, b, tmp, length_a, length_b);
    renormalizationalgorithm(tmp, length_a + length_b, s, length_result);

    return;
}

/**Implementation of FP exansion multiplication with k terms
 * Input a and b of length k
 * Output r of length k
 */
void multiplication(double *a, double *b, double *r, const int sizea, const int sizeb, const int sizer)
{

    int k = sizea;
    double *err = (double *)calloc((sizea * sizea + 2 * sizea+1), sizeof(double));
    double *r_ext = (double *)calloc((sizea * sizea), sizeof(double));

    twoMultFMA(a[0], b[0], &(r_ext[0]), &(err[0]));

    for (int n = 1; n <= (k - 1); n++)
    {
        double *e_tmp = (double *)calloc((n+1), sizeof(double));
        double *p = (double *)calloc((n + 1), sizeof(double));

        for (int i = 0; i <= n; i++)
        {
            twoMultFMA(a[i], b[n - i], &(p[i]), &(e_tmp[i]));
            //------------------------------------------------
            
            //-------------------------------- twoMultFMa inline
        }
        double *tmp = &err[-n - 1];
        double *tmp1 = (double *)calloc((n * n + n+1), sizeof(double));

        //-------------------------------- vecsum inline
        int length = (n * n + n);
        double *s = (double *)alloca(length * sizeof(double));
        s[length - 1] = err[n * n - 1];
        for (int i = length - 2; i >= 0; i--)
        {
            if (i <= (n))
            {
                tmp = p;
            }
            double s_tmp, e_tmplocal;
            
            double sl = tmp[i] + s[i + 1];
            double t = sl - s[i + 1];
            double e = (tmp[i] - t) + (s[i + 1] - (sl - t));
            s_tmp = sl;
            e_tmplocal = e;
            //--------------------------------
            s[i] = s_tmp;
            tmp1[i + 1] = e_tmplocal;
        }
        r_ext[n] = s[0];

        //--------------------------------
        /* write  tmp1 into r_ext[n], e[0:n^2 +n-1]*/
       
        for (int i = 0; i <= (n * n + n - 1); i++)
        {
            err[i] = tmp1[i + 1];
        }
        // writes err[0:n^2 +n-1],e_tmp[0:n] into err[0:n^2 +2n-1]
        for (int i = (n * n+n); i <= (n*n + 2*n -1); i++)
        {
            err[i] = e_tmp[i - (n * n+n) ];
           
        }
    }

    for (int i = 1; i <= (k - 1); i++)
    {

        r_ext[k] += a[i] * b[k - i];
    }

    for (int i = 0; i <= (k * k - 1); i++)
    {

        r_ext[i] += err[i];
    }

    renormalizationalgorithm(r_ext, k + 1, r, sizea);
    return;
}


void multiplicationsimdunrolled(double *a, double *b, double *r, const int sizea, const int sizeb, const int sizer)
{

    int k = sizea;
    double *err = (double *)calloc((sizea * sizea + 2 * sizea+1), sizeof(double));
    double *r_ext = (double *)calloc((sizea * sizea), sizeof(double));

    double pi = a[0] *b[0];
    double e = fma(a[0], b[0], -pi);
    r_ext[0] = pi;
    err[0] = e;

    for (int n = 1; n <= (k - 1); n++)
    {
        double *e_tmp = (double *)calloc((n+1), sizeof(double));
        double *p = (double *)calloc((n + 1), sizeof(double));
        int i =0;
        for (; i <= n-32; i+=32)
        {
            __m256d a1 = _mm256_loadu_pd(&a[i]);
            __m256d b1 = _mm256_loadu_pd(&b[n - i-3]);
            b1 =  _mm256_permute4x64_pd(b1, 0b00011011);
            __m256d pi1 = _mm256_mul_pd(a1, b1);
            __m256d e1 = _mm256_fmsub_pd(a1, b1, pi1);
            _mm256_storeu_pd(&p[i], pi1);
            _mm256_storeu_pd(&e_tmp[i], e1);

            __m256d a2 = _mm256_loadu_pd(&a[i+4]);
            __m256d b2 = _mm256_loadu_pd(&b[n - i-7]);
            b2 =  _mm256_permute4x64_pd(b2, 0b00011011);
            __m256d pi2 = _mm256_mul_pd(a2, b2);
            __m256d e2 = _mm256_fmsub_pd(a2, b2, pi2);
            _mm256_storeu_pd(&p[i+4], pi2);
            _mm256_storeu_pd(&e_tmp[i+4], e2);

            __m256d a3 = _mm256_loadu_pd(&a[i+8]);
            __m256d b3 = _mm256_loadu_pd(&b[n - i-11]);
            b3 =  _mm256_permute4x64_pd(b3, 0b00011011);
            __m256d pi3 = _mm256_mul_pd(a3, b3);
            __m256d e3 = _mm256_fmsub_pd(a3, b3, pi3);
            _mm256_storeu_pd(&p[i+8], pi3);
            _mm256_storeu_pd(&e_tmp[i+8], e3);

            __m256d a4 = _mm256_loadu_pd(&a[i+12]);
            __m256d b4 = _mm256_loadu_pd(&b[n - i-15]);
            b4 =  _mm256_permute4x64_pd(b4, 0b00011011);
            __m256d pi4 = _mm256_mul_pd(a4, b4);
            __m256d e4 = _mm256_fmsub_pd(a4, b4, pi4);
            _mm256_storeu_pd(&p[i+12], pi4);
            _mm256_storeu_pd(&e_tmp[i+12], e4);

            __m256d a5 = _mm256_loadu_pd(&a[i+16]);
            __m256d b5 = _mm256_loadu_pd(&b[n - i-19]);
            b5 =  _mm256_permute4x64_pd(b5, 0b00011011);
            __m256d pi5 = _mm256_mul_pd(a5, b5);
            __m256d e5 = _mm256_fmsub_pd(a5, b5, pi5);
            _mm256_storeu_pd(&p[i+16], pi5);
            _mm256_storeu_pd(&e_tmp[i+16], e5);

            __m256d a6 = _mm256_loadu_pd(&a[i+20]);
            __m256d b6 = _mm256_loadu_pd(&b[n - i-23]);
            b6 =  _mm256_permute4x64_pd(b6, 0b00011011);
            __m256d pi6 = _mm256_mul_pd(a6, b6);
            __m256d e6 = _mm256_fmsub_pd(a6, b6, pi6);
            _mm256_storeu_pd(&p[i+20], pi6);
            _mm256_storeu_pd(&e_tmp[i+20], e6);

            __m256d a7 = _mm256_loadu_pd(&a[i+24]);
            __m256d b7 = _mm256_loadu_pd(&b[n - i-27]);
            b7 =  _mm256_permute4x64_pd(b7, 0b00011011);
            __m256d pi7 = _mm256_mul_pd(a7, b7);
            __m256d e7 = _mm256_fmsub_pd(a7, b7, pi7);
            _mm256_storeu_pd(&p[i+24], pi7);
            _mm256_storeu_pd(&e_tmp[i+24], e7);

            __m256d a8 = _mm256_loadu_pd(&a[i+28]);
            __m256d b8 = _mm256_loadu_pd(&b[n - i-31]);
            b8 =  _mm256_permute4x64_pd(b8, 0b00011011);
            __m256d pi8 = _mm256_mul_pd(a8, b8);
            __m256d e8 = _mm256_fmsub_pd(a8, b8, pi8);
            _mm256_storeu_pd(&p[i+28], pi8);
            _mm256_storeu_pd(&e_tmp[i+28], e8);
            

           
        }
        for (; i <= n; i++)
        {
            //twoMultFMA(a[i], b[n - i], &(e_tmp[i]), &(e_tmp[i]));
            double pi = a[i] *b[n - i];
            double e = fma(a[i], b[n - i], -pi);
            p[i] = pi;
            e_tmp[i] = e;
        }
        double *tmp = &err[-n - 1];
        double *tmp1 = (double *)calloc((n * n + n+1), sizeof(double));

        //-------------------------------- vecsum inline
        int length = (n * n + n);
        double *s = (double *)alloca(length * sizeof(double));
        s[length - 1] = err[n * n - 1];
        for (int i = length - 2; i >= n+1; i-=1)
        {
            
            double s_tmp, e_tmplocal;
           

            //-------------------------------- two sum inline
            double sl = tmp[i] + s[i + 1];
            double t = sl - s[i + 1];
            double e = (tmp[i] - t) + (s[i + 1] - (sl - t));
            s_tmp = sl;
            e_tmplocal = e;
            //--------------------------------
            s[i] = s_tmp;
            tmp1[i + 1] = e_tmplocal;

            

        }
        tmp = p;
        for (int i = n; i >= 0; i--)
        {
            
            
            double s_tmp, e_tmplocal;
           

            //-------------------------------- two sum inline
            double sl = tmp[i] + s[i + 1];
            double t = sl - s[i + 1];
            double e = (tmp[i] - t) + (s[i + 1] - (sl - t));
            s_tmp = sl;
            e_tmplocal = e;
            //--------------------------------
            s[i] = s_tmp;
            tmp1[i + 1] = e_tmplocal;
        }
        r_ext[n] = s[0];

        //--------------------------------
        /* write  tmp1 into r_ext[n], e[0:n^2 +n-1]*/
       i = 0;
       for(; i <= (n * n + n - 1)-32; i+=32)
        {
            __m256d tmp1_vec = _mm256_loadu_pd(&tmp1[i + 1]);
            _mm256_storeu_pd(&err[i], tmp1_vec);
            __m256d tmp1_vec1 = _mm256_loadu_pd(&tmp1[i + 1+4]);
            _mm256_storeu_pd(&err[i+4], tmp1_vec1);
            __m256d tmp1_vec2 = _mm256_loadu_pd(&tmp1[i + 1+8]);
            _mm256_storeu_pd(&err[i+8], tmp1_vec2);
            __m256d tmp1_vec3 = _mm256_loadu_pd(&tmp1[i + 1+12]);
            _mm256_storeu_pd(&err[i+12], tmp1_vec3);
            __m256d tmp1_vec4 = _mm256_loadu_pd(&tmp1[i + 1+16]);
            _mm256_storeu_pd(&err[i+16], tmp1_vec4);
            __m256d tmp1_vec5 = _mm256_loadu_pd(&tmp1[i + 1+20]);
            _mm256_storeu_pd(&err[i+20], tmp1_vec5);
            __m256d tmp1_vec6 = _mm256_loadu_pd(&tmp1[i + 1+24]);
            _mm256_storeu_pd(&err[i+24], tmp1_vec6);
            __m256d tmp1_vec7 = _mm256_loadu_pd(&tmp1[i + 1+28]);
            _mm256_storeu_pd(&err[i+28], tmp1_vec7);
          
           
            //err[i] = tmp1[i + 1];
        }
        for (; i <= (n * n + n - 1); i++)
        {
            err[i] = tmp1[i + 1];
        }
        //---------------------------------------------

        // writes err[0:n^2 +n-1],e_tmp[0:n] into err[0:n^2 +2n-1]
        i=(n * n+n); 
        for (; i <= (n*n + 2*n -1)-32; i+=32)
        {
            __m256d tmp1_vec = _mm256_loadu_pd(&tmp1[i - (n * n+n) ]);
            _mm256_storeu_pd(&err[i], tmp1_vec);
            __m256d tmp1_vec1 = _mm256_loadu_pd(&tmp1[i +4- (n * n+n) ]);
            _mm256_storeu_pd(&err[i+4], tmp1_vec1);
            __m256d tmp1_vec2 = _mm256_loadu_pd(&tmp1[i +8- (n * n+n) ]);
            _mm256_storeu_pd(&err[i+8], tmp1_vec2);
            __m256d tmp1_vec3 = _mm256_loadu_pd(&tmp1[i +12- (n * n+n) ]);
            _mm256_storeu_pd(&err[i+12], tmp1_vec3);
            __m256d tmp1_vec4 = _mm256_loadu_pd(&tmp1[i +16- (n * n+n) ]);
            _mm256_storeu_pd(&err[i+16], tmp1_vec4);
            __m256d tmp1_vec5 = _mm256_loadu_pd(&tmp1[i +20- (n * n+n) ]);
            _mm256_storeu_pd(&err[i+20], tmp1_vec5);
            __m256d tmp1_vec6 = _mm256_loadu_pd(&tmp1[i +24- (n * n+n) ]);
            _mm256_storeu_pd(&err[i+24], tmp1_vec6);
            __m256d tmp1_vec7 = _mm256_loadu_pd(&tmp1[i +28- (n * n+n) ]);
            _mm256_storeu_pd(&err[i+28], tmp1_vec7);
          
           
        }
        for (; i <= (n*n + 2*n -1); i++)
        {
            err[i] = e_tmp[i - (n * n+n) ];
           
        }
    }
    int i = 1;
    for (; i <= (k - 1)-32; i+=32)
    {
        __m256d a = _mm256_loadu_pd(&a[i]);
        __m256d b = _mm256_loadu_pd(&b[k - i-3]);
        __m256d r_ext = _mm256_loadu_pd(&r_ext[k]);
        __m256d tmp1_vec = _mm256_permute4x64_pd(b, 0b00011011);
        __m256d  mul = _mm256_mul_pd(a, tmp1_vec);
        r_ext = _mm256_add_pd(r_ext, mul);
      
       
       __m256d a1 = _mm256_loadu_pd(&a[i+4]);
        __m256d b1 = _mm256_loadu_pd(&b[k - i-3-4]);
        __m256d tmp1_vec1 = _mm256_permute4x64_pd(b1, 0b00011011);
        __m256d  mul1 = _mm256_mul_pd(a1, tmp1_vec1);
        r_ext = _mm256_add_pd(r_ext, mul1);

        __m256d a2 = _mm256_loadu_pd(&a[i+8]);
        __m256d b2 = _mm256_loadu_pd(&b[k - i-3-8]);
        __m256d tmp1_vec2 = _mm256_permute4x64_pd(b2, 0b00011011);
        __m256d  mul2 = _mm256_mul_pd(a2, tmp1_vec2);
        r_ext = _mm256_add_pd(r_ext, mul2);

        __m256d a3 = _mm256_loadu_pd(&a[i+12]);
        __m256d b3 = _mm256_loadu_pd(&b[k - i-3-12]);
        __m256d tmp1_vec3 = _mm256_permute4x64_pd(b3, 0b00011011);
        __m256d  mul3 = _mm256_mul_pd(a3, tmp1_vec3);
        r_ext = _mm256_add_pd(r_ext, mul3);

        __m256d a4 = _mm256_loadu_pd(&a[i+16]);
        __m256d b4 = _mm256_loadu_pd(&b[k - i-3-16]);
        __m256d tmp1_vec4 = _mm256_permute4x64_pd(b4, 0b00011011);
        __m256d  mul4 = _mm256_mul_pd(a4, tmp1_vec4);
        r_ext = _mm256_add_pd(r_ext, mul4);

        __m256d a5 = _mm256_loadu_pd(&a[i+20]);
        __m256d b5 = _mm256_loadu_pd(&b[k - i-3-20]);
        __m256d tmp1_vec5 = _mm256_permute4x64_pd(b5, 0b00011011);
        __m256d  mul5 = _mm256_mul_pd(a5, tmp1_vec5);
        r_ext = _mm256_add_pd(r_ext, mul5);

        __m256d a6 = _mm256_loadu_pd(&a[i+24]);
        __m256d b6 = _mm256_loadu_pd(&b[k - i-3-24]);
        __m256d tmp1_vec6 = _mm256_permute4x64_pd(b6, 0b00011011);
        __m256d  mul6 = _mm256_mul_pd(a6, tmp1_vec6);
        r_ext = _mm256_add_pd(r_ext, mul6);

        __m256d a7 = _mm256_loadu_pd(&a[i+28]);
        __m256d b7 = _mm256_loadu_pd(&b[k - i-3-28]);
        __m256d tmp1_vec7 = _mm256_permute4x64_pd(b7, 0b00011011);
        __m256d  mul7 = _mm256_mul_pd(a7, tmp1_vec7);
        r_ext = _mm256_add_pd(r_ext, mul7);
        
       _mm256_storeu_pd(&r_ext[k],r_ext);
        
        

    }
    for (; i <= (k - 1); i++)
    {

        r_ext[k] += a[i] * b[k - i];
    }

    i=0;
    for (; i <= (k * k - 1)-32; i+=32)
    {
        __m256d tmp1_vec = _mm256_loadu_pd(&err[i]);
        __m256d tmp2_vec = _mm256_loadu_pd(&r_ext[i]);
        _mm256_storeu_pd(&r_ext[i], _mm256_add_pd(tmp1_vec, tmp2_vec));
        __m256d tmp1_vec1 = _mm256_loadu_pd(&err[i+4]);
        __m256d tmp2_vec1 = _mm256_loadu_pd(&r_ext[i+4]);
        _mm256_storeu_pd(&r_ext[i+4], _mm256_add_pd(tmp1_vec1, tmp2_vec1));
        __m256d tmp1_vec2 = _mm256_loadu_pd(&err[i+8]);
        __m256d tmp2_vec2 = _mm256_loadu_pd(&r_ext[i+8]);
        _mm256_storeu_pd(&r_ext[i+8], _mm256_add_pd(tmp1_vec2, tmp2_vec2));
        __m256d tmp1_vec3 = _mm256_loadu_pd(&err[i+12]);
        __m256d tmp2_vec3 = _mm256_loadu_pd(&r_ext[i+12]);
        _mm256_storeu_pd(&r_ext[i+12], _mm256_add_pd(tmp1_vec3, tmp2_vec3));
        __m256d tmp1_vec4 = _mm256_loadu_pd(&err[i+16]);
        __m256d tmp2_vec4 = _mm256_loadu_pd(&r_ext[i+16]);
        _mm256_storeu_pd(&r_ext[i+16], _mm256_add_pd(tmp1_vec4, tmp2_vec4));
        __m256d tmp1_vec5 = _mm256_loadu_pd(&err[i+20]);
        __m256d tmp2_vec5 = _mm256_loadu_pd(&r_ext[i+20]);
        _mm256_storeu_pd(&r_ext[i+20], _mm256_add_pd(tmp1_vec5, tmp2_vec5));
        __m256d tmp1_vec6 = _mm256_loadu_pd(&err[i+24]);
        __m256d tmp2_vec6 = _mm256_loadu_pd(&r_ext[i+24]);
        _mm256_storeu_pd(&r_ext[i+24], _mm256_add_pd(tmp1_vec6, tmp2_vec6));
        __m256d tmp1_vec7 = _mm256_loadu_pd(&err[i+28]);
        __m256d tmp2_vec7 = _mm256_loadu_pd(&r_ext[i+28]);
        _mm256_storeu_pd(&r_ext[i+28], _mm256_add_pd(tmp1_vec7, tmp2_vec7));

    }
     for (; i <= (k * k - 1); i++)
    {

        r_ext[i] += err[i];
    }

    renormalizationalgorithm(r_ext, k + 1, r, sizea);
    return;
}

// helper
int exponent(double d)
{
    int result;
    frexp(d, &result);
    return result;
}

void deposit(double *p, double *b1, double *b2)
{
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

void accumulate(double p, double e, double *b, int sh, int l)
{
    int c = dbl_prec - binSize - 1;
    if (l < binSize - 2 * c - 1)
    {
        deposit(&p, &b[sh], &b[sh + 1]);
        deposit(&e, &b[sh + 1], &b[sh + 2]);
    }
    else if (l < binSize - c)
    {
        deposit(&p, &b[sh], &b[sh + 1]);
        twoSum(b[sh + 1], e, &b[sh + 1], &e);
        deposit(&e, &b[sh + 2], &b[sh + 3]);
    }
    else
    {
        twoSum(b[sh], p, &b[sh], &p);
        deposit(&p, &b[sh + 1], &b[sh + 2]);
        deposit(&e, &b[sh + 2], &b[sh + 3]);
    }
}

/**
Implementation of Renormalize algorithm (Algorithm 3, Paper 2)

Input: 	vector<double> x (size n)
        int n, k
Output: vector<double> r (ulp-nonoverlapping) (size k)
**/

void renormalize(double *x, double *r, int n, int k)
{
    double eps = x[0];

    int j = 0;
    int i = 1;
    while (i < n && j < k)
    {
        twoSum(eps, x[i], &r[j], &eps);
        if (eps == 0)
        { // no overflow
            eps = r[j];
        }
        else
        {
            j++;
        }
        i++;
    }
    if (eps != 0 && j < k)
    {
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

void mult2(double *x, double *y, double *pi, int n, int m, int r)
{

    double *B = (double *)alloca((r * dbl_prec / binSize + 2) * sizeof(double));
    // get sum of first exponents
    int e = exponent(x[0]) + exponent(y[0]);

    // initialize each Bin with starting value
    for (int i = 0; i < r * dbl_prec / binSize + 2; i++)
    {
        B[i] = ldexp(1.5, e - (i + 1) * binSize + dbl_prec - 1); // 1.5*2^(e-(i+1)b+p-1)
    }
    int j, l, sh;
    double p, err;
    for (int i = 0; i <= fmin(n - 1, r); i++)
    {
        for (j = 0; j <= fmin(m - 1, r - 1 - i); j++)
        {
            // p = two_prod(x[i], y[j], &err); this leads to more similar result as CAMPARY
            twoMultFMA(x[i], y[j], &p, &err);
            l = e - exponent(x[i]) - exponent(y[j]);
            sh = floor(l / binSize);      // bin of the first pair
            l = l - sh * binSize;         // number of leading bits
            accumulate(p, err, B, sh, l); // add to correct bins
        }
        j -= 1;
        if (j < m - 1)
        { //  I don't get what this part does
            p = x[i] * y[j];
            l = e - exponent(x[i]) - exponent(y[j]);
            sh = floor(l / binSize);
            l = l - sh * binSize;
            accumulate(p, 0., B, sh, l);
        }
    }
    for (int i = 0; i < r * dbl_prec / binSize + 2; i++)
    {
        B[i] = B[i] - ldexp(1.5, e - (i + 1) * binSize + dbl_prec - 1); // B_i - 1.5*2^(e-(i+1)b+p-1)
    }

    renormalize(B, pi, r * dbl_prec / binSize + 2, r);
}
