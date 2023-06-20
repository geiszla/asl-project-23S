#include <immintrin.h>
#include <cstring>
/* multiplies 4 doubles at a time
it's important that the input sizes of all doubles and the outputsizes are the same
no cehck for this consition is performed


*/
void fourtimesmultiplicationversion0(double *a0, double *b0, double *a1, double *b1, double *a2, double *b2, double *a3, double *b3, double *r0, double *r1, double *r2, double *r3, const int sizea, const int sizeb, const int sizer)
{
    multiplication(a0, b0, r0, sizea, sizeb, sizer);
    multiplication(a1, b1, r1, sizea, sizeb, sizer);
    multiplication(a2, b2, r2, sizea, sizeb, sizer);
    multiplication(a3, b3, r3, sizea, sizeb, sizer);
    return;
}

void fourtimesmultiplicationversion1(double *a0, double *b0, double *a1, double *b1, double *a2, double *b2, double *a3, double *b3, double *r0, double *r1, double *r2, double *r3, const int sizea, const int sizeb, const int sizer)
{

    int k = sizea;
    double *err = (double *)alloca(2 * (sizea * sizea + 3 * sizea) * sizeof(double));
    for (int i = 0; i < 2 * sizea * sizea + 3 * sizea; i++)
    {
        err[i] = 0;
    }
    double *r_ext = (double *)alloca(2 * (sizea * sizea + sizea) * sizeof(double));
    for (int i = 0; i < 2 * (sizea * sizea); i++)
    {
        r_ext[i] = 0.0;
    }

    twoMultFMA(a0[0], b0[0], &(r_ext[0]), &(err[0]));

    for (int n = 1; n <= (k - 1); n++)
    {
        double *e_tmp = (double *)alloca(2 * (n + 1) * sizeof(double));
        double *p = (double *)alloca(2 * (n + 1) * sizeof(double));
        for (int i = 0; i < 2 * (n + 1); i++)
        {
            e_tmp[i] = 0.0;
            p[i] = 0.0;
        }
        for (int i = 0; i <= n; i++)
        {
            twoMultFMA(a0[i], b0[n - i], &(p[i]), &(e_tmp[i]));
        }
        double *tmp = (double *)alloca((n * n + n + 1) * sizeof(double));
        double *tmp1 = (double *)alloca((n * n + n + 1) * sizeof(double));
        for (int i = 0; i < (n * n + n + 1); i++)
        {
            tmp1[i] = 0;
        }
        // generate p[0:n], e[0:n^2-1] into tmp
        for (int i = 0; i <= n; i++)
        {
            tmp[i] = p[i];
        }
        for (int i = 0; i <= (n * n - 1); i++)
        {
            tmp[n + 1 + i] = err[i];
        }
        vecSum(tmp, tmp1, (n * n + n));
        /* write  tmp1 into r_ext[n], e[0:n^2 +n-1]*/
        r_ext[n] = tmp1[0];
        for (int i = 0; i <= (n * n + n - 1); i++)
        {
            err[i] = tmp1[i + 1];
        }
        // writes err[0:n^2 +n-1],e_tmp[0:n] into err[0:n^2 +2n-1]
        for (int i = (n * n); i <= ((n ^ 2) + 2 * n - 1); i++)
        {
            err[i] = e_tmp[i - n * n];
        }
    }

    for (int i = 1; i <= (k - 1); i++)
    {

        r_ext[k] += a0[i] * b0[k - i];
    }

    for (int i = 0; i <= (k * k - 1); i++)
    {

        r_ext[i] += err[i];
    }

    renormalizationalgorithm(r_ext, k + 1, r0, sizea);

    // multiplication(a1, b1, r1, sizea, sizeb, sizer); -------------------------------------------------------------

    k = sizea;
    err = (double *)alloca(2 * (sizea * sizea + 3 * sizea) * sizeof(double));
    for (int i = 0; i < 2 * sizea * sizea + 3 * sizea; i++)
    {
        err[i] = 0;
    }
    r_ext = (double *)alloca(2 * (sizea * sizea + sizea) * sizeof(double));
    for (int i = 0; i < 2 * (sizea * sizea); i++)
    {
        r_ext[i] = 0.0;
    }

    twoMultFMA(a1[0], b1[0], &(r_ext[0]), &(err[0]));

    for (int n = 1; n <= (k - 1); n++)
    {
        double *e_tmp = (double *)alloca(2 * (n + 1) * sizeof(double));
        double *p = (double *)alloca(2 * (n + 1) * sizeof(double));
        for (int i = 0; i < 2 * (n + 1); i++)
        {
            e_tmp[i] = 0.0;
            p[i] = 0.0;
        }
        for (int i = 0; i <= n; i++)
        {
            twoMultFMA(a1[i], b1[n - i], &(p[i]), &(e_tmp[i]));
        }
        double *tmp = (double *)alloca((n * n + n + 1) * sizeof(double));
        double *tmp1 = (double *)alloca((n * n + n + 1) * sizeof(double));
        for (int i = 0; i < (n * n + n + 1); i++)
        {
            tmp1[i] = 0;
        }
        // generate p[0:n], e[0:n^2-1] into tmp
        for (int i = 0; i <= n; i++)
        {
            tmp[i] = p[i];
        }
        for (int i = 0; i <= (n * n - 1); i++)
        {
            tmp[n + 1 + i] = err[i];
        }
        vecSum(tmp, tmp1, (n * n + n));
        /* write  tmp1 into r_ext[n], e[0:n^2 +n-1]*/
        r_ext[n] = tmp1[0];
        for (int i = 0; i <= (n * n + n - 1); i++)
        {
            err[i] = tmp1[i + 1];
        }
        // writes err[0:n^2 +n-1],e_tmp[0:n] into err[0:n^2 +2n-1]
        for (int i = (n * n); i <= ((n ^ 2) + 2 * n - 1); i++)
        {
            err[i] = e_tmp[i - n * n];
        }
    }

    for (int i = 1; i <= (k - 1); i++)
    {

        r_ext[k] += a1[i] * b1[k - i];
    }

    for (int i = 0; i <= (k * k - 1); i++)
    {

        r_ext[i] += err[i];
    }

    renormalizationalgorithm(r_ext, k + 1, r1, sizea);
    // multiplication(a2, b2, r2, sizea, sizeb, sizer); ----------------------------------------------------------------------------------------------------------
    k = sizea;
    err = (double *)alloca(2 * (sizea * sizea + 3 * sizea) * sizeof(double));
    for (int i = 0; i < 2 * sizea * sizea + 3 * sizea; i++)
    {
        err[i] = 0;
    }
    r_ext = (double *)alloca(2 * (sizea * sizea + sizea) * sizeof(double));
    for (int i = 0; i < 2 * (sizea * sizea); i++)
    {
        r_ext[i] = 0.0;
    }

    twoMultFMA(a2[0], b2[0], &(r_ext[0]), &(err[0]));

    for (int n = 1; n <= (k - 1); n++)
    {
        double *e_tmp = (double *)alloca(2 * (n + 1) * sizeof(double));
        double *p = (double *)alloca(2 * (n + 1) * sizeof(double));
        for (int i = 0; i < 2 * (n + 1); i++)
        {
            e_tmp[i] = 0.0;
            p[i] = 0.0;
        }
        for (int i = 0; i <= n; i++)
        {
            twoMultFMA(a2[i], b2[n - i], &(p[i]), &(e_tmp[i]));
        }
        double *tmp = (double *)alloca((n * n + n + 1) * sizeof(double));
        double *tmp1 = (double *)alloca((n * n + n + 1) * sizeof(double));
        for (int i = 0; i < (n * n + n + 1); i++)
        {
            tmp1[i] = 0;
        }
        // generate p[0:n], e[0:n^2-1] into tmp
        for (int i = 0; i <= n; i++)
        {
            tmp[i] = p[i];
        }
        for (int i = 0; i <= (n * n - 1); i++)
        {
            tmp[n + 1 + i] = err[i];
        }
        vecSum(tmp, tmp1, (n * n + n));
        /* write  tmp1 into r_ext[n], e[0:n^2 +n-1]*/
        r_ext[n] = tmp1[0];
        for (int i = 0; i <= (n * n + n - 1); i++)
        {
            err[i] = tmp1[i + 1];
        }
        // writes err[0:n^2 +n-1],e_tmp[0:n] into err[0:n^2 +2n-1]
        for (int i = (n * n); i <= ((n ^ 2) + 2 * n - 1); i++)
        {
            err[i] = e_tmp[i - n * n];
        }
    }

    for (int i = 1; i <= (k - 1); i++)
    {

        r_ext[k] += a2[i] * b2[k - i];
    }

    for (int i = 0; i <= (k * k - 1); i++)
    {

        r_ext[i] += err[i];
    }

    renormalizationalgorithm(r_ext, k + 1, r2, sizea);

    // multiplication(a3, b3, r3, sizea, sizeb, sizer); -----------------------------------------------------------------------
    k = sizea;
    err = (double *)alloca(2 * (sizea * sizea + 3 * sizea) * sizeof(double));
    for (int i = 0; i < 2 * sizea * sizea + 3 * sizea; i++)
    {
        err[i] = 0;
    }
    r_ext = (double *)alloca(2 * (sizea * sizea + sizea) * sizeof(double));
    for (int i = 0; i < 2 * (sizea * sizea); i++)
    {
        r_ext[i] = 0.0;
    }

    twoMultFMA(a1[0], b1[0], &(r_ext[0]), &(err[0]));

    for (int n = 1; n <= (k - 1); n++)
    {
        double *e_tmp = (double *)alloca(2 * (n + 1) * sizeof(double));
        double *p = (double *)alloca(2 * (n + 1) * sizeof(double));
        for (int i = 0; i < 2 * (n + 1); i++)
        {
            e_tmp[i] = 0.0;
            p[i] = 0.0;
        }
        for (int i = 0; i <= n; i++)
        {
            twoMultFMA(a3[i], b3[n - i], &(p[i]), &(e_tmp[i]));
        }
        double *tmp = (double *)alloca((n * n + n + 1) * sizeof(double));
        double *tmp1 = (double *)alloca((n * n + n + 1) * sizeof(double));
        for (int i = 0; i < (n * n + n + 1); i++)
        {
            tmp1[i] = 0;
        }
        // generate p[0:n], e[0:n^2-1] into tmp
        for (int i = 0; i <= n; i++)
        {
            tmp[i] = p[i];
        }
        for (int i = 0; i <= (n * n - 1); i++)
        {
            tmp[n + 1 + i] = err[i];
        }
        vecSum(tmp, tmp1, (n * n + n));
        /* write  tmp1 into r_ext[n], e[0:n^2 +n-1]*/
        r_ext[n] = tmp1[0];
        for (int i = 0; i <= (n * n + n - 1); i++)
        {
            err[i] = tmp1[i + 1];
        }
        // writes err[0:n^2 +n-1],e_tmp[0:n] into err[0:n^2 +2n-1]
        for (int i = (n * n); i <= ((n ^ 2) + 2 * n - 1); i++)
        {
            err[i] = e_tmp[i - n * n];
        }
    }

    for (int i = 1; i <= (k - 1); i++)
    {

        r_ext[k] += a3[i] * b3[k - i];
    }

    for (int i = 0; i <= (k * k - 1); i++)
    {

        r_ext[i] += err[i];
    }

    renormalizationalgorithm(r_ext, k + 1, r3, sizea);
    return;
}

void fourtimesmultiplicationversion3(double *a0, double *b0, double *a1, double *b1, double *a2, double *b2, double *a3, double *b3, double *r0, double *r1, double *r2, double *r3, const int sizea, const int sizeb, const int sizer)
{

    int k = sizea;

    double *err0 = (double *)_mm_malloc(2 * (sizea * sizea + 3 * sizea) * sizeof(double), 64);
    double *err1 = (double *)_mm_malloc(2 * (sizea * sizea + 3 * sizea) * sizeof(double), 64);
    double *err2 = (double *)_mm_malloc(2 * (sizea * sizea + 3 * sizea) * sizeof(double), 64);
    double *err3 = (double *)_mm_malloc(2 * (sizea * sizea + 3 * sizea) * sizeof(double), 64);

    __m256d zero = _mm256_setzero_pd();
    __m256d zero1 = _mm256_setzero_pd();
    for (int i = 0; i < 2 * sizea * sizea + 3 * sizea; i += 4)
    {
        _mm256_storeu_pd(&err0[i], zero);
        _mm256_storeu_pd(&err1[i], zero1);
        _mm256_storeu_pd(&err2[i], zero);
        _mm256_storeu_pd(&err3[i], zero1);
    }

    double *r_ext0 = (double *)_mm_malloc(2 * (sizea * sizea + sizea) * sizeof(double), 64);
    double *r_ext1 = (double *)_mm_malloc(2 * (sizea * sizea + sizea) * sizeof(double), 64);
    double *r_ext2 = (double *)_mm_malloc(2 * (sizea * sizea + sizea) * sizeof(double), 64);
    double *r_ext3 = (double *)_mm_malloc(2 * (sizea * sizea + sizea) * sizeof(double), 64);

    for (int i = 0; i < 2 * (sizea * sizea); i += 4)
    {
        _mm256_storeu_pd(&r_ext0[i], zero);
        _mm256_storeu_pd(&r_ext1[i], zero1);
        _mm256_storeu_pd(&r_ext2[i], zero);
        _mm256_storeu_pd(&r_ext3[i], zero1);
    }

    // twoMultFMA(a0[0], b0[0], &(r_ext0[0]), &(err0[0]));
    // twoMultFMA(a1[0], b1[0], &(r_ext1[0]), &(err1[0]));
    // twoMultFMA(a2[0], b2[0], &(r_ext2[0]), &(err2[0]));
    // twoMultFMA(a3[0], b3[0], &(r_ext3[0]), &(err3[0]));

    double pi0 = a0[0] * b0[0];
    double pi1 = a1[0] * b1[0];
    double pi2 = a2[0] * b2[0];
    double pi3 = a3[0] * b3[0];
    double e0 = fma(a0[0], b0[0], -pi0);
    double e1 = fma(a1[0], b1[0], -pi1);
    double e2 = fma(a2[0], b2[0], -pi2);
    double e3 = fma(a3[0], b3[0], -pi3);
    r_ext0[0] = pi0;
    r_ext1[0] = pi1;
    r_ext2[0] = pi2;
    r_ext3[0] = pi3;
    err0[0] = e0;
    err1[0] = e1;
    err2[0] = e2;
    err3[0] = e3;

    for (int n = 1; n <= (k - 1); n++)
    {
        double *e_tmp0 = (double *)alloca(2 * (n + 1) * sizeof(double));
        double *p0 = (double *)alloca(2 * (n + 1) * sizeof(double));
        double *e_tmp1 = (double *)alloca(2 * (n + 1) * sizeof(double));
        double *p1 = (double *)alloca(2 * (n + 1) * sizeof(double));
        double *e_tmp2 = (double *)alloca(2 * (n + 1) * sizeof(double));
        double *p2 = (double *)alloca(2 * (n + 1) * sizeof(double));
        double *e_tmp3 = (double *)alloca(2 * (n + 1) * sizeof(double));
        double *p3 = (double *)alloca(2 * (n + 1) * sizeof(double));

        for (int i = 0; i < 2 * (n + 1); i++)
        {
            e_tmp0[i] = 0.0;
            p0[i] = 0.0;
            e_tmp1[i] = 0.0;
            p1[i] = 0.0;
            e_tmp2[i] = 0.0;
            p2[i] = 0.0;
            e_tmp3[i] = 0.0;
            p3[i] = 0.0;
        }
        for (int i = 0; i <= n; i++)
        {
            
            double a[4] = {a0[i], a1[i], a2[i], a3[i]};
            double b[4] = {b0[n - i], b1[n - i], b2[n - i], b3[n - i]};
            double p[4] = {};
            double e[4] = {};
            __m256d av = _mm256_load_pd(a);
            __m256d bv = _mm256_load_pd(b);
            __m256d pv = _mm256_mul_pd(av, bv);
            __m256d ev = _mm256_fmsub_pd(av, bv, pv);
            _mm256_store_pd(p, pv);
            _mm256_store_pd(e, ev);
            p0[i] = p[0];
            p1[i] = p[1];
            p2[i] = p[2];
            p3[i] = p[3];
            e_tmp0[i] = e[0];
            e_tmp1[i] = e[1];
            e_tmp2[i] = e[2];
            e_tmp3[i] = e[3];


        }
        double *tmp_0 = (double *)alloca((n * n + n + 1) * sizeof(double));
        double *tmp1_0 = (double *)alloca((n * n + n + 1) * sizeof(double));
        double *tmp_1 = (double *)alloca((n * n + n + 1) * sizeof(double));
        double *tmp1_1 = (double *)alloca((n * n + n + 1) * sizeof(double));
        double *tmp_2 = (double *)alloca((n * n + n + 1) * sizeof(double));
        double *tmp1_2 = (double *)alloca((n * n + n + 1) * sizeof(double));
        double *tmp_3 = (double *)alloca((n * n + n + 1) * sizeof(double));
        double *tmp1_3 = (double *)alloca((n * n + n + 1) * sizeof(double));

        tmp1_0[n * n + n] = 0;
        tmp1_1[n * n + n] = 0;
        tmp1_2[n * n + n] = 0;
        tmp1_3[n * n + n] = 0;

        // generate p[0:n], e[0:n^2-1] into tmp
        int i = 0;
        for (; i <= n - 4; i += 4)
        {
            __m256d p0_vec = _mm256_load_pd(&p0[i]);
            __m256d p1_vec = _mm256_load_pd(&p1[i]);
            __m256d p2_vec = _mm256_load_pd(&p2[i]);
            __m256d p3_vec = _mm256_load_pd(&p3[i]);

            _mm256_store_pd(&tmp_0[i], p0_vec);
            _mm256_store_pd(&tmp_1[i], p1_vec);
            _mm256_store_pd(&tmp_2[i], p2_vec);
            _mm256_store_pd(&tmp_3[i], p3_vec);
        }
        for (; i <= n; i++)
        {
            tmp_0[i] = p0[i];
            tmp_1[i] = p1[i];
            tmp_2[i] = p2[i];
            tmp_3[i] = p3[i];
        }
        i = 0;
        for (; i <= (n * n - 1) - 5; i += 4)
        {
            __m256d err0_vec = _mm256_load_pd(&err0[i]);
            __m256d err1_vec = _mm256_load_pd(&err1[i]);
            __m256d err2_vec = _mm256_load_pd(&err2[i]);
            __m256d err3_vec = _mm256_load_pd(&err3[i]);

            _mm256_store_pd(&tmp_0[n + 1 + i], err0_vec);
            _mm256_store_pd(&tmp_1[n + 1 + i], err1_vec);
            _mm256_store_pd(&tmp_2[n + 1 + i], err2_vec);
            _mm256_store_pd(&tmp_3[n + 1 + i], err3_vec);
        }
        for (; i <= (n * n - 1); i++)
        {
            tmp_0[n + 1 + i] = err0[i];
            tmp_1[n + 1 + i] = err1[i];
            tmp_2[n + 1 + i] = err2[i];
            tmp_3[n + 1 + i] = err3[i];
        }
  

        int length = (n * n + n);
        double *s0 = (double *)alloca(length * sizeof(double));
        double *s1 = (double *)alloca(length * sizeof(double));
        double *s2 = (double *)alloca(length * sizeof(double));
        double *s3 = (double *)alloca(length * sizeof(double));
        s0[length - 1] = tmp_0[length - 1];
        s1[length - 1] = tmp_1[length - 1];
        s2[length - 1] = tmp_2[length - 1];
        s3[length - 1] = tmp_3[length - 1];
        for (int i = length - 2; i >= 0; i--)
        {
            double s_tmp0, e_tmp0;
            double s_tmp1, e_tmp1;
            double s_tmp2, e_tmp2;
            double s_tmp3, e_tmp3;
            double ssum0 = tmp_0[i] + s0[i + 1];
            double ssum1 = tmp_1[i] + s1[i + 1];
            double ssum2 = tmp_2[i] + s2[i + 1];
            double ssum3 = tmp_3[i] + s3[i + 1];
            double t0 = ssum0 - s0[i + 1];
            double t1 = ssum1 - s1[i + 1];
            double t2 = ssum2 - s2[i + 1];
            double t3 = ssum3 - s3[i + 1];

            e_tmp0 = (tmp_0[i] - t0) + (s0[i + 1] - (ssum0 - t0));
            e_tmp1 = (tmp_1[i] - t1) + (s1[i + 1] - (ssum1 - t1));
            e_tmp2 = (tmp_2[i] - t2) + (s2[i + 1] - (ssum2 - t2));
            e_tmp3 = (tmp_3[i] - t3) + (s3[i + 1] - (ssum3 - t3));

            s_tmp0 = ssum0;
            s_tmp1 = ssum1;
            s_tmp2 = ssum2;
            s_tmp3 = ssum3;
            
            s0[i] = s_tmp0;
            s1[i] = s_tmp1;
            s2[i] = s_tmp2;
            s3[i] = s_tmp3;
            tmp1_0[i + 1] = e_tmp0;
            tmp1_1[i + 1] = e_tmp1;
            tmp1_2[i + 1] = e_tmp2;
            tmp1_3[i + 1] = e_tmp3;
        }


        
        tmp1_0[0] = s0[0];
        tmp1_1[0] = s1[0];
        tmp1_2[0] = s2[0];
        tmp1_3[0] = s3[0];


        /* write  tmp1 into r_ext[n], e[0:n^2 +n-1]*/
        r_ext0[n] = tmp1_0[0];
        r_ext1[n] = tmp1_1[0];
        r_ext2[n] = tmp1_2[0];
        r_ext3[n] = tmp1_3[0];
        i = 0;
        for (; i <= (n * n + n - 1) - 5; i += 4)
        {
            __m256d tmp1_0_vec = _mm256_loadu_pd(&tmp1_0[i + 1]);
            __m256d tmp1_1_vec = _mm256_loadu_pd(&tmp1_1[i + 1]);
            __m256d tmp1_2_vec = _mm256_loadu_pd(&tmp1_2[i + 1]);
            __m256d tmp1_3_vec = _mm256_loadu_pd(&tmp1_3[i + 1]);

            _mm256_storeu_pd(&err0[i], tmp1_0_vec);
            _mm256_storeu_pd(&err1[i], tmp1_1_vec);
            _mm256_storeu_pd(&err2[i], tmp1_2_vec);
            _mm256_storeu_pd(&err3[i], tmp1_3_vec);
        }
        for (; i <= (n * n + n - 1); i++)
        {
            err0[i] = tmp1_0[i + 1];
            err1[i] = tmp1_1[i + 1];
            err2[i] = tmp1_2[i + 1];
            err3[i] = tmp1_3[i + 1];
        }
        // writes err[0:n^2 +n-1],e_tmp[0:n] into err[0:n^2 +2n-1]
        for (int i = (n * n); i <= ((n ^ 2) + 2 * n - 1); i++)
        {
            err0[i] = e_tmp0[i - n * n];
            err1[i] = e_tmp1[i - n * n];
            err2[i] = e_tmp2[i - n * n];
            err3[i] = e_tmp3[i - n * n];
        }
    }

    for (int i = 1; i <= (k - 1); i++)
    {

        r_ext0[k] += a0[i] * b0[k - i];
        r_ext1[k] += a1[i] * b1[k - i];
        r_ext2[k] += a2[i] * b2[k - i];
        r_ext3[k] += a3[i] * b3[k - i];
    }

    for (int i = 0; i <= (k * k - 1); i++)
    {

        r_ext0[i] += err0[i];
        r_ext1[i] += err1[i];
        r_ext2[i] += err2[i];
        r_ext3[i] += err3[i];
    }

    renormalizationalgorithm(r_ext0, k + 1, r0, sizea);
    renormalizationalgorithm(r_ext1, k + 1, r1, sizea);
    renormalizationalgorithm(r_ext2, k + 1, r2, sizea);
    renormalizationalgorithm(r_ext3, k + 1, r3, sizea);

    // multiplication(a1, b1, r1, sizea, sizeb, sizer); -------------------------------------------------------------
    _mm_free(err0);
    _mm_free(err1);
    _mm_free(err2);
    _mm_free(err3);
    _mm_free(r_ext0);
    _mm_free(r_ext1);
    _mm_free(r_ext2);
    _mm_free(r_ext3);
}
