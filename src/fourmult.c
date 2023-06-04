#include <immintrin.h>
/* multiplies 4 doubles at a time 
it's important that the input sizes of all doubles and the outputsizes are the same 
no cehck for this consition is performed*/
void fourtimesmultiplicationversion0(double *a0, double *b0, double *a1, double *b1, double *a2, double *b2,double *a3, double *b3,double *r0, double *r1, double *r2, double *r3, const int sizea, const int sizeb, const int sizer)
{
    multiplication(a0,b0,r0,sizea,sizeb,sizer);
    multiplication(a1,b1,r1,sizea,sizeb,sizer);
    multiplication(a2,b2,r2,sizea,sizeb,sizer);
    multiplication(a3,b3,r3,sizea,sizeb,sizer);
    return;
}
// unrolled version of fourtimesmultiplicationversion0
void fourtimesmultiplicationversion1(double *a0, double *b0, double *a1, double *b1, double *a2, double *b2,double *a3, double *b3,double *r0, double *r1, double *r2, double *r3, const int sizea, const int sizeb, const int sizer)
{  
    int k = sizea;
    double *err0 = (double *)calloc((sizea * sizea + 2 * sizea+1), sizeof(double));
    double *r_ext0 = (double *)calloc((sizea * sizea), sizeof(double));
    double *err1 = (double *)calloc((sizea * sizea + 2 * sizea+1), sizeof(double));
    double *r_ext1 = (double *)calloc((sizea * sizea), sizeof(double));
    double *err2 = (double *)calloc((sizea * sizea + 2 * sizea+1), sizeof(double));
    double *r_ext2 = (double *)calloc((sizea * sizea), sizeof(double));
    double *err3 = (double *)calloc((sizea * sizea + 2 * sizea+1), sizeof(double));
    double *r_ext3 = (double *)calloc((sizea * sizea), sizeof(double));

    twoMultFMA(a0[0], b0[0], &(r_ext0[0]), &(err0[0]));
    twoMultFMA(a1[0], b1[0], &(r_ext1[0]), &(err1[0]));
    twoMultFMA(a2[0], b2[0], &(r_ext2[0]), &(err2[0]));
    twoMultFMA(a3[0], b3[0], &(r_ext3[0]), &(err3[0]));

    for (int n = 1; n <= (k - 1); n++)
    {
        double *e_tmp0 = (double *)calloc((n+1), sizeof(double));
        double *p0 = (double *)calloc((n + 1), sizeof(double));
        double *e_tmp1 = (double *)calloc((n+1), sizeof(double));
        double *p1 = (double *)calloc((n + 1), sizeof(double));
        double *e_tmp2 = (double *)calloc((n+1), sizeof(double));
        double *p2 = (double *)calloc((n + 1), sizeof(double));
        double *e_tmp3 = (double *)calloc((n+1), sizeof(double));
        double *p3 = (double *)calloc((n + 1), sizeof(double));

        for (int i = 0; i <= n; i++)
        {
            twoMultFMA(a0[i], b0[n - i], &(p0[i]), &(e_tmp0[i]));
            twoMultFMA(a1[i], b1[n - i], &(p1[i]), &(e_tmp1[i]));
            twoMultFMA(a2[i], b2[n - i], &(p2[i]), &(e_tmp2[i]));
            twoMultFMA(a3[i], b3[n - i], &(p3[i]), &(e_tmp3[i]));
            //------------------------------------------------
            
            //-------------------------------- twoMultFMa inline
        }
        double *tmp0 = &err0[-n - 1];
        double *tmp1_0 = (double *)calloc((n * n + n+1), sizeof(double));
        double *tmp1 = &err1[-n - 1];
        double *tmp1_1 = (double *)calloc((n * n + n+1), sizeof(double));
        double *tmp2 = &err2[-n - 1];
        double *tmp1_2 = (double *)calloc((n * n + n+1), sizeof(double));
        double *tmp3 = &err3[-n - 1];
        double *tmp1_3 = (double *)calloc((n * n + n+1), sizeof(double));
        


        //-------------------------------- vecsum inline
        int length = (n * n + n);
        double *s0 = (double *)alloca(length * sizeof(double));
        double *s1 = (double *)alloca(length * sizeof(double));
        double *s2 = (double *)alloca(length * sizeof(double));
        double *s3 = (double *)alloca(length * sizeof(double));
        s0[length - 1] = err0[n * n - 1];
        s1[length - 1] = err1[n * n - 1];
        s2[length - 1] = err2[n * n - 1];
        s3[length - 1] = err3[n * n - 1];
        for (int i = length - 2; i >= 0; i--)
        {
            if (i <= (n))
            {
                tmp0 = p0;
                tmp1 = p1;
                tmp2 = p2;
                tmp3 = p3;
            }
            double s_tmp0, e_tmplocal0;
            double s_tmp1, e_tmplocal1;
            double s_tmp2, e_tmplocal2;
            double s_tmp3, e_tmplocal3;
            
            double sl0 = tmp0[i] + s0[i + 1];
            double sl1 = tmp1[i] + s1[i + 1];
            double sl2 = tmp2[i] + s2[i + 1];
            double sl3 = tmp3[i] + s3[i + 1];

            double t0 = sl0 - s0[i + 1];
            double t1 = sl1 - s1[i + 1];
            double t2 = sl2 - s2[i + 1];
            double t3 = sl3 - s3[i + 1];
            double e0 = (tmp0[i] - t0) + (s0[i + 1] - (sl0 - t0));
            double e1 = (tmp1[i] - t1) + (s1[i + 1] - (sl1 - t1));
            double e2 = (tmp2[i] - t2) + (s2[i + 1] - (sl2 - t2));
            double e3 = (tmp3[i] - t3) + (s3[i + 1] - (sl3 - t3));
            s_tmp0 = sl0;
            s_tmp1 = sl1;
            s_tmp2 = sl2;
            s_tmp3 = sl3;
            e_tmplocal0 = e0;
            e_tmplocal1 = e1;
            e_tmplocal2 = e2;
            e_tmplocal3 = e3;
            //--------------------------------
            s0[i] = s_tmp0;
            s1[i] = s_tmp1;
            s2[i] = s_tmp2;
            s3[i] = s_tmp3;
            tmp1_0[i + 1] = e_tmplocal0;
            tmp1_1[i + 1] = e_tmplocal1;
            tmp1_2[i + 1] = e_tmplocal2;
            tmp1_3[i + 1] = e_tmplocal3;
        }
        r_ext0[n] = s0[0];
        r_ext1[n] = s1[0];
        r_ext2[n] = s2[0];
        r_ext3[n] = s3[0];

        //--------------------------------
        /* write  tmp1 into r_ext[n], e[0:n^2 +n-1]*/
       
        for (int i = 0; i <= (n * n + n - 1); i++)
        {
            err0[i] = tmp1_0[i + 1];
            err1[i] = tmp1_1[i + 1];
            err2[i] = tmp1_2[i + 1];
            err3[i] = tmp1_3[i + 1];
        }
        // writes err[0:n^2 +n-1],e_tmp[0:n] into err[0:n^2 +2n-1]
        for (int i = (n * n+n); i <= (n*n + 2*n -1); i++)
        {
            err0[i] = e_tmp0[i - (n * n+n) ];
            err1[i] = e_tmp1[i - (n * n+n) ];
            err2[i] = e_tmp2[i - (n * n+n) ];
            err3[i] = e_tmp3[i - (n * n+n) ];
           
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
    return;
}



void fourtimesmultiplicationversion2(double *a0, double *b0, double *a1, double *b1, double *a2, double *b2,double *a3, double *b3,double *r0, double *r1, double *r2, double *r3, const int sizea, const int sizeb, const int sizer)
{  
    int k = sizea;
    double *err0 = (double *)calloc((sizea * sizea + 2 * sizea+1), sizeof(double));
    double *r_ext0 = (double *)calloc((sizea * sizea), sizeof(double));
    double *err1 = (double *)calloc((sizea * sizea + 2 * sizea+1), sizeof(double));
    double *r_ext1 = (double *)calloc((sizea * sizea), sizeof(double));
    double *err2 = (double *)calloc((sizea * sizea + 2 * sizea+1), sizeof(double));
    double *r_ext2 = (double *)calloc((sizea * sizea), sizeof(double));
    double *err3 = (double *)calloc((sizea * sizea + 2 * sizea+1), sizeof(double));
    double *r_ext3 = (double *)calloc((sizea * sizea), sizeof(double));

    twoMultFMA(a0[0], b0[0], &(r_ext0[0]), &(err0[0]));
    twoMultFMA(a1[0], b1[0], &(r_ext1[0]), &(err1[0]));
    twoMultFMA(a2[0], b2[0], &(r_ext2[0]), &(err2[0]));
    twoMultFMA(a3[0], b3[0], &(r_ext3[0]), &(err3[0]));


    for (int n = 1; n <= (k - 1); n++)
    {
        double *e_tmp0 = (double *)calloc((n+1), sizeof(double));
        double *p0 = (double *)calloc((n + 1), sizeof(double));
        double *e_tmp1 = (double *)calloc((n+1), sizeof(double));
        double *p1 = (double *)calloc((n + 1), sizeof(double));
        double *e_tmp2 = (double *)calloc((n+1), sizeof(double));
        double *p2 = (double *)calloc((n + 1), sizeof(double));
        double *e_tmp3 = (double *)calloc((n+1), sizeof(double));
        double *p3 = (double *)calloc((n + 1), sizeof(double));

        for (int i = 0; i <= n; i++)
        {
            twoMultFMA(a0[i], b0[n - i], &(p0[i]), &(e_tmp0[i]));
            twoMultFMA(a1[i], b1[n - i], &(p1[i]), &(e_tmp1[i]));
            twoMultFMA(a2[i], b2[n - i], &(p2[i]), &(e_tmp2[i]));
            twoMultFMA(a3[i], b3[n - i], &(p3[i]), &(e_tmp3[i]));
         
        }
        double *tmp0 = &err0[-n - 1];
        double *tmp1_0 = (double *)calloc((n * n + n+1), sizeof(double));
        double *tmp1 = &err1[-n - 1];
        double *tmp1_1 = (double *)calloc((n * n + n+1), sizeof(double));
        double *tmp2 = &err2[-n - 1];
        double *tmp1_2 = (double *)calloc((n * n + n+1), sizeof(double));
        double *tmp3 = &err3[-n - 1];
        double *tmp1_3 = (double *)calloc((n * n + n+1), sizeof(double));
        


        //-------------------------------- vecsum inline
        int length = (n * n + n);
        double *s0 = (double *)alloca(length * sizeof(double));
        double *s1 = (double *)alloca(length * sizeof(double));
        double *s2 = (double *)alloca(length * sizeof(double));
        double *s3 = (double *)alloca(length * sizeof(double));
        s0[length - 1] = err0[n * n - 1];
        s1[length - 1] = err1[n * n - 1];
        s2[length - 1] = err2[n * n - 1];
        s3[length - 1] = err3[n * n - 1];
        for (int i = length - 2; i >= 0; i--)
        {
            if (i <= (n))
            {
                tmp0 = p0;
                tmp1 = p1;
                tmp2 = p2;
                tmp3 = p3;
            }
            double s_tmp0, e_tmplocal0;
            double s_tmp1, e_tmplocal1;
            double s_tmp2, e_tmplocal2;
            double s_tmp3, e_tmplocal3;
            
            double sl0 = tmp0[i] + s0[i + 1];
            double sl1 = tmp1[i] + s1[i + 1];
            double sl2 = tmp2[i] + s2[i + 1];
            double sl3 = tmp3[i] + s3[i + 1];

            double t0 = sl0 - s0[i + 1];
            double t1 = sl1 - s1[i + 1];
            double t2 = sl2 - s2[i + 1];
            double t3 = sl3 - s3[i + 1];
            double e0 = (tmp0[i] - t0) + (s0[i + 1] - (sl0 - t0));
            double e1 = (tmp1[i] - t1) + (s1[i + 1] - (sl1 - t1));
            double e2 = (tmp2[i] - t2) + (s2[i + 1] - (sl2 - t2));
            double e3 = (tmp3[i] - t3) + (s3[i + 1] - (sl3 - t3));
            s_tmp0 = sl0;
            s_tmp1 = sl1;
            s_tmp2 = sl2;
            s_tmp3 = sl3;
            e_tmplocal0 = e0;
            e_tmplocal1 = e1;
            e_tmplocal2 = e2;
            e_tmplocal3 = e3;
            //--------------------------------
            s0[i] = s_tmp0;
            s1[i] = s_tmp1;
            s2[i] = s_tmp2;
            s3[i] = s_tmp3;
            tmp1_0[i + 1] = e_tmplocal0;
            tmp1_1[i + 1] = e_tmplocal1;
            tmp1_2[i + 1] = e_tmplocal2;
            tmp1_3[i + 1] = e_tmplocal3;
        }
        r_ext0[n] = s0[0];
        r_ext1[n] = s1[0];
        r_ext2[n] = s2[0];
        r_ext3[n] = s3[0];

        //--------------------------------
        /* write  tmp1 into r_ext[n], e[0:n^2 +n-1]*/
       
        for (int i = 0; i <= (n * n + n - 1); i++)
        {
            err0[i] = tmp1_0[i + 1];
            err1[i] = tmp1_1[i + 1];
            err2[i] = tmp1_2[i + 1];
            err3[i] = tmp1_3[i + 1];
        }
        // writes err[0:n^2 +n-1],e_tmp[0:n] into err[0:n^2 +2n-1]
        for (int i = (n * n+n); i <= (n*n + 2*n -1); i++)
        {
            err0[i] = e_tmp0[i - (n * n+n) ];
            err1[i] = e_tmp1[i - (n * n+n) ];
            err2[i] = e_tmp2[i - (n * n+n) ];
            err3[i] = e_tmp3[i - (n * n+n) ];
           
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
    return;
}
