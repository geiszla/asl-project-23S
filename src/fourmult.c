
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
{  int k = sizea;
    double *err0 = (double *)alloca(2*(sizea * sizea +1*sizea +1) * sizeof(double));
    double *r_ext0 = (double *)alloca((sizea  + 1) * sizeof(double));
    double *err1 = (double *)alloca(2*(sizea * sizea +1*sizea +1) * sizeof(double));
    double *r_ext1 = (double *)alloca((sizea  + 1) * sizeof(double));
    double *err2 = (double *)alloca(2*(sizea * sizea +1*sizea +1) * sizeof(double));
    double *r_ext2 = (double *)alloca((sizea  + 1) * sizeof(double));
    double *err3 = (double *)alloca(2*(sizea * sizea +1*sizea +1) * sizeof(double));
    double *r_ext3 = (double *)alloca((sizea  + 1) * sizeof(double));
   
    

    twoMultFMA(a0[0], b0[0], &(r_ext0[0]), &(err0[0]));
    twoMultFMA(a1[0], b1[0], &(r_ext1[0]), &(err1[0]));
    twoMultFMA(a2[0], b2[0], &(r_ext2[0]), &(err2[0]));
    twoMultFMA(a3[0], b3[0], &(r_ext3[0]), &(err3[0]));

    for (int n = 1; n <= (k - 1); n++)
    {
        double *e_tmp0 = (double *)alloca((n+1) * sizeof(double));
        double *p0 = (double *)alloca((n+1)  * sizeof(double));
        double *e_tmp1 = (double *)alloca((n+1) * sizeof(double));
        double *p1 = (double *)alloca((n+1)  * sizeof(double));
        double *e_tmp2 = (double *)alloca((n+1) * sizeof(double));
        double *p2 = (double *)alloca((n+1)  * sizeof(double));
        double *e_tmp3 = (double *)alloca((n+1) * sizeof(double));
        double *p3 = (double *)alloca((n+1)  * sizeof(double));

        for (int i = 0; i < (n+1); i++)
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
            twoMultFMA(a0[i], b0[n - i], &(p0[i]), &(e_tmp0[i]));
            twoMultFMA(a1[i], b1[n - i], &(p1[i]), &(e_tmp1[i]));
            twoMultFMA(a2[i], b2[n - i], &(p2[i]), &(e_tmp2[i]));
            twoMultFMA(a3[i], b3[n - i], &(p3[i]), &(e_tmp3[i]));
        }
        double *tmp_0 = (double *)alloca((n * n + 3*n) * sizeof(double));
        double *tmp1_0 = (double *)alloca(3*n * sizeof(double));
        double *tmp_1 = (double *)alloca((n * n + 3*n) * sizeof(double));
        double *tmp1_1 = (double *)alloca(3*n * sizeof(double));
        double *tmp_2 = (double *)alloca((n * n + 3*n) * sizeof(double));
        double *tmp1_2 = (double *)alloca(3*n * sizeof(double));
        double *tmp_3 = (double *)alloca((n * n + 3*n) * sizeof(double));
        double *tmp1_3 = (double *)alloca(3*n * sizeof(double));
        for (int i = 0; i <= 3*n; i++)
        {
            tmp1_0[i] = 0.0;
            tmp1_1[i] = 0.0;
            tmp1_2[i] = 0.0;
            tmp1_3[i] = 0.0;

        }
        for (int i = 0; i <= n * n + 3*n ; i++)
        {
            tmp_0[i] = 0.0;
            tmp_1[i] = 0.0;
            tmp_2[i] = 0.0;
            tmp_3[i] = 0.0;
        }
        // generate p[0:n], e[0:n^2-1] into tmp
        for (int i = 0; i <= n; i++)
        {
            tmp_0[i] = p0[i];
            tmp_1[i] = p1[i];
            tmp_2[i] = p2[i];
            tmp_3[i] = p3[i];
        }
        for (int i = 0; i <= n * n - 1; i++)
        {
            tmp_0[n + i] = err0[i];
            tmp_1[n + i] = err1[i];
            tmp_2[n + i] = err2[i];
            tmp_3[n + i] = err3[i];
        }
        vecSum(tmp_0, tmp1_0, (n * n + n));
        vecSum(tmp_1, tmp1_1, (n * n + n));
        vecSum(tmp_2, tmp1_2, (n * n + n));
        vecSum(tmp_3, tmp1_3, (n * n + n));
        /* write  tmp1 into r_ext[n], e[0:n^2 +n-1]*/
        r_ext0[n] = tmp1_0[0];
        r_ext1[n] = tmp1_1[0];
        r_ext2[n] = tmp1_2[0];
        r_ext3[n] = tmp1_3[0];
        for (int i = 0; i <= n * n + n - 1; i++)
        {
            err0[i] = tmp1_0[i + 1];
            err1[i] = tmp1_1[i + 1];
            err2[i] = tmp1_2[i + 1];
            err3[i] = tmp1_3[i + 1];
        }
        // writes err[0:n^2 +n-1],e_tmp[0:n] into err[0:n^2 +2n-1]
        for (int i = 0; i <= n * n - 1; i++)
        {
            err0[n * n + n + i] = err0[i];
            err1[n * n + n + i] = err1[i];
            err2[n * n + n + i] = err2[i];
            err3[n * n + n + i] = err3[i];
        }
        for (int i = 0; i <= n; i++)
        {
            err0[i] = e_tmp0[i];
            err1[i] = e_tmp1[i];
            err2[i] = e_tmp2[i];
            err3[i] = e_tmp3[i];
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
      
        r_ext0[i] +=err0[ i];
        r_ext1[i] +=err1[ i];
        r_ext2[i] +=err2[ i];
        r_ext3[i] +=err3[ i];
    }

    renormalizationalgorithm(r_ext0, k + 1, r0, sizea);
    renormalizationalgorithm(r_ext1, k + 1, r1, sizea);
    renormalizationalgorithm(r_ext2, k + 1, r2, sizea);
    renormalizationalgorithm(r_ext3, k + 1, r3, sizea);
    return;
}