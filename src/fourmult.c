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

inline void renormalization_fast2(double x[], int size_of_x, double f[], int m)
{
  double *temp_array = (double *)alloca((3 * size_of_x - 3) * sizeof(double));

  if (m >= 20)
  {
    memset(temp_array, 0, (3 * size_of_x - 3) * sizeof(double));
  }

  double temp = x[size_of_x - 1];
  temp_array[size_of_x - 1] = temp;

  if (size_of_x < 11)
  {
    double st;
    double s3;
    double s2;
    double s1;

    int i;

    for (i = size_of_x - 2; i >= 3; i -= 4)
    {
      st = temp;
      s3 = x[i] + temp;
      s2 = x[i - 1] + s3;
      s1 = x[i - 2] + s2;
      temp = x[i - 3] + s1;

      __m256d s = _mm256_set_pd(s3, s2, s1, temp);
      __m256d b = _mm256_set_pd(st, s3, s2, s1);

      __m256d a = _mm256_loadu_pd(&x[i - 3]);

      __m256d e = _mm256_sub_pd(b, _mm256_sub_pd(s, a));

      _mm256_storeu_pd(&temp_array[i - 2], e);
    }

    for (; i >= 0; i--)
    {
      double a = x[i];
      double b = temp;

      temp = a + b;
      double e = b - (temp - a);

      temp_array[i + 1] = e;
    }
  }
  else if (size_of_x < 40)
  {
    for (int i = size_of_x - 2; i >= 0; i--)
    {
      double a = x[i];
      double b = temp;

      temp = a + temp;
      double e = b - (temp - a);

      temp_array[i + 1] = e;
    }
  }
  else
  {
    double s7;
    double s6;
    double s5;
    double s4;

    double s3;
    double s2;
    double s1;

    temp_array[size_of_x - 1] = temp;

    int i;

    for (i = size_of_x - 2; i >= 7; i -= 8)
    {
      s7 = x[i] + temp;
      temp_array[i] = s7;

      s6 = x[i - 1] + s7;
      temp_array[i - 1] = s6;

      s5 = x[i - 2] + s6;
      temp_array[i - 2] = s5;

      s4 = x[i - 3] + s5;
      temp_array[i - 3] = s4;

      s3 = x[i - 4] + s4;
      temp_array[i - 4] = s3;

      s2 = x[i - 5] + s3;
      temp_array[i - 5] = s2;

      s1 = x[i - 6] + s2;
      temp_array[i - 6] = s1;

      temp = x[i - 7] + s1;
      temp_array[i - 7] = temp;
    }

    for (i = size_of_x - 2; i >= 7; i -= 8)
    {
      int start_index = i - 3;
      int next_start_index = start_index - 4;

      __m256d a0 = _mm256_loadu_pd(&x[start_index]);
      __m256d a1 = _mm256_loadu_pd(&x[next_start_index]);

      __m256d b0 = _mm256_loadu_pd(&temp_array[start_index + 1]);
      __m256d b1 = _mm256_loadu_pd(&temp_array[next_start_index + 1]);

      __m256d s0 = _mm256_loadu_pd(&temp_array[start_index]);
      __m256d s1 = _mm256_loadu_pd(&temp_array[next_start_index]);

      __m256d e0 = _mm256_sub_pd(b0, _mm256_sub_pd(s0, a0));
      __m256d e1 = _mm256_sub_pd(b1, _mm256_sub_pd(s1, a1));

      _mm256_storeu_pd(&temp_array[start_index + 1], e0);
      _mm256_storeu_pd(&temp_array[next_start_index + 1], e1);
    }

    for (; i >= 0; i--)
    {
      double a = x[i];
      double b = temp;

      temp = a + b;
      double e = b - (temp - a);

      temp_array[i + 1] = e;
    }
  }

  temp_array[0] = temp;

  // vecSumErrBranch
  int j = 0;

  double last_error = temp_array[1];
  if (last_error != 0)
  {
    temp = last_error;
    j++;
  }

  double j_limit = m < size_of_x ? m : m - 1;

  for (int i = 2; i <= size_of_x - 1; i++)
  {
    double b = temp_array[i];

    double s = temp + b;
    temp = b - (s - temp);

    temp_array[j] = s;

    if (temp == 0)
    {
      temp = s;
    }
    else if (j == j_limit)
    {
      break;
    }
    else
    {
      j++;
    }
  }

  if (temp != 0 && j < j_limit)
  {
    temp_array[j] = temp;
    j++;
  }

  for (; j <= j_limit; j++)
  {
    temp_array[j] = 0;
  }

  int vec_sum_limit = m < size_of_x ? m - 1 : size_of_x - 2;

  int i = 0;

  if (m >= 20)
  {
    double *s_array = (double *)alloca(4 * sizeof(double));

    for (; i < (vec_sum_limit / 3) - 3; i += 4)
    {
      double e0 = temp_array[i];
      double b0 = temp_array[i + 1];

      double s0 = e0 + b0;
      e0 = b0 - (s0 - e0);

      f[i] = s0;

      b0 = temp_array[i + 2];

      double s1 = e0 + b0;
      e0 = b0 - (s1 - e0);

      b0 = temp_array[i + 3];

      double s2 = e0 + b0;
      e0 = b0 - (s2 - e0);

      __m256d e = _mm256_set_pd(0, 0, s1, e0);
      __m256d b = _mm256_set_pd(0, 0, s2, temp_array[i + 4]);

      __m256d s = _mm256_add_pd(e, b);
      e = _mm256_sub_pd(b, _mm256_sub_pd(s, e));

      _mm256_storeu_pd(s_array, s);
      f[i + 1] = s_array[1];

      b = _mm256_unpacklo_pd(_mm256_set_pd(0, 0, 0, temp_array[i + 5]), s);

      s = _mm256_add_pd(e, b);
      e = _mm256_sub_pd(b, _mm256_sub_pd(s, e));

      b = _mm256_unpacklo_pd(_mm256_set_pd(0, 0, 0, temp_array[i + 6]), s);

      __m256d st = _mm256_permute4x64_pd(s, 0b11011111);
      s = _mm256_add_pd(e, b);

      e = _mm256_sub_pd(b, _mm256_sub_pd(s, e));
      e = _mm256_add_pd(e, st);

      st = _mm256_permute4x64_pd(s, 0b10010011);
      b = _mm256_add_pd(st, _mm256_loadu_pd(&temp_array[i + 7]));

      s = _mm256_add_pd(e, b);
      e = _mm256_sub_pd(b, _mm256_sub_pd(s, e));

      _mm256_storeu_pd(s_array, s);
      f[i + 2] = s_array[2];

      st = _mm256_permute4x64_pd(s, 0b10010011);
      b = _mm256_add_pd(st, _mm256_loadu_pd(&temp_array[i + 8]));

      s = _mm256_add_pd(e, b);
      e = _mm256_sub_pd(b, _mm256_sub_pd(s, e));

      st = _mm256_permute4x64_pd(s, 0b10010011);
      b = _mm256_blend_pd(st, _mm256_loadu_pd(&temp_array[i + 9]), 0b0001);

      st = _mm256_permute4x64_pd(s, 0b10111111);
      s = _mm256_add_pd(e, b);
      e = _mm256_sub_pd(b, _mm256_sub_pd(s, e));
      e = _mm256_blend_pd(e, st, 0b1000);

      st = _mm256_permute4x64_pd(s, 0b10010011);
      b = _mm256_blend_pd(st, _mm256_loadu_pd(&temp_array[i + 10]), 0b0001);

      s = _mm256_add_pd(e, b);
      e = _mm256_sub_pd(b, _mm256_sub_pd(s, e));

      _mm256_storeu_pd(s_array, s);
      f[i + 3] = s_array[3];

      int offset = i + 11;

      for (int j = i + 4; j < vec_sum_limit; j++)
      {
        st = _mm256_permute4x64_pd(s, 0b10010011);
        b = _mm256_blend_pd(st, _mm256_loadu_pd(&temp_array[offset]), 0b0001);

        s = _mm256_add_pd(e, b);
        e = _mm256_sub_pd(b, _mm256_sub_pd(s, e));

        st = _mm256_permute4x64_pd(s, 0b00000011);
        _mm256_storeu_pd(&temp_array[j], st);

        offset++;
      }
    }
  }

  for (; i < vec_sum_limit; i++)
  {
    // vecSumErr
    double e = temp_array[i];

    for (int j = i; j < vec_sum_limit; j++)
    {
      double b = temp_array[j + 1];

      double s = e + b;
      e = b - (s - e);

      temp_array[j] = s;
    }

    temp_array[vec_sum_limit] = e;

    f[i] = temp_array[i];
  }

  f[vec_sum_limit] = temp_array[vec_sum_limit];

  if (m == size_of_x)
  {
    f[size_of_x - 1] = temp_array[size_of_x - 1];
  }
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
        i = length - 2;
       
        
        for (; i >= 0; i--)
        {
            double s_tmp0, e_tmp0;
            double s_tmp1, e_tmp1;
            double s_tmp2, e_tmp2;
            double s_tmp3, e_tmp3;
            s0[i] = tmp_0[i] + s0[i + 1];
            s1[i] = tmp_1[i] + s1[i + 1];
            s2[i] = tmp_2[i] + s2[i + 1];
            s3[i] = tmp_3[i] + s3[i + 1];
            double t0 = s0[i] - s0[i + 1];
            double t1 = s1[i] - s1[i + 1];
            double t2 = s2[i] - s2[i + 1];
            double t3 = s3[i] - s3[i + 1];

            e_tmp0 = (tmp_0[i] - t0) + (s0[i + 1] - (s0[i] - t0));
            e_tmp1 = (tmp_1[i] - t1) + (s1[i + 1] - (s1[i] - t1));
            e_tmp2 = (tmp_2[i] - t2) + (s2[i + 1] - (s2[i] - t2));
            e_tmp3 = (tmp_3[i] - t3) + (s3[i + 1] - (s3[i] - t3));

           
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

    renormalization_fast2(r_ext0, k + 1, r0, sizea);
    renormalization_fast2(r_ext1, k + 1, r1, sizea);
    renormalization_fast2(r_ext2, k + 1, r2, sizea);
    renormalization_fast2(r_ext3, k + 1, r3, sizea);

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
