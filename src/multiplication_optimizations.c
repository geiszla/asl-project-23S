#include <immintrin.h>

void multiplication2(double *a, double *b, double *r, const int sizea, const int sizeb,
                     const int sizer)
{
  int k = sizea;

  int err_size = k * k;
  double *err = (double *)alloca(err_size * sizeof(double));

  for (int i = 0; i < err_size; i++)
  {
    err[i] = 0;
  }

  int ext_size = 2 * (k * k);
  double *r_ext = (double *)alloca(ext_size * sizeof(double));

  for (int i = 0; i < ext_size; i++)
  {
    r_ext[i] = 0;
  }

  double x = a[0];
  double y = b[0];

  double pi = x * y;

  err[0] = fma(x, y, -pi);
  r_ext[0] = pi;

  for (int n = 1; n <= (k - 1); n++)
  {
    int tmp_size = 2 * (n + 1);

    double *e_tmp = (double *)alloca(tmp_size * sizeof(double));
    double *p = (double *)alloca(tmp_size * sizeof(double));

    for (int i = 0; i < tmp_size; i++)
    {
      e_tmp[i] = 0;
      p[i] = 0;
    }

    for (int i = 0; i <= n; i++)
    {
      x = a[i];
      y = b[n - i];

      pi = x * y;

      e_tmp[i] = fma(x, y, -pi);
      p[i] = pi;
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

    for (int i = 0; i <= n * n - 1; i++)
    {
      tmp[n + 1 + i] = err[i];
    }

    vecSum5(tmp, tmp1, (n * n + n));

    /* write  tmp1 into r_ext[n], e[0:n^2 +n-1]*/
    r_ext[n] = tmp1[0];

    for (int i = 0; i <= n * n + n - 1; i++)
    {
      err[i] = tmp1[i + 1];
    }

    // writes err[0:n^2 +n-1],e_tmp[0:n] into err[0:n^2 +2n-1]
    for (int i = n * n; i <= ((n * n) + 2 * n - 1); i++)
    {
      err[i] = e_tmp[i - n * n];
    }
  }

  for (int i = 1; i <= k - 1; i++)
  {
    r_ext[k] += a[i] * b[k - i];
  }

  for (int i = 0; i <= k * k - 1; i++)
  {
    r_ext[i] += err[i];
  }

  renormalization4(r_ext, k + 1, r, k);

  return;
}

void multiplication3(double *a, double *b, double *r, const int sizea, const int sizeb, const int sizer)
{

  int k = sizea;
  double *err = (double *)calloc((sizea * sizea + 2 * sizea + 1), sizeof(double));
  double *r_ext = (double *)calloc((sizea * sizea), sizeof(double));

  double pi = a[0] * b[0];
  double e = fma(a[0], b[0], -pi);
  r_ext[0] = pi;
  err[0] = e;

  for (int n = 1; n <= (k - 1); n++)
  {
    double *e_tmp = (double *)calloc((n + 1), sizeof(double));
    double *p = (double *)calloc((n + 1), sizeof(double));
    int i = 0;
    for (; i <= n - 32; i += 32)
    {
      __m256d a1 = _mm256_loadu_pd(&a[i]);
      __m256d b1 = _mm256_loadu_pd(&b[n - i - 3]);
      b1 = _mm256_permute4x64_pd(b1, 0b00011011);
      __m256d pi1 = _mm256_mul_pd(a1, b1);
      __m256d e1 = _mm256_fmsub_pd(a1, b1, pi1);
      _mm256_storeu_pd(&p[i], pi1);
      _mm256_storeu_pd(&e_tmp[i], e1);

      __m256d a2 = _mm256_loadu_pd(&a[i + 4]);
      __m256d b2 = _mm256_loadu_pd(&b[n - i - 7]);
      b2 = _mm256_permute4x64_pd(b2, 0b00011011);
      __m256d pi2 = _mm256_mul_pd(a2, b2);
      __m256d e2 = _mm256_fmsub_pd(a2, b2, pi2);
      _mm256_storeu_pd(&p[i + 4], pi2);
      _mm256_storeu_pd(&e_tmp[i + 4], e2);

      __m256d a3 = _mm256_loadu_pd(&a[i + 8]);
      __m256d b3 = _mm256_loadu_pd(&b[n - i - 11]);
      b3 = _mm256_permute4x64_pd(b3, 0b00011011);
      __m256d pi3 = _mm256_mul_pd(a3, b3);
      __m256d e3 = _mm256_fmsub_pd(a3, b3, pi3);
      _mm256_storeu_pd(&p[i + 8], pi3);
      _mm256_storeu_pd(&e_tmp[i + 8], e3);

      __m256d a4 = _mm256_loadu_pd(&a[i + 12]);
      __m256d b4 = _mm256_loadu_pd(&b[n - i - 15]);
      b4 = _mm256_permute4x64_pd(b4, 0b00011011);
      __m256d pi4 = _mm256_mul_pd(a4, b4);
      __m256d e4 = _mm256_fmsub_pd(a4, b4, pi4);
      _mm256_storeu_pd(&p[i + 12], pi4);
      _mm256_storeu_pd(&e_tmp[i + 12], e4);

      __m256d a5 = _mm256_loadu_pd(&a[i + 16]);
      __m256d b5 = _mm256_loadu_pd(&b[n - i - 19]);
      b5 = _mm256_permute4x64_pd(b5, 0b00011011);
      __m256d pi5 = _mm256_mul_pd(a5, b5);
      __m256d e5 = _mm256_fmsub_pd(a5, b5, pi5);
      _mm256_storeu_pd(&p[i + 16], pi5);
      _mm256_storeu_pd(&e_tmp[i + 16], e5);

      __m256d a6 = _mm256_loadu_pd(&a[i + 20]);
      __m256d b6 = _mm256_loadu_pd(&b[n - i - 23]);
      b6 = _mm256_permute4x64_pd(b6, 0b00011011);
      __m256d pi6 = _mm256_mul_pd(a6, b6);
      __m256d e6 = _mm256_fmsub_pd(a6, b6, pi6);
      _mm256_storeu_pd(&p[i + 20], pi6);
      _mm256_storeu_pd(&e_tmp[i + 20], e6);

      __m256d a7 = _mm256_loadu_pd(&a[i + 24]);
      __m256d b7 = _mm256_loadu_pd(&b[n - i - 27]);
      b7 = _mm256_permute4x64_pd(b7, 0b00011011);
      __m256d pi7 = _mm256_mul_pd(a7, b7);
      __m256d e7 = _mm256_fmsub_pd(a7, b7, pi7);
      _mm256_storeu_pd(&p[i + 24], pi7);
      _mm256_storeu_pd(&e_tmp[i + 24], e7);

      __m256d a8 = _mm256_loadu_pd(&a[i + 28]);
      __m256d b8 = _mm256_loadu_pd(&b[n - i - 31]);
      b8 = _mm256_permute4x64_pd(b8, 0b00011011);
      __m256d pi8 = _mm256_mul_pd(a8, b8);
      __m256d e8 = _mm256_fmsub_pd(a8, b8, pi8);
      _mm256_storeu_pd(&p[i + 28], pi8);
      _mm256_storeu_pd(&e_tmp[i + 28], e8);
    }
    for (; i <= n; i++)
    {
      // twoMultFMA(a[i], b[n - i], &(e_tmp[i]), &(e_tmp[i]));
      double pi = a[i] * b[n - i];
      double e = fma(a[i], b[n - i], -pi);
      p[i] = pi;
      e_tmp[i] = e;
    }
    double *tmp = &err[-n - 1];
    double *tmp1 = (double *)calloc((n * n + n + 1), sizeof(double));

    //-------------------------------- vecsum inline
    int length = (n * n + n);
    double *s = (double *)alloca(length * sizeof(double));
    s[length - 1] = err[n * n - 1];
    for (int i = length - 2; i >= n + 1; i -= 1)
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
    for (; i <= (n * n + n - 1) - 32; i += 32)
    {
      __m256d tmp1_vec = _mm256_loadu_pd(&tmp1[i + 1]);
      _mm256_storeu_pd(&err[i], tmp1_vec);
      __m256d tmp1_vec1 = _mm256_loadu_pd(&tmp1[i + 1 + 4]);
      _mm256_storeu_pd(&err[i + 4], tmp1_vec1);
      __m256d tmp1_vec2 = _mm256_loadu_pd(&tmp1[i + 1 + 8]);
      _mm256_storeu_pd(&err[i + 8], tmp1_vec2);
      __m256d tmp1_vec3 = _mm256_loadu_pd(&tmp1[i + 1 + 12]);
      _mm256_storeu_pd(&err[i + 12], tmp1_vec3);
      __m256d tmp1_vec4 = _mm256_loadu_pd(&tmp1[i + 1 + 16]);
      _mm256_storeu_pd(&err[i + 16], tmp1_vec4);
      __m256d tmp1_vec5 = _mm256_loadu_pd(&tmp1[i + 1 + 20]);
      _mm256_storeu_pd(&err[i + 20], tmp1_vec5);
      __m256d tmp1_vec6 = _mm256_loadu_pd(&tmp1[i + 1 + 24]);
      _mm256_storeu_pd(&err[i + 24], tmp1_vec6);
      __m256d tmp1_vec7 = _mm256_loadu_pd(&tmp1[i + 1 + 28]);
      _mm256_storeu_pd(&err[i + 28], tmp1_vec7);

      // err[i] = tmp1[i + 1];
    }
    for (; i <= (n * n + n - 1); i++)
    {
      err[i] = tmp1[i + 1];
    }
    //---------------------------------------------

    // writes err[0:n^2 +n-1],e_tmp[0:n] into err[0:n^2 +2n-1]
    i = (n * n + n);
    for (; i <= (n * n + 2 * n - 1) - 32; i += 32)
    {
      __m256d tmp1_vec = _mm256_loadu_pd(&tmp1[i - (n * n + n)]);
      _mm256_storeu_pd(&err[i], tmp1_vec);
      __m256d tmp1_vec1 = _mm256_loadu_pd(&tmp1[i + 4 - (n * n + n)]);
      _mm256_storeu_pd(&err[i + 4], tmp1_vec1);
      __m256d tmp1_vec2 = _mm256_loadu_pd(&tmp1[i + 8 - (n * n + n)]);
      _mm256_storeu_pd(&err[i + 8], tmp1_vec2);
      __m256d tmp1_vec3 = _mm256_loadu_pd(&tmp1[i + 12 - (n * n + n)]);
      _mm256_storeu_pd(&err[i + 12], tmp1_vec3);
      __m256d tmp1_vec4 = _mm256_loadu_pd(&tmp1[i + 16 - (n * n + n)]);
      _mm256_storeu_pd(&err[i + 16], tmp1_vec4);
      __m256d tmp1_vec5 = _mm256_loadu_pd(&tmp1[i + 20 - (n * n + n)]);
      _mm256_storeu_pd(&err[i + 20], tmp1_vec5);
      __m256d tmp1_vec6 = _mm256_loadu_pd(&tmp1[i + 24 - (n * n + n)]);
      _mm256_storeu_pd(&err[i + 24], tmp1_vec6);
      __m256d tmp1_vec7 = _mm256_loadu_pd(&tmp1[i + 28 - (n * n + n)]);
      _mm256_storeu_pd(&err[i + 28], tmp1_vec7);
    }
    for (; i <= (n * n + 2 * n - 1); i++)
    {
      err[i] = e_tmp[i - (n * n + n)];
    }
  }
  int i = 1;
  for (; i <= (k - 1) - 32; i += 32)
  {
    __m256d a = _mm256_loadu_pd(&a[i]);
    __m256d b = _mm256_loadu_pd(&b[k - i - 3]);
    __m256d r_ext = _mm256_loadu_pd(&r_ext[k]);
    __m256d tmp1_vec = _mm256_permute4x64_pd(b, 0b00011011);
    __m256d mul = _mm256_mul_pd(a, tmp1_vec);
    r_ext = _mm256_add_pd(r_ext, mul);

    __m256d a1 = _mm256_loadu_pd(&a[i + 4]);
    __m256d b1 = _mm256_loadu_pd(&b[k - i - 3 - 4]);
    __m256d tmp1_vec1 = _mm256_permute4x64_pd(b1, 0b00011011);
    __m256d mul1 = _mm256_mul_pd(a1, tmp1_vec1);
    r_ext = _mm256_add_pd(r_ext, mul1);

    __m256d a2 = _mm256_loadu_pd(&a[i + 8]);
    __m256d b2 = _mm256_loadu_pd(&b[k - i - 3 - 8]);
    __m256d tmp1_vec2 = _mm256_permute4x64_pd(b2, 0b00011011);
    __m256d mul2 = _mm256_mul_pd(a2, tmp1_vec2);
    r_ext = _mm256_add_pd(r_ext, mul2);

    __m256d a3 = _mm256_loadu_pd(&a[i + 12]);
    __m256d b3 = _mm256_loadu_pd(&b[k - i - 3 - 12]);
    __m256d tmp1_vec3 = _mm256_permute4x64_pd(b3, 0b00011011);
    __m256d mul3 = _mm256_mul_pd(a3, tmp1_vec3);
    r_ext = _mm256_add_pd(r_ext, mul3);

    __m256d a4 = _mm256_loadu_pd(&a[i + 16]);
    __m256d b4 = _mm256_loadu_pd(&b[k - i - 3 - 16]);
    __m256d tmp1_vec4 = _mm256_permute4x64_pd(b4, 0b00011011);
    __m256d mul4 = _mm256_mul_pd(a4, tmp1_vec4);
    r_ext = _mm256_add_pd(r_ext, mul4);

    __m256d a5 = _mm256_loadu_pd(&a[i + 20]);
    __m256d b5 = _mm256_loadu_pd(&b[k - i - 3 - 20]);
    __m256d tmp1_vec5 = _mm256_permute4x64_pd(b5, 0b00011011);
    __m256d mul5 = _mm256_mul_pd(a5, tmp1_vec5);
    r_ext = _mm256_add_pd(r_ext, mul5);

    __m256d a6 = _mm256_loadu_pd(&a[i + 24]);
    __m256d b6 = _mm256_loadu_pd(&b[k - i - 3 - 24]);
    __m256d tmp1_vec6 = _mm256_permute4x64_pd(b6, 0b00011011);
    __m256d mul6 = _mm256_mul_pd(a6, tmp1_vec6);
    r_ext = _mm256_add_pd(r_ext, mul6);

    __m256d a7 = _mm256_loadu_pd(&a[i + 28]);
    __m256d b7 = _mm256_loadu_pd(&b[k - i - 3 - 28]);
    __m256d tmp1_vec7 = _mm256_permute4x64_pd(b7, 0b00011011);
    __m256d mul7 = _mm256_mul_pd(a7, tmp1_vec7);
    r_ext = _mm256_add_pd(r_ext, mul7);

    _mm256_storeu_pd(&r_ext[k], r_ext);
  }
  for (; i <= (k - 1); i++)
  {

    r_ext[k] += a[i] * b[k - i];
  }

  i = 0;
  for (; i <= (k * k - 1) - 32; i += 32)
  {
    __m256d tmp1_vec = _mm256_loadu_pd(&err[i]);
    __m256d tmp2_vec = _mm256_loadu_pd(&r_ext[i]);
    _mm256_storeu_pd(&r_ext[i], _mm256_add_pd(tmp1_vec, tmp2_vec));
    __m256d tmp1_vec1 = _mm256_loadu_pd(&err[i + 4]);
    __m256d tmp2_vec1 = _mm256_loadu_pd(&r_ext[i + 4]);
    _mm256_storeu_pd(&r_ext[i + 4], _mm256_add_pd(tmp1_vec1, tmp2_vec1));
    __m256d tmp1_vec2 = _mm256_loadu_pd(&err[i + 8]);
    __m256d tmp2_vec2 = _mm256_loadu_pd(&r_ext[i + 8]);
    _mm256_storeu_pd(&r_ext[i + 8], _mm256_add_pd(tmp1_vec2, tmp2_vec2));
    __m256d tmp1_vec3 = _mm256_loadu_pd(&err[i + 12]);
    __m256d tmp2_vec3 = _mm256_loadu_pd(&r_ext[i + 12]);
    _mm256_storeu_pd(&r_ext[i + 12], _mm256_add_pd(tmp1_vec3, tmp2_vec3));
    __m256d tmp1_vec4 = _mm256_loadu_pd(&err[i + 16]);
    __m256d tmp2_vec4 = _mm256_loadu_pd(&r_ext[i + 16]);
    _mm256_storeu_pd(&r_ext[i + 16], _mm256_add_pd(tmp1_vec4, tmp2_vec4));
    __m256d tmp1_vec5 = _mm256_loadu_pd(&err[i + 20]);
    __m256d tmp2_vec5 = _mm256_loadu_pd(&r_ext[i + 20]);
    _mm256_storeu_pd(&r_ext[i + 20], _mm256_add_pd(tmp1_vec5, tmp2_vec5));
    __m256d tmp1_vec6 = _mm256_loadu_pd(&err[i + 24]);
    __m256d tmp2_vec6 = _mm256_loadu_pd(&r_ext[i + 24]);
    _mm256_storeu_pd(&r_ext[i + 24], _mm256_add_pd(tmp1_vec6, tmp2_vec6));
    __m256d tmp1_vec7 = _mm256_loadu_pd(&err[i + 28]);
    __m256d tmp2_vec7 = _mm256_loadu_pd(&r_ext[i + 28]);
    _mm256_storeu_pd(&r_ext[i + 28], _mm256_add_pd(tmp1_vec7, tmp2_vec7));
  }
  for (; i <= (k * k - 1); i++)
  {

    r_ext[i] += err[i];
  }

  renormalizationalgorithm(r_ext, k + 1, r, sizea);
  return;
}
