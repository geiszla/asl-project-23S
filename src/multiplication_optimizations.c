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

void multiplication3(double *a, double *b, double *r, const int sizea, const int sizeb,
                     const int sizer)
{
  int k = sizea;

  int err_size = k * k;
  double *err = (double *)alloca(err_size * sizeof(double));

  int i = 0;
  __m256d err_vec = _mm256_setzero_pd();
  for (; i < err_size - 32; i += 32)
  {
    _mm256_storeu_pd(&err[i], err_vec);
    _mm256_storeu_pd(&err[i + 4], err_vec);
    _mm256_storeu_pd(&err[i + 8], err_vec);
    _mm256_storeu_pd(&err[i + 12], err_vec);
    _mm256_storeu_pd(&err[i + 16], err_vec);
    _mm256_storeu_pd(&err[i + 20], err_vec);
    _mm256_storeu_pd(&err[i + 24], err_vec);
    _mm256_storeu_pd(&err[i + 28], err_vec);
  }
  for (; i < err_size; i += 1)
  {
    err[i] = 0;
  }

  int ext_size = 2 * (k * k);
  double *r_ext = (double *)alloca(ext_size * sizeof(double));

  int runner = 0;
  for (; runner < ext_size - 32; runner += 32)
  {
    _mm256_storeu_pd(&r_ext[runner], err_vec);
    _mm256_storeu_pd(&r_ext[runner + 4], err_vec);
    _mm256_storeu_pd(&r_ext[runner + 8], err_vec);
    _mm256_storeu_pd(&r_ext[runner + 12], err_vec);
    _mm256_storeu_pd(&r_ext[runner + 16], err_vec);
    _mm256_storeu_pd(&r_ext[runner + 20], err_vec);
    _mm256_storeu_pd(&r_ext[runner + 24], err_vec);
    _mm256_storeu_pd(&r_ext[runner + 28], err_vec);
  }
  for (; runner < ext_size; runner++)
  {
    r_ext[runner] = 0;
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

    int runner1 = 0;
    for (; runner1 < tmp_size - 32; runner1 += 32)
    {
      _mm256_storeu_pd(&e_tmp[runner1], err_vec);
      _mm256_storeu_pd(&p[runner1], err_vec);
      _mm256_storeu_pd(&e_tmp[runner1 + 4], err_vec);
      _mm256_storeu_pd(&p[runner1 + 4], err_vec);
      _mm256_storeu_pd(&e_tmp[runner1 + 8], err_vec);
      _mm256_storeu_pd(&p[runner1 + 8], err_vec);
      _mm256_storeu_pd(&e_tmp[runner1 + 12], err_vec);
      _mm256_storeu_pd(&p[runner1 + 12], err_vec);
      _mm256_storeu_pd(&e_tmp[runner1 + 16], err_vec);
      _mm256_storeu_pd(&p[runner1 + 16], err_vec);
      _mm256_storeu_pd(&e_tmp[runner1 + 20], err_vec);
      _mm256_storeu_pd(&p[runner1 + 20], err_vec);
      _mm256_storeu_pd(&e_tmp[runner1 + 24], err_vec);
      _mm256_storeu_pd(&p[runner1 + 24], err_vec);
      _mm256_storeu_pd(&e_tmp[runner1 + 28], err_vec);
      _mm256_storeu_pd(&p[runner1 + 28], err_vec);
    }
    for (; runner1 < tmp_size; runner1++)
    {
      e_tmp[runner1] = 0;
      p[runner1] = 0;
    }

    int runner2 = 0;
    for (; runner2 <= n - 12; runner2 += 4)
    {
      __m256d x_vec = _mm256_loadu_pd(&a[runner2]);
      __m256d b_vec = _mm256_loadu_pd(&b[n - runner2 - 3]);
      __m256d y_vec = _mm256_permute4x64_pd(b_vec, 0b00011011);
      __m256d pi_vec = _mm256_mul_pd(x_vec, y_vec);
      __m256d e_tmp_vec = _mm256_fmsub_pd(x_vec, y_vec, pi_vec);
      _mm256_storeu_pd(&e_tmp[runner2], e_tmp_vec);
      _mm256_storeu_pd(&p[runner2], pi_vec);
    }

    for (; runner2 <= n; runner2++)
    {
      x = a[runner2];
      y = b[n - runner2];
      pi = x * y;
      e_tmp[runner2] = fma(x, y, -pi);
      p[runner2] = pi;
    }

    double *tmp = (double *)alloca((n * n + n + 1) * sizeof(double));
    double *tmp1 = (double *)alloca((n * n + n + 1) * sizeof(double));

    int runner3 = 0;
    for (; runner3 <= n; runner3++)
    {
      tmp[runner3] = p[runner3];
      tmp[n + 1 + runner3] = err[runner3];
    }

    for (; runner3 <= n * n - 1; runner3++)
    {
      tmp[n + 1 + runner3] = err[runner3];
    }

    vecSum5(tmp, tmp1, (n * n + n));

    r_ext[n] = tmp1[0];

    for (int i = 0; i <= n * n + n - 1; i++)
    {
      err[i] = tmp1[i + 1];
    }

    for (int i = n * n; i <= ((n * n) + 2 * n - 1); i++)
    {
      err[i] = e_tmp[i - n * n];
    }
  }

  double r_ext_sum[4];
  double err_sum[4];
  r_ext_sum[0] = 0;
  err_sum[0] = 0;
  int runner4 = 0;

  __m256d r_ext_sum_vec = _mm256_setzero_pd();
  for (; runner4 <= k - 17; runner4 += 16)
  {
    __m256d x_vec = _mm256_load_pd(&a[runner4]);
    __m256d b_vec = _mm256_load_pd(&b[k - runner4 - 3]);

    __m256d x_vec1 = _mm256_load_pd(&a[runner4 + 4]);
    __m256d b_vec1 = _mm256_load_pd(&b[k - runner4 - 3 - 4]);
    __m256d x_vec2 = _mm256_load_pd(&a[runner4 + 8]);
    __m256d b_vec2 = _mm256_load_pd(&b[k - runner4 - 3 - 8]);
    __m256d x_vec3 = _mm256_load_pd(&a[runner4 + 12]);
    __m256d b_vec3 = _mm256_load_pd(&b[k - runner4 - 3 - 12]);

    __m256d y_vec = _mm256_permute4x64_pd(b_vec, 0b00011011);
    __m256d r = _mm256_mul_pd(x_vec, y_vec);
    r_ext_sum_vec = _mm256_add_pd(r_ext_sum_vec, r);

    __m256d y_vec1 = _mm256_permute4x64_pd(b_vec1, 0b00011011);
    __m256d r1 = _mm256_mul_pd(x_vec1, y_vec1);
    r_ext_sum_vec = _mm256_add_pd(r_ext_sum_vec, r1);

    __m256d y_vec2 = _mm256_permute4x64_pd(b_vec2, 0b00011011);
    __m256d r2 = _mm256_mul_pd(x_vec2, y_vec2);
    r_ext_sum_vec = _mm256_add_pd(r_ext_sum_vec, r2);

    __m256d y_vec3 = _mm256_permute4x64_pd(b_vec3, 0b00011011);
    __m256d r3 = _mm256_mul_pd(x_vec3, y_vec3);
    r_ext_sum_vec = _mm256_add_pd(r_ext_sum_vec, r3);
  }

  r_ext_sum_vec = _mm256_hadd_pd(r_ext_sum_vec, r_ext_sum_vec);
  r_ext_sum_vec = _mm256_hadd_pd(r_ext_sum_vec, r_ext_sum_vec);

  _mm256_storeu_pd(&r_ext_sum[0], r_ext_sum_vec);
  err_sum[0] += k * err[k - 1];
  for (; runner4 <= k - 1; runner4++)
  {
    r_ext_sum[0] += a[runner4] * b[k - runner4];
  }

  r_ext[k] += r_ext_sum[0];
  r_ext[k - 1] += err_sum[0];

  renormalization4(r_ext, k + 1, r, k);

  return;
}
