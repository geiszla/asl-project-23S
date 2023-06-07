#include <immintrin.h>
#include <cstring>

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

  int err_size = k * k;
  int ext_size = 2 * (k * k);

  double *err = (double *)alloca((err_size) * sizeof(double));
  double *r_ext = (double *)alloca((ext_size) * sizeof(double));

  memset(err, 0, (err_size) * sizeof(double));
  memset(r_ext, 0, (ext_size) * sizeof(double));

  int i = 0;
  int runner = 0;
  double r_ext_sum = 0;
  double err_sum = 0;

  double x = a[0];
  double y = b[0];
  double pi = x * y;

  err[0] = fma(x, y, -pi);
  r_ext[0] = pi;

  for (int n = 1; n <= (k - 1); n++)
  {
    int tmp_size = 2 * (n + 1);

    double *e_tmp = (double *)alloca((tmp_size) * sizeof(double));
    memset(e_tmp, 0, (tmp_size) * sizeof(double));

    double *p = (double *)alloca((tmp_size) * sizeof(double));
    memset(p, 0, (tmp_size) * sizeof(double));

    int runner1 = 0;
    int runner2 = 0;

    for (; runner2 <= n - 8; runner2 += 4)
    {
      // load b_vec first and then permute so it is faster
      __m256d b_vec = _mm256_loadu_pd(&b[n - runner2 - 3]);
      __m256d x_vec = _mm256_loadu_pd(&a[runner2]);
      __m256d y_vec = _mm256_permute4x64_pd(b_vec, 0b00011011);
      __m256d pi_vec = _mm256_mul_pd(x_vec, y_vec);
      __m256d e_tmp_vec = _mm256_fmsub_pd(x_vec, y_vec, pi_vec);

      _mm256_store_pd(&e_tmp[runner2], e_tmp_vec);
      _mm256_store_pd(&p[runner2], pi_vec);
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
    memset(tmp, 0, (n * n + n + 1) * sizeof(double));

    double *tmp1 = (double *)alloca((n * n + n + 1) * sizeof(double));
    memset(tmp1, 0, (n * n + n + 1) * sizeof(double));

    int run2 = 0;

    for (; run2 <= n - 16; run2 += 16)
    {
      __m256d x_vec = _mm256_loadu_pd(&p[run2]);
      _mm256_store_pd(&tmp[run2], x_vec);
      __m256d x_vec1 = _mm256_loadu_pd(&p[run2 + 4]);
      _mm256_store_pd(&tmp[run2 + 4], x_vec1);
      __m256d x_vec2 = _mm256_loadu_pd(&p[run2 + 8]);
      _mm256_store_pd(&tmp[run2 + 8], x_vec2);
      __m256d x_vec3 = _mm256_loadu_pd(&p[run2 + 12]);
      _mm256_store_pd(&tmp[run2 + 12], x_vec3);
    }

    for (; run2 <= n; run2++)
    {
      tmp[run2] = p[run2];
    }

    int run3 = 0;

    for (; run3 <= n * n - 5; run3 += 4)
    {
      __m256d x_vec = _mm256_loadu_pd(&err[run3]);
      _mm256_storeu_pd(&tmp[n + 1 + run3], x_vec);
      // tmp[n + 1 + run3] = err[run3];
    }

    for (; run3 <= n * n - 1; run3++)
    {
      tmp[n + 1 + run3] = err[run3];
    }

    vecSum5(tmp, tmp1, (n * n + n));

    r_ext[n] = tmp1[0];

    int run5 = 0;

    for (; run5 <= n * n + n - 5; run5 += 4)
    {
      // err[run5] = tmp1[run5 + 1];
      __m256d x_vec = _mm256_loadu_pd(&tmp1[run5 + 1]);
      _mm256_storeu_pd(&err[run5], x_vec);
    }

    for (; run5 <= n * n + n - 1; run5++)
    {
      err[run5] = tmp1[run5 + 1];
    }

    int run7 = n * n;

    for (; run7 <= ((n * n) + 2 * n - 5); run7++)
    {
      // err[run7] = e_tmp[run7 - n * n];
      __m256d x_vec = _mm256_loadu_pd(&e_tmp[run7 - n * n]);
      _mm256_storeu_pd(&err[run7], x_vec);
    }

    for (; run7 <= ((n * n) + 2 * n - 1); run7++)
    {
      err[run7] = e_tmp[run7 - n * n];
    }

    r_ext_sum += a[n] * b[k - n];
    err_sum += err[k - 1];
  }

  /////////////////////////////////////////////////

  r_ext[k] += r_ext_sum;
  r_ext[k - 1] += err_sum;

  int run8 = k * k;

  for (; run8 <= k * k - 5; run8 += 4)
  {
    // r_ext[run8] += err[run8];
    __m256d x_vec = _mm256_loadu_pd(&r_ext[run8]);
    __m256d y_vec = _mm256_loadu_pd(&err[run8]);
    __m256d z_vec = _mm256_add_pd(x_vec, y_vec);
    _mm256_storeu_pd(&r_ext[run8], z_vec);
  }

  for (; run8 <= k * k - 1; run8++)
  {
    r_ext[run8] += err[run8];
  }

  renormalization4(r_ext, k + 1, r, k);

  return;
}
