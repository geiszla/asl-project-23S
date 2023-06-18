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

void multiplication_fast(double *a, double *b, double *r, const int sizea, const int sizeb, const int sizer)
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

    int size_of_x = n * n + n;

    double temp = tmp[size_of_x - 1];
    tmp1[size_of_x - 1] = temp;

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
        s3 = tmp[i] + temp;
        s2 = tmp[i - 1] + s3;
        s1 = tmp[i - 2] + s2;
        temp = tmp[i - 3] + s1;

        __m256d s = _mm256_set_pd(s3, s2, s1, temp);
        __m256d b = _mm256_set_pd(st, s3, s2, s1);

        __m256d a = _mm256_loadu_pd(&tmp[i - 3]);

        __m256d e = _mm256_sub_pd(b, _mm256_sub_pd(s, a));

        _mm256_storeu_pd(&tmp1[i - 2], e);
      }

      for (; i >= 0; i--)
      {
        double a = tmp[i];
        double b = temp;

        temp = a + b;
        double e = b - (temp - a);

        tmp1[i + 1] = e;
      }
    }
    else if (size_of_x < 40)
    {
      for (int i = size_of_x - 2; i >= 0; i--)
      {
        double a = tmp[i];
        double b = temp;

        temp = a + temp;
        double e = b - (temp - a);

        tmp1[i + 1] = e;
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

      tmp1[size_of_x - 1] = temp;

      int i;

      for (i = size_of_x - 2; i >= 7; i -= 8)
      {
        s7 = tmp[i] + temp;
        tmp1[i] = s7;

        s6 = tmp[i - 1] + s7;
        tmp1[i - 1] = s6;

        s5 = tmp[i - 2] + s6;
        tmp1[i - 2] = s5;

        s4 = tmp[i - 3] + s5;
        tmp1[i - 3] = s4;

        s3 = tmp[i - 4] + s4;
        tmp1[i - 4] = s3;

        s2 = tmp[i - 5] + s3;
        tmp1[i - 5] = s2;

        s1 = tmp[i - 6] + s2;
        tmp1[i - 6] = s1;

        temp = tmp[i - 7] + s1;
        tmp1[i - 7] = temp;
      }

      for (i = size_of_x - 2; i >= 7; i -= 8)
      {
        int start_index = i - 3;
        int next_start_index = start_index - 4;

        __m256d a0 = _mm256_loadu_pd(&tmp[start_index]);
        __m256d a1 = _mm256_loadu_pd(&tmp[next_start_index]);

        __m256d b0 = _mm256_loadu_pd(&tmp1[start_index + 1]);
        __m256d b1 = _mm256_loadu_pd(&tmp1[next_start_index + 1]);

        __m256d s0 = _mm256_loadu_pd(&tmp1[start_index]);
        __m256d s1 = _mm256_loadu_pd(&tmp1[next_start_index]);

        __m256d e0 = _mm256_sub_pd(b0, _mm256_sub_pd(s0, a0));
        __m256d e1 = _mm256_sub_pd(b1, _mm256_sub_pd(s1, a1));

        _mm256_storeu_pd(&tmp1[start_index + 1], e0);
        _mm256_storeu_pd(&tmp1[next_start_index + 1], e1);
      }

      for (; i >= 0; i--)
      {
        double a = tmp[i];
        double b = temp;

        temp = a + b;
        double e = b - (temp - a);

        tmp1[i + 1] = e;
      }
    }

    tmp1[0] = temp;
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

  int size_of_x = k + 1;

  if (k >= 20)
  {
    memset(err, 0, (3 * size_of_x - 3) * sizeof(double));
  }

  double temp = r_ext[size_of_x - 1];
  err[size_of_x - 1] = temp;

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
      s3 = r_ext[i] + temp;
      s2 = r_ext[i - 1] + s3;
      s1 = r_ext[i - 2] + s2;
      temp = r_ext[i - 3] + s1;

      __m256d s = _mm256_set_pd(s3, s2, s1, temp);
      __m256d b = _mm256_set_pd(st, s3, s2, s1);

      __m256d a = _mm256_loadu_pd(&r_ext[i - 3]);

      __m256d e = _mm256_sub_pd(b, _mm256_sub_pd(s, a));

      _mm256_storeu_pd(&err[i - 2], e);
    }

    for (; i >= 0; i--)
    {
      double a = r_ext[i];
      double b = temp;

      temp = a + b;
      double e = b - (temp - a);

      err[i + 1] = e;
    }
  }
  else if (size_of_x < 40)
  {
    for (int i = size_of_x - 2; i >= 0; i--)
    {
      double a = r_ext[i];
      double b = temp;

      temp = a + temp;
      double e = b - (temp - a);

      err[i + 1] = e;
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

    err[size_of_x - 1] = temp;

    int i;

    for (i = size_of_x - 2; i >= 7; i -= 8)
    {
      s7 = r_ext[i] + temp;
      err[i] = s7;

      s6 = r_ext[i - 1] + s7;
      err[i - 1] = s6;

      s5 = r_ext[i - 2] + s6;
      err[i - 2] = s5;

      s4 = r_ext[i - 3] + s5;
      err[i - 3] = s4;

      s3 = r_ext[i - 4] + s4;
      err[i - 4] = s3;

      s2 = r_ext[i - 5] + s3;
      err[i - 5] = s2;

      s1 = r_ext[i - 6] + s2;
      err[i - 6] = s1;

      temp = r_ext[i - 7] + s1;
      err[i - 7] = temp;
    }

    for (i = size_of_x - 2; i >= 7; i -= 8)
    {
      int start_index = i - 3;
      int next_start_index = start_index - 4;

      __m256d a0 = _mm256_loadu_pd(&r_ext[start_index]);
      __m256d a1 = _mm256_loadu_pd(&r_ext[next_start_index]);

      __m256d b0 = _mm256_loadu_pd(&err[start_index + 1]);
      __m256d b1 = _mm256_loadu_pd(&err[next_start_index + 1]);

      __m256d s0 = _mm256_loadu_pd(&err[start_index]);
      __m256d s1 = _mm256_loadu_pd(&err[next_start_index]);

      __m256d e0 = _mm256_sub_pd(b0, _mm256_sub_pd(s0, a0));
      __m256d e1 = _mm256_sub_pd(b1, _mm256_sub_pd(s1, a1));

      _mm256_storeu_pd(&err[start_index + 1], e0);
      _mm256_storeu_pd(&err[next_start_index + 1], e1);
    }

    for (; i >= 0; i--)
    {
      double a = r_ext[i];
      double b = temp;

      temp = a + b;
      double e = b - (temp - a);

      err[i + 1] = e;
    }
  }

  err[0] = temp;

  // vecSumErrBranch
  int j = 0;

  double last_error = err[1];
  if (last_error != 0)
  {
    temp = last_error;
    j++;
  }

  double j_limit = k < size_of_x ? k : k - 1;

  for (int i = 2; i <= size_of_x - 1; i++)
  {
    double b = err[i];

    double s = temp + b;
    temp = b - (s - temp);

    err[j] = s;

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
    err[j] = temp;
    j++;
  }

  for (; j <= j_limit; j++)
  {
    err[j] = 0;
  }

  int vec_sum_limit = k < size_of_x ? k - 1 : size_of_x - 2;

  i = 0;

  if (k >= 20)
  {
    double *s_array = (double *)alloca(4 * sizeof(double));

    for (; i < (vec_sum_limit / 3) - 3; i += 4)
    {
      double e0 = err[i];
      double b0 = err[i + 1];

      double s0 = e0 + b0;
      e0 = b0 - (s0 - e0);

      r[i] = s0;

      b0 = err[i + 2];

      double s1 = e0 + b0;
      e0 = b0 - (s1 - e0);

      b0 = err[i + 3];

      double s2 = e0 + b0;
      e0 = b0 - (s2 - e0);

      __m256d e = _mm256_set_pd(0, 0, s1, e0);
      __m256d b = _mm256_set_pd(0, 0, s2, err[i + 4]);

      __m256d s = _mm256_add_pd(e, b);
      e = _mm256_sub_pd(b, _mm256_sub_pd(s, e));

      _mm256_storeu_pd(s_array, s);
      r[i + 1] = s_array[1];

      b = _mm256_unpacklo_pd(_mm256_set_pd(0, 0, 0, err[i + 5]), s);

      s = _mm256_add_pd(e, b);
      e = _mm256_sub_pd(b, _mm256_sub_pd(s, e));

      b = _mm256_unpacklo_pd(_mm256_set_pd(0, 0, 0, err[i + 6]), s);

      __m256d st = _mm256_permute4x64_pd(s, 0b11011111);
      s = _mm256_add_pd(e, b);

      e = _mm256_sub_pd(b, _mm256_sub_pd(s, e));
      e = _mm256_add_pd(e, st);

      st = _mm256_permute4x64_pd(s, 0b10010011);
      b = _mm256_add_pd(st, _mm256_loadu_pd(&err[i + 7]));

      s = _mm256_add_pd(e, b);
      e = _mm256_sub_pd(b, _mm256_sub_pd(s, e));

      _mm256_storeu_pd(s_array, s);
      r[i + 2] = s_array[2];

      st = _mm256_permute4x64_pd(s, 0b10010011);
      b = _mm256_add_pd(st, _mm256_loadu_pd(&err[i + 8]));

      s = _mm256_add_pd(e, b);
      e = _mm256_sub_pd(b, _mm256_sub_pd(s, e));

      st = _mm256_permute4x64_pd(s, 0b10010011);
      b = _mm256_blend_pd(st, _mm256_loadu_pd(&err[i + 9]), 0b0001);

      st = _mm256_permute4x64_pd(s, 0b10111111);
      s = _mm256_add_pd(e, b);
      e = _mm256_sub_pd(b, _mm256_sub_pd(s, e));
      e = _mm256_blend_pd(e, st, 0b1000);

      st = _mm256_permute4x64_pd(s, 0b10010011);
      b = _mm256_blend_pd(st, _mm256_loadu_pd(&err[i + 10]), 0b0001);

      s = _mm256_add_pd(e, b);
      e = _mm256_sub_pd(b, _mm256_sub_pd(s, e));

      _mm256_storeu_pd(s_array, s);
      r[i + 3] = s_array[3];

      int offset = i + 11;

      for (int j = i + 4; j < vec_sum_limit; j++)
      {
        st = _mm256_permute4x64_pd(s, 0b10010011);
        b = _mm256_blend_pd(st, _mm256_loadu_pd(&err[offset]), 0b0001);

        s = _mm256_add_pd(e, b);
        e = _mm256_sub_pd(b, _mm256_sub_pd(s, e));

        st = _mm256_permute4x64_pd(s, 0b00000011);
        _mm256_storeu_pd(&err[j], st);

        offset++;
      }
    }
  }

  for (; i < vec_sum_limit; i++)
  {
    // vecSumErr
    double e = err[i];

    for (int j = i; j < vec_sum_limit; j++)
    {
      double b = err[j + 1];

      double s = e + b;
      e = b - (s - e);

      err[j] = s;
    }

    err[vec_sum_limit] = e;

    r[i] = err[i];
  }

  r[vec_sum_limit] = err[vec_sum_limit];

  if (k == size_of_x)
  {
    r[size_of_x - 1] = err[size_of_x - 1];
  }

  return;
}
