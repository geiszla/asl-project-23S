#include <immintrin.h>

inline void vecSum2(double *x, double *e_res, int in_out_size)
{
  double *s = (double *)alloca(in_out_size * sizeof(double));

  s[in_out_size - 1] = x[in_out_size - 1];

  for (int i = in_out_size - 2; i >= 0; i--)
  {
    double a = x[i];
    double b = s[i + 1];

    double s2 = a + b;
    double t = s2 - b;

    double e = (a - t) + (b - (s2 - t));

    s[i] = s2;
    e_res[i + 1] = e;
  }

  e_res[0] = s[0];
}

inline void vecSum3(double *x, double *e_res, int in_out_size)
{
  double s = x[in_out_size - 1];

  for (int i = in_out_size - 2; i >= 0; i--)
  {
    double a = x[i];
    double b = s;

    s = a + b;
    double t = s - b;

    double e = (a - t) + (b - (s - t));

    e_res[i + 1] = e;
  }

  e_res[0] = s;
}

// Best for input size < 20
inline void vecSum4(double *x, double *e_res, int in_out_size)
{
  double st;
  double s3;
  double s2;
  double s1;
  double s0 = x[in_out_size - 1];

  int i;

  for (i = in_out_size - 2; i >= 3; i -= 4)
  {
    st = s0;
    s3 = x[i] + s0;
    s2 = x[i - 1] + s3;
    s1 = x[i - 2] + s2;
    s0 = x[i - 3] + s1;

    __m256d s = _mm256_set_pd(s3, s2, s1, s0);
    __m256d b = _mm256_set_pd(st, s3, s2, s1);

    __m256d a = _mm256_loadu_pd(&x[i - 3]);

    __m256d t = _mm256_sub_pd(s, b);
    __m256d e = _mm256_add_pd(_mm256_sub_pd(a, t), _mm256_sub_pd(b, _mm256_sub_pd(s, t)));

    _mm256_storeu_pd(&e_res[i - 2], e);
  }

  for (; i >= 0; i--)
  {
    double a = x[i];
    double b = s0;

    s0 = a + b;
    double t = s0 - b;

    double e = (a - t) + (b - (s0 - t));

    e_res[i + 1] = e;
  }

  e_res[0] = s0;
}

// Best for input size >= 21, < 33
inline void vecSum5(double *x, double *e_res, int in_out_size)
{
  double *s_array = (double *)alloca(in_out_size * sizeof(double));

  double s3;
  double s2;
  double s1;

  double s0 = x[in_out_size - 1];
  s_array[in_out_size - 1] = s0;

  for (int i = in_out_size - 2; i >= 3; i -= 4)
  {
    s3 = x[i] + s0;
    s_array[i] = s3;

    s2 = x[i - 1] + s3;
    s_array[i - 1] = s2;

    s1 = x[i - 2] + s2;
    s_array[i - 2] = s1;

    s0 = x[i - 3] + s1;
    s_array[i - 3] = s0;
  }

  int i;

  for (i = in_out_size - 2; i >= 3; i -= 4)
  {
    int start_index = i - 3;

    __m256d s = _mm256_loadu_pd(&s_array[start_index]);
    __m256d b = _mm256_loadu_pd(&s_array[start_index + 1]);

    __m256d a = _mm256_loadu_pd(&x[start_index]);

    __m256d t = _mm256_sub_pd(s, b);
    __m256d e = _mm256_add_pd(_mm256_sub_pd(a, t), _mm256_sub_pd(b, _mm256_sub_pd(s, t)));

    _mm256_storeu_pd(&e_res[start_index + 1], e);
  }

  for (; i >= 0; i--)
  {
    double a = x[i];
    double b = s0;

    s0 = a + b;
    double t = s0 - b;

    double e = (a - t) + (b - (s0 - t));

    e_res[i + 1] = e;
  }

  e_res[0] = s0;
}

// Best for input size >= 33
inline void vecSum6(double *x, double *e_res, int in_out_size)
{
  double *s_array = (double *)alloca(in_out_size * sizeof(double));

  double s7;
  double s6;
  double s5;
  double s4;

  double s3;
  double s2;
  double s1;

  double s0 = x[in_out_size - 1];
  s_array[in_out_size - 1] = s0;

  int i;

  for (i = in_out_size - 2; i >= 7; i -= 8)
  {
    s7 = x[i] + s0;
    s_array[i] = s7;

    s6 = x[i - 1] + s7;
    s_array[i - 1] = s6;

    s5 = x[i - 2] + s6;
    s_array[i - 2] = s5;

    s4 = x[i - 3] + s5;
    s_array[i - 3] = s4;

    s3 = x[i - 4] + s4;
    s_array[i - 4] = s3;

    s2 = x[i - 5] + s3;
    s_array[i - 5] = s2;

    s1 = x[i - 6] + s2;
    s_array[i - 6] = s1;

    s0 = x[i - 7] + s1;
    s_array[i - 7] = s0;
  }

  for (i = in_out_size - 2; i >= 7; i -= 8)
  {
    // Note: calculating all s[c] first might make more ilp for vector instructions possible
    // Note: unrolling might also be possible
    int start_index = i - 3;
    int next_start_index = start_index - 4;

    __m256d s0 = _mm256_loadu_pd(&s_array[start_index]);
    __m256d b0 = _mm256_loadu_pd(&s_array[start_index + 1]);

    __m256d s1 = _mm256_loadu_pd(&s_array[next_start_index]);
    __m256d b1 = _mm256_loadu_pd(&s_array[next_start_index + 1]);

    __m256d a0 = _mm256_loadu_pd(&x[start_index]);
    __m256d a1 = _mm256_loadu_pd(&x[next_start_index]);

    __m256d t0 = _mm256_sub_pd(s0, b0);
    __m256d t01 = _mm256_sub_pd(s0, t0);
    __m256d t02 = _mm256_sub_pd(b0, t01);
    __m256d t03 = _mm256_sub_pd(a0, t0);

    __m256d e0 = _mm256_add_pd(t03, t02);

    __m256d t1 = _mm256_sub_pd(s1, b1);
    __m256d t11 = _mm256_sub_pd(s1, t1);
    __m256d t12 = _mm256_sub_pd(b1, t11);
    __m256d t13 = _mm256_sub_pd(a1, t1);

    __m256d e1 = _mm256_add_pd(t13, t12);

    _mm256_storeu_pd(&e_res[start_index + 1], e0);
    _mm256_storeu_pd(&e_res[next_start_index + 1], e1);
  }

  for (; i >= 0; i--)
  {
    double a = x[i];
    double b = s0;

    s0 = a + b;
    double t = s0 - b;

    double e = (a - t) + (b - (s0 - t));

    e_res[i + 1] = e;
  }

  e_res[0] = s0;
}

inline void vecSum3_fast(double *x, double *e_res, int in_out_size)
{
  double s = x[in_out_size - 1];

  for (int i = in_out_size - 2; i >= 0; i--)
  {
    double a = x[i];
    double b = s;

    s = a + b;
    double e = b - (s - a);

    e_res[i + 1] = e;
  }

  e_res[0] = s;
}

inline void vecSum4_fast(double *x, double *e_res, int in_out_size)
{
  double st;
  double s3;
  double s2;
  double s1;
  double s0 = x[in_out_size - 1];

  int i;

  for (i = in_out_size - 2; i >= 3; i -= 4)
  {
    st = s0;
    s3 = x[i] + s0;
    s2 = x[i - 1] + s3;
    s1 = x[i - 2] + s2;
    s0 = x[i - 3] + s1;

    __m256d a = _mm256_loadu_pd(&x[i - 3]);
    __m256d b = _mm256_set_pd(st, s3, s2, s1);

    __m256d s = _mm256_set_pd(s3, s2, s1, s0);

    __m256d e = _mm256_sub_pd(b, _mm256_sub_pd(s, a));

    _mm256_storeu_pd(&e_res[i - 3], e);
  }

  for (; i >= 0; i--)
  {
    double a = x[i];
    double b = s0;

    s0 = a + b;
    double e = b - (s0 - a);

    e_res[i + 1] = e;
  }

  e_res[0] = s0;
}

inline void vecSum5_fast(double *x, double *e_res, int in_out_size)
{
  double *s_array = (double *)alloca(in_out_size * sizeof(double));

  double s3;
  double s2;
  double s1;

  double s0 = x[in_out_size - 1];
  s_array[in_out_size - 1] = s0;

  for (int i = in_out_size - 2; i >= 3; i -= 4)
  {
    s3 = x[i] + s0;
    s_array[i] = s3;

    s2 = x[i - 1] + s3;
    s_array[i - 1] = s2;

    s1 = x[i - 2] + s2;
    s_array[i - 2] = s1;

    s0 = x[i - 3] + s1;
    s_array[i - 3] = s0;
  }

  int i;

  for (i = in_out_size - 2; i >= 3; i -= 4)
  {
    int start_index = i - 3;

    __m256d a = _mm256_loadu_pd(&x[start_index]);
    __m256d b = _mm256_loadu_pd(&s_array[start_index + 1]);

    __m256d s = _mm256_loadu_pd(&s_array[start_index]);

    __m256d e = _mm256_sub_pd(b, _mm256_sub_pd(s, a));

    _mm256_storeu_pd(&e_res[start_index + 1], e);
  }

  for (; i >= 0; i--)
  {
    double a = x[i];
    double b = s0;

    s0 = a + b;
    double e = b - (s0 - a);

    e_res[i + 1] = e;
  }

  e_res[0] = s0;
}

inline void vecSum6_fast(double *x, double *e_res, int in_out_size)
{
  double *s_array = (double *)alloca(in_out_size * sizeof(double));

  double s7;
  double s6;
  double s5;
  double s4;

  double s3;
  double s2;
  double s1;

  double s0 = x[in_out_size - 1];
  s_array[in_out_size - 1] = s0;

  int i;

  for (i = in_out_size - 2; i >= 7; i -= 8)
  {
    s7 = x[i] + s0;
    s_array[i] = s7;

    s6 = x[i - 1] + s7;
    s_array[i - 1] = s6;

    s5 = x[i - 2] + s6;
    s_array[i - 2] = s5;

    s4 = x[i - 3] + s5;
    s_array[i - 3] = s4;

    s3 = x[i - 4] + s4;
    s_array[i - 4] = s3;

    s2 = x[i - 5] + s3;
    s_array[i - 5] = s2;

    s1 = x[i - 6] + s2;
    s_array[i - 6] = s1;

    s0 = x[i - 7] + s1;
    s_array[i - 7] = s0;
  }

  for (i = in_out_size - 2; i >= 7; i -= 8)
  {
    int start_index = i - 3;
    int next_start_index = start_index - 4;

    __m256d a0 = _mm256_loadu_pd(&x[start_index]);
    __m256d a1 = _mm256_loadu_pd(&x[next_start_index]);

    __m256d b0 = _mm256_loadu_pd(&s_array[start_index + 1]);
    __m256d b1 = _mm256_loadu_pd(&s_array[next_start_index + 1]);

    __m256d s0 = _mm256_loadu_pd(&s_array[start_index]);
    __m256d s1 = _mm256_loadu_pd(&s_array[next_start_index]);

    __m256d e0 = _mm256_sub_pd(b0, _mm256_sub_pd(s0, a0));
    __m256d e1 = _mm256_sub_pd(b1, _mm256_sub_pd(s1, a1));

    _mm256_storeu_pd(&e_res[start_index + 1], e0);
    _mm256_storeu_pd(&e_res[next_start_index + 1], e1);
  }

  for (; i >= 0; i--)
  {
    double a = x[i];
    double b = s0;

    s0 = a + b;
    double e = b - (s0 - a);

    e_res[i + 1] = e;
  }

  e_res[0] = s0;
}

void renormalizationalgorithm2(double x[], int size_of_x, double f[], int m)
{
  double *err = (double *)alloca((size_of_x) * sizeof(double));
  double *f_tmp = (double *)alloca((m + 1) * sizeof(double));

  for (int i = 0; i <= m; i++)
  {
    f_tmp[i] = 0;
  }

  vecSum6_fast(x, err, size_of_x);
  vecSumErrBranch(err, size_of_x, m + 1, f_tmp);

  for (int i = 0; i <= (m - 2); i++)
  {
    vecSumErr(&(f_tmp[i]), m - i + 1, &(f_tmp[i]));
    f[i] = f_tmp[i];
  }

  f[m - 1] = f_tmp[m - 1];

  return;
}
