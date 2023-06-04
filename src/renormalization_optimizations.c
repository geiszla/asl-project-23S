inline void renormalization2(double x[], int size_of_x, double f[], int m)
{
  double *err = (double *)alloca((size_of_x) * sizeof(double));
  double *f_tmp = (double *)alloca((m + 1) * sizeof(double));

  for (int i = 0; i <= m; i++)
  {
    f_tmp[i] = 0;
  }

  if (size_of_x < 18)
  {
    vecSum3(x, err, size_of_x);
  }
  else
  {
    vecSum5(x, err, size_of_x);
  }

  vecSumErrBranch2(err, size_of_x, m + 1, f_tmp);

  for (int i = 0; i <= (m - 2); i++)
  {
    vecSumErr2(&(f_tmp[i]), m - i + 1, &(f_tmp[i]));

    f[i] = f_tmp[i];
  }

  f[m - 1] = f_tmp[m - 1];
}

inline void renormalization3(double x[], int size_of_x, double f[], int m)
{
  double *err = (double *)alloca((size_of_x) * sizeof(double));
  double *f_tmp = (double *)alloca((m + 1) * sizeof(double));

  for (int i = 0; i <= m; i++)
  {
    f_tmp[i] = 0;
  }

  if (size_of_x < 18)
  {
    vecSum3(x, err, size_of_x);
  }
  else
  {
    vecSum5(x, err, size_of_x);
  }

  // vecSumErrBranch
  double err2 = err[0];

  int j = 0;

  for (int i = 1; i <= size_of_x - 1; i++)
  {
    double b = err[i];

    double s = err2 + b;
    double t = s - b;

    err2 = (err2 - t) + (b - (s - t));

    f_tmp[j] = s;

    if (err2 == 0)
    {
      err2 = s;
    }
    else if (j >= m)
    {
      break;
    }
    else
    {
      j++;
    }
  }

  if (err2 != 0 && j < m)
  {
    f_tmp[j] = err2;
  }

  for (int i = 0; i < m - 1; i++)
  {
    // vecSumErr
    double e = f_tmp[i];

    for (int j = i; j < m; j++)
    {
      double b = f_tmp[j + 1];

      double s = e + b;
      double t = s - b;

      e = (e - t) + (b - (s - t));

      f_tmp[j] = s;
    }

    f_tmp[m] = e;

    f[i] = f_tmp[i];
  }

  f[m - 1] = f_tmp[m - 1];
}

inline void renormalization4(double x[], int size_of_x, double f[], int m)
{
  double *temp_array = (double *)alloca((size_of_x + 1) * sizeof(double));

  double temp = x[size_of_x - 1];

  if (size_of_x < 18)
  {
    for (int i = size_of_x - 2; i >= 0; i--)
    {
      double a = x[i];
      double b = temp;

      temp = a + b;
      double t = temp - b;

      double e = (a - t) + (b - (temp - t));

      temp_array[i + 1] = e;
    }
  }
  else
  {
    double *s_array = (double *)alloca(size_of_x * sizeof(double));

    double s3;
    double s2;
    double s1;

    s_array[size_of_x - 1] = temp;

    for (int i = size_of_x - 2; i >= 3; i -= 4)
    {
      s3 = x[i] + temp;
      s_array[i] = s3;

      s2 = x[i - 1] + s3;
      s_array[i - 1] = s2;

      s1 = x[i - 2] + s2;
      s_array[i - 2] = s1;

      temp = x[i - 3] + s1;
      s_array[i - 3] = temp;
    }

    int i;

    for (i = size_of_x - 2; i >= 3; i -= 4)
    {
      int start_index = i - 3;

      __m256d s = _mm256_loadu_pd(&s_array[start_index]);
      __m256d b = _mm256_loadu_pd(&s_array[start_index + 1]);

      __m256d a = _mm256_loadu_pd(&x[start_index]);

      __m256d t = _mm256_sub_pd(s, b);
      __m256d e = _mm256_add_pd(_mm256_sub_pd(a, t), _mm256_sub_pd(b, _mm256_sub_pd(s, t)));

      _mm256_storeu_pd(&temp_array[start_index + 1], e);
    }

    for (; i >= 0; i--)
    {
      double a = x[i];
      double b = temp;

      temp = a + b;
      double t = temp - b;

      double e = (a - t) + (b - (temp - t));

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
    double t = s - b;

    temp = (temp - t) + (b - (s - t));

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

  for (int i = 0; i < vec_sum_limit; i++)
  {
    // vecSumErr
    double e = temp_array[i];

    for (int j = i; j < vec_sum_limit; j++)
    {
      double b = temp_array[j + 1];

      double s = e + b;
      double t = s - b;

      e = (e - t) + (b - (s - t));

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
