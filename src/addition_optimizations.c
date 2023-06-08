void addition2(double *a, double *b, double *s, int length_a, int length_b, int length_result)
{
  double *tmp = (double *)alloca((length_a + length_b) * sizeof(double));

  merge(a, b, tmp, length_a, length_b);
  renormalizationalgorithm(tmp, length_a + length_b, s, length_result);

  return;
}

void addition3(double *a, double *b, double *s, int length_a, int length_b, int length_result)
{
  double *tmp = (double *)alloca((length_a + length_b) * sizeof(double));

  int i = 0;
  int j = 0;

  for (int n = 0; n < length_a + length_b; n++)
  {
    if (i == length_a || (j < length_b && fabs(b[j]) > fabs(a[i])))
    {
      tmp[n] = b[j];
      j++;
    }
    else
    {
      tmp[n] = a[i];
      i++;
    }
  }

  renormalizationalgorithm(tmp, length_a + length_b, s, length_result);
}

void addition4(double *a, double *b, double *f, int length_a, int length_b, int length_result)
{
  int size_of_x = length_a + length_b;

  double *merged_array = (double *)alloca((size_of_x + 1) * sizeof(double));
  double *temp_array = (double *)alloca((size_of_x + 1) * sizeof(double));

  int i = 0;
  int j = 0;

  for (int n = 0; n < size_of_x; n++)
  {
    double b_j = b[j];
    double a_i = a[i];

    double b_abs = b_j > 0 ? b_j : -b_j;
    double a_abs = a_i > 0 ? a_i : -a_i;

    if (i == length_a || (j < length_b && b_abs > a_abs))
    {
      merged_array[n] = b[j];
      j++;
    }
    else
    {
      merged_array[n] = a[i];
      i++;
    }
  }

  double temp = merged_array[size_of_x - 1];
  temp_array[size_of_x - 1] = temp;

  if (size_of_x < 20)
  {
    double st;
    double s3;
    double s2;
    double s1;

    int i;

    for (i = size_of_x - 2; i >= 3; i -= 4)
    {
      st = temp;
      s3 = merged_array[i] + temp;
      s2 = merged_array[i - 1] + s3;
      s1 = merged_array[i - 2] + s2;
      temp = merged_array[i - 3] + s1;

      __m256d s = _mm256_set_pd(s3, s2, s1, temp);
      __m256d b = _mm256_set_pd(st, s3, s2, s1);

      __m256d a = _mm256_loadu_pd(&merged_array[i - 3]);

      __m256d t = _mm256_sub_pd(s, b);
      __m256d e = _mm256_add_pd(_mm256_sub_pd(a, t), _mm256_sub_pd(b, _mm256_sub_pd(s, t)));

      _mm256_storeu_pd(&temp_array[i - 2], e);
    }

    for (; i >= 0; i--)
    {
      double a = merged_array[i];
      double b = temp;

      temp = a + b;
      double t = temp - b;

      double e = (a - t) + (b - (temp - t));

      temp_array[i + 1] = e;
    }
  }
  else if (size_of_x < 33)
  {
    double s3;
    double s2;
    double s1;

    temp_array[size_of_x - 1] = temp;

    for (int i = size_of_x - 2; i >= 3; i -= 4)
    {
      s3 = merged_array[i] + temp;
      temp_array[i] = s3;

      s2 = merged_array[i - 1] + s3;
      temp_array[i - 1] = s2;

      s1 = merged_array[i - 2] + s2;
      temp_array[i - 2] = s1;

      temp = merged_array[i - 3] + s1;
      temp_array[i - 3] = temp;
    }

    int i;

    for (i = size_of_x - 2; i >= 3; i -= 4)
    {
      int start_index = i - 3;

      __m256d s = _mm256_loadu_pd(&temp_array[start_index]);
      __m256d b = _mm256_loadu_pd(&temp_array[start_index + 1]);

      __m256d a = _mm256_loadu_pd(&merged_array[start_index]);

      __m256d t = _mm256_sub_pd(s, b);
      __m256d e = _mm256_add_pd(_mm256_sub_pd(a, t), _mm256_sub_pd(b, _mm256_sub_pd(s, t)));

      _mm256_storeu_pd(&temp_array[start_index + 1], e);
    }

    for (; i >= 0; i--)
    {
      double a = merged_array[i];
      double b = temp;

      temp = a + b;
      double t = temp - b;

      double e = (a - t) + (b - (temp - t));

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
      s7 = merged_array[i] + temp;
      temp_array[i] = s7;

      s6 = merged_array[i - 1] + s7;
      temp_array[i - 1] = s6;

      s5 = merged_array[i - 2] + s6;
      temp_array[i - 2] = s5;

      s4 = merged_array[i - 3] + s5;
      temp_array[i - 3] = s4;

      s3 = merged_array[i - 4] + s4;
      temp_array[i - 4] = s3;

      s2 = merged_array[i - 5] + s3;
      temp_array[i - 5] = s2;

      s1 = merged_array[i - 6] + s2;
      temp_array[i - 6] = s1;

      temp = merged_array[i - 7] + s1;
      temp_array[i - 7] = temp;
    }

    for (i = size_of_x - 2; i >= 7; i -= 8)
    {
      // Note: calculating all s[c] first might make more ilp for vector instructions possible
      // Note: unrolling might also be possible
      int start_index = i - 3;
      int next_start_index = start_index - 4;

      __m256d temp = _mm256_loadu_pd(&temp_array[start_index]);
      __m256d b0 = _mm256_loadu_pd(&temp_array[start_index + 1]);

      __m256d s1 = _mm256_loadu_pd(&temp_array[next_start_index]);
      __m256d b1 = _mm256_loadu_pd(&temp_array[next_start_index + 1]);

      __m256d a0 = _mm256_loadu_pd(&merged_array[start_index]);
      __m256d a1 = _mm256_loadu_pd(&merged_array[next_start_index]);

      __m256d t0 = _mm256_sub_pd(temp, b0);
      __m256d t01 = _mm256_sub_pd(temp, t0);
      __m256d t02 = _mm256_sub_pd(b0, t01);
      __m256d t03 = _mm256_sub_pd(a0, t0);

      __m256d e0 = _mm256_add_pd(t03, t02);

      __m256d t1 = _mm256_sub_pd(s1, b1);
      __m256d t11 = _mm256_sub_pd(s1, t1);
      __m256d t12 = _mm256_sub_pd(b1, t11);
      __m256d t13 = _mm256_sub_pd(a1, t1);

      __m256d e1 = _mm256_add_pd(t13, t12);

      _mm256_storeu_pd(&temp_array[start_index + 1], e0);
      _mm256_storeu_pd(&temp_array[next_start_index + 1], e1);
    }

    for (; i >= 0; i--)
    {
      double a = merged_array[i];
      double b = temp;

      temp = a + b;
      double t = temp - b;

      double e = (a - t) + (b - (temp - t));

      temp_array[i + 1] = e;
    }
  }

  temp_array[0] = temp;

  // vecSumErrBranch
  j = 0;

  double last_error = temp_array[1];
  if (last_error != 0)
  {
    temp = last_error;
    j++;
  }

  double j_limit = length_result < size_of_x ? length_result : length_result - 1;

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

  int vec_sum_limit = length_result < size_of_x ? length_result - 1 : size_of_x - 2;

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

  if (length_result == size_of_x)
  {
    f[size_of_x - 1] = temp_array[size_of_x - 1];
  }
}

void addition_fast(double *a, double *b, double *f, int length_a, int length_b, int length_result)
{
  int size_of_x = length_a + length_b;

  double *merged_array = (double *)alloca((size_of_x + 1) * sizeof(double));
  double *temp_array = (double *)alloca((size_of_x + 1) * sizeof(double));

  int i = 0;
  int j = 0;

  for (int n = 0; n < size_of_x; n++)
  {
    double b_j = b[j];
    double a_i = a[i];

    double b_abs = b_j > 0 ? b_j : -b_j;
    double a_abs = a_i > 0 ? a_i : -a_i;

    if (i == length_a || (j < length_b && b_abs > a_abs))
    {
      merged_array[n] = b[j];
      j++;
    }
    else
    {
      merged_array[n] = a[i];
      i++;
    }
  }

  double temp = merged_array[size_of_x - 1];
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
      s3 = merged_array[i] + temp;
      s2 = merged_array[i - 1] + s3;
      s1 = merged_array[i - 2] + s2;
      temp = merged_array[i - 3] + s1;

      __m256d s = _mm256_set_pd(s3, s2, s1, temp);
      __m256d b = _mm256_set_pd(st, s3, s2, s1);

      __m256d a = _mm256_loadu_pd(&merged_array[i - 3]);

      __m256d e = _mm256_sub_pd(b, _mm256_sub_pd(s, a));

      _mm256_storeu_pd(&temp_array[i - 2], e);
    }

    for (; i >= 0; i--)
    {
      double a = merged_array[i];
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
      double a = merged_array[i];
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
      s7 = merged_array[i] + temp;
      temp_array[i] = s7;

      s6 = merged_array[i - 1] + s7;
      temp_array[i - 1] = s6;

      s5 = merged_array[i - 2] + s6;
      temp_array[i - 2] = s5;

      s4 = merged_array[i - 3] + s5;
      temp_array[i - 3] = s4;

      s3 = merged_array[i - 4] + s4;
      temp_array[i - 4] = s3;

      s2 = merged_array[i - 5] + s3;
      temp_array[i - 5] = s2;

      s1 = merged_array[i - 6] + s2;
      temp_array[i - 6] = s1;

      temp = merged_array[i - 7] + s1;
      temp_array[i - 7] = temp;
    }

    for (i = size_of_x - 2; i >= 7; i -= 8)
    {
      int start_index = i - 3;
      int next_start_index = start_index - 4;

      __m256d a0 = _mm256_loadu_pd(&merged_array[start_index]);
      __m256d a1 = _mm256_loadu_pd(&merged_array[next_start_index]);

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
      double a = merged_array[i];
      double b = temp;

      temp = a + b;
      double e = b - (temp - a);

      temp_array[i + 1] = e;
    }
  }

  temp_array[0] = temp;

  // vecSumErrBranch
  j = 0;

  double last_error = temp_array[1];
  if (last_error != 0)
  {
    temp = last_error;
    j++;
  }

  double j_limit = length_result < size_of_x ? length_result : length_result - 1;

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

  int vec_sum_limit = length_result < size_of_x ? length_result - 1 : size_of_x - 2;

  for (int i = 0; i < vec_sum_limit; i++)
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

  if (length_result == size_of_x)
  {
    f[size_of_x - 1] = temp_array[size_of_x - 1];
  }
}
