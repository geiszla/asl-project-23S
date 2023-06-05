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
      temp_array[n] = b[j];
      j++;
    }
    else
    {
      temp_array[n] = a[i];
      i++;
    }
  }

  double temp = temp_array[size_of_x - 1];

  if (size_of_x < 18)
  {
    for (int i = size_of_x - 2; i >= 0; i--)
    {
      double a = temp_array[i];
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
      s3 = temp_array[i] + temp;
      s_array[i] = s3;

      s2 = temp_array[i - 1] + s3;
      s_array[i - 1] = s2;

      s1 = temp_array[i - 2] + s2;
      s_array[i - 2] = s1;

      temp = temp_array[i - 3] + s1;
      s_array[i - 3] = temp;
    }

    int i;

    for (i = size_of_x - 2; i >= 3; i -= 4)
    {
      int start_index = i - 3;

      __m256d s = _mm256_loadu_pd(&s_array[start_index]);
      __m256d b = _mm256_loadu_pd(&s_array[start_index + 1]);

      __m256d a = _mm256_loadu_pd(&temp_array[start_index]);

      __m256d t = _mm256_sub_pd(s, b);
      __m256d e = _mm256_add_pd(_mm256_sub_pd(a, t), _mm256_sub_pd(b, _mm256_sub_pd(s, t)));

      _mm256_storeu_pd(&temp_array[start_index + 1], e);
    }

    for (; i >= 0; i--)
    {
      double a = temp_array[i];
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
