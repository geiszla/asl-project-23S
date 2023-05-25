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

  return;
}

void addition4(double *a, double *b, double *f, int length_a, int length_b, int length_result)
{
  int size_of_x = length_a + length_b;
  double *tmp = (double *)alloca(size_of_x * sizeof(double));

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
      tmp[n] = b[j];
      j++;
    }
    else
    {
      tmp[n] = a[i];
      i++;
    }
  }

  double *err = (double *)alloca((size_of_x) * sizeof(double));
  double *f_tmp = (double *)alloca((length_result + 1) * sizeof(double));

  for (int i = 0; i <= length_result; i++)
  {
    f_tmp[i] = 0;
  }

  // vecSum
  double s = tmp[size_of_x - 1];

  // Note: could we parallelize the error calculations? (probably not)
  // Or do it recursively?
  for (int i = size_of_x - 2; i >= 0; i--)
  {
    double a = tmp[i];
    double b = s;

    s = a + b;
    double t = s - b;

    double e = (a - t) + (b - (s - t));

    err[i + 1] = e;
  }

  // vecSumErrBranch
  double err2 = s;

  j = 0;

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
    else if (j >= length_result)
    {
      break;
    }
    else
    {
      j++;
    }
  }

  if (err2 != 0 && j < length_result)
  {
    f_tmp[j] = err2;
  }

  for (int i = 0; i < length_result - 1; i++)
  {
    // vecSumErr
    double e = f_tmp[i];
    int k = i;

    for (int j = 0; j < length_result - i; j++)
    {
      double b = f_tmp[k + 1];

      double s = e + b;
      double t = s - b;

      e = (e - t) + (b - (s - t));

      f_tmp[k] = s;
      k++;
    }

    f_tmp[length_result - i] = e;

    f[i] = f_tmp[i];
  }

  f[length_result - 1] = f_tmp[length_result - 1];

  return;
}
