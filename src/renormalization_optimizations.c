inline void renormalization2(double x[], int size_of_x, double f[], int m)
{
  double *err = (double *)alloca((size_of_x) * sizeof(double));
  double *f_tmp = (double *)alloca((m + 1) * sizeof(double));

  for (int i = 0; i <= m; i++)
  {
    f_tmp[i] = 0;
  }

  vecSum3(x, err, size_of_x);
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

  // vecSum
  double s = x[size_of_x - 1];

  for (int i = size_of_x - 2; i >= 0; i--)
  {
    double a = x[i];
    double b = s;

    s = a + b;
    double t = s - b;

    double e = (a - t) + (b - (s - t));

    err[i + 1] = e;
  }

  // vecSumErrBranch
  double err2 = s;

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

    for (int j = 0; j < m - i; j++)
    {
      double b = f_tmp[i + j + 1];

      double s = e + b;
      double t = s - b;

      e = (e - t) + (b - (s - t));

      f_tmp[i + j] = s;
    }

    f_tmp[m - i] = e;

    f[i] = f_tmp[i];
  }

  f[m - 1] = f_tmp[m - 1];
}

inline void renormalization4(double x[], int size_of_x, double f[], int m)
{
  double *err = (double *)alloca((size_of_x) * sizeof(double));
  double *f_tmp = (double *)alloca((m + 1) * sizeof(double));

  for (int i = 0; i <= m; i++)
  {
    f_tmp[i] = 0;
  }

  // vecSum
  double s = x[size_of_x - 1];

  // Note: could we parallelize the error calculations? (probably not)
  // Or do it recursively?
  for (int i = size_of_x - 2; i >= 0; i--)
  {
    double a = x[i];
    double b = s;

    s = a + b;
    double t = s - b;

    double e = (a - t) + (b - (s - t));

    err[i + 1] = e;
  }

  // vecSumErrBranch
  double err2 = s;

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
    int k = i;

    for (int j = 0; j < m - i; j++)
    {
      double b = f_tmp[k + 1];

      double s = e + b;
      double t = s - b;

      e = (e - t) + (b - (s - t));

      f_tmp[k] = s;
      k++;
    }

    f_tmp[m - i] = e;

    f[i] = f_tmp[i];
  }

  f[m - 1] = f_tmp[m - 1];
}
