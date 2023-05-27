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

void renormalizationalgorithm2(double x[], int size_of_x, double f[], int m)
{
  double *err = (double *)alloca((size_of_x) * sizeof(double));
  double *f_tmp = (double *)alloca((m + 1) * sizeof(double));

  for (int i = 0; i <= m; i++)
  {
    f_tmp[i] = 0;
  }

  vecSum3(x, err, size_of_x);
  vecSumErrBranch(err, size_of_x, m + 1, f_tmp);

  for (int i = 0; i <= (m - 2); i++)
  {
    vecSumErr(&(f_tmp[i]), m - i + 1, &(f_tmp[i]));
    f[i] = f_tmp[i];
  }

  f[m - 1] = f_tmp[m - 1];

  return;
}
