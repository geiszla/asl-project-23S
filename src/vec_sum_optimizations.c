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
    e_res[i] = e;
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

    e_res[i] = e;
  }

  e_res[0] = s;
}

void renormalizationalgorithm2(double x[], int size_of_x, double f[], int m)
{
  int length = size_of_x;

  double *e = (double *)alloca(length * sizeof(double));
  double *temp = (double *)alloca(length * sizeof(double));
  double *g = (double *)alloca(m * sizeof(double));

  vecSum3(x, e, size_of_x);
  vecSumErrBranch(e, length, m + 1, temp);

  for (int i = 0; i <= m - 2; i++)
  {
    vecSumErr(&temp[i], m, g);
    for (int b = 0; b < m; b++)
    {
      double tmp = g[b];
      temp[b + i] = tmp;
    }
    double t = temp[i];
    f[i] = t;
  }
  f[m - 1] = temp[m - 1];
}
