inline void vecSumErrBranch2(double *e, int n, int m, double *f)
{
  double err = e[0];

  int l = m - 1;
  int j = 0;

  for (int i = 1; i <= n - 1; i++)
  {
    double b = e[i];

    double s = err + b;
    double t = s - b;

    err = (err - t) + (b - (s - t));

    f[j] = s;

    if (err == 0)
    {
      err = s;
    }
    else if (j >= l)
    {
      return;
    }
    else
    {
      j++;
    }
  }

  if (err != 0 && j < m)
  {
    f[j] = err;
  }
}

void renormalizationalgorithm4(double x[], int size_of_x, double f[], int m)
{
  int length = size_of_x;

  double *e = (double *)alloca(length * sizeof(double));
  double *temp = (double *)alloca(length * sizeof(double));
  double *g = (double *)alloca(m * sizeof(double));

  vecSum(x, e, size_of_x);
  vecSumErrBranch2(e, length, m + 1, temp);

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
