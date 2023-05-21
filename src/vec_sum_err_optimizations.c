inline void vecSumErr2(double *f, int n, double *g)
{
  int m = n - 1;

  double e = f[0];

  for (int i = 0; i < m; i++)
  {
    double b = f[i + 1];

    double s = e + b;
    double t = s - b;

    e = (e - t) + (b - (s - t));

    g[i] = s;
  }

  g[m] = e;
}

void renormalizationalgorithm3(double x[], int size_of_x, double f[], int m)
{
  int length = size_of_x;

  double *e = (double *)alloca(length * sizeof(double));
  double *temp = (double *)alloca(length * sizeof(double));
  double *g = (double *)alloca(m * sizeof(double));

  vecSum(x, e, size_of_x);
  vecSumErrBranch(e, length, m + 1, temp);

  for (int i = 0; i <= m - 2; i++)
  {
    vecSumErr2(&temp[i], m, g);
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
