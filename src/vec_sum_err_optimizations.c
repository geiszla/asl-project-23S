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

inline void vecSumErr_fast(double *f, int n, double *g)
{
  int m = n - 1;

  double e = f[0];

  for (int i = 0; i < m; i++)
  {
    double b = f[i + 1];

    double s = e + b;
    e = b - (s - e);

    g[i] = s;
  }

  g[m] = e;
}
