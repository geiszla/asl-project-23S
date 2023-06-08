inline void vecSumErrBranch2(double *e, int n, int m, double *f)
{
  int l = m - 1;
  int j = 0;

  double err = e[0];

  for (int i = 1; i < n; i++)
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

inline void vecSumErrBranch_fast(double *e, int n, int m, double *f)
{
  int l = m - 1;
  int j = 0;

  double err = e[0];

  for (int i = 1; i < n; i++)
  {
    double b = e[i];

    double s = err + b;
    err = b - (s - err);

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
