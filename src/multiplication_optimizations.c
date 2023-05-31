void multiplication2(double *a, double *b, double *r, const int sizea, const int sizeb, const int sizer)
{
  int k = sizea;

  int err_size = k * k;
  double *err = (double *)alloca(err_size * sizeof(double));

  for (int i = 0; i < err_size; i++)
  {
    err[i] = 0;
  }

  int ext_size = 2 * (k * k);
  double *r_ext = (double *)alloca(ext_size * sizeof(double));

  for (int i = 0; i < ext_size; i++)
  {
    r_ext[i] = 0;
  }

  double x = a[0];
  double y = b[0];

  double pi = x * y;

  err[0] = fma(x, y, -pi);
  r_ext[0] = pi;

  for (int n = 1; n <= (k - 1); n++)
  {
    int tmp_size = 2 * (n + 1);

    double *e_tmp = (double *)alloca(tmp_size * sizeof(double));
    double *p = (double *)alloca(tmp_size * sizeof(double));

    for (int i = 0; i < tmp_size; i++)
    {
      e_tmp[i] = 0;
      p[i] = 0;
    }

    for (int i = 0; i <= n; i++)
    {
      x = a[i];
      y = b[n - i];

      pi = x * y;

      e_tmp[i] = fma(x, y, -pi);
      p[i] = pi;
    }

    double *tmp = (double *)alloca((n * n + n + 1) * sizeof(double));
    double *tmp1 = (double *)alloca((n * n + n + 1) * sizeof(double));

    for (int i = 0; i < (n * n + n + 1); i++)
    {
      tmp1[i] = 0;
    }

    // generate p[0:n], e[0:n^2-1] into tmp
    for (int i = 0; i <= n; i++)
    {
      tmp[i] = p[i];
    }

    for (int i = 0; i <= n * n - 1; i++)
    {
      tmp[n + 1 + i] = err[i];
    }

    vecSum(tmp, tmp1, (n * n + n));

    /* write  tmp1 into r_ext[n], e[0:n^2 +n-1]*/
    r_ext[n] = tmp1[0];

    for (int i = 0; i <= n * n + n - 1; i++)
    {
      err[i] = tmp1[i + 1];
    }

    // writes err[0:n^2 +n-1],e_tmp[0:n] into err[0:n^2 +2n-1]
    for (int i = n * n; i <= ((n * n) + 2 * n - 1); i++)
    {
      err[i] = e_tmp[i - n * n];
    }
  }

  for (int i = 1; i <= k - 1; i++)
  {
    r_ext[k] += a[i] * b[k - i];
  }

  for (int i = 0; i <= k * k - 1; i++)
  {
    r_ext[i] += err[i];
  }

  renormalization4(r_ext, k + 1, r, k);

  return;
}
