void multiplication(double *a, double *b, double *r, const int sizea, const int sizeb, const int sizer)
{
  int k = sizea;

  double *err = (double *)alloca(2 * (k * k + 3 * k) * sizeof(double));

  for (int i = 0; i < 2 * k * k + 3 * k; i++)
  {
    err[i] = 0;
  }

  int ext_size = 2 * (k * k);
  double *r_ext = (double *)alloca(ext_size * sizeof(double));

  for (int i = 0; i < ext_size; i++)
  {
    r_ext[i] = 0;
  }

  twoMultFMA(a[0], b[0], &(r_ext[0]), &(err[0]));

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
      twoMultFMA(a[i], b[n - i], &(p[i]), &(e_tmp[i]));
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
      tmp[n + i] = err[i];
    }

    vecSum(tmp, tmp1, (n * n + n));

    /* write  tmp1 into r_ext[n], e[0:n^2 +n-1]*/
    r_ext[n] = tmp1[0];
    for (int i = 0; i <= n * n + n - 1; i++)
    {
      err[i] = tmp1[i + 1];
    }

    // writes err[0:n^2 +n-1],e_tmp[0:n] into err[0:n^2 +2n-1]
    for (int i = 0; i <= n * n - 1; i++)
    {
      err[n * n + n + i] = err[i];
    }

    for (int i = 0; i <= n; i++)
    {
      err[i] = e_tmp[i];
    }
  }

  double r_ext_k = r_ext[k];
  for (int i = 1; i < k; i++)
  {
    r_ext_k += a[i] * b[k - i];
  }
  r_ext[k] = r_ext_k;

  for (int i = 0; i < k * k; i++)
  {
    r_ext[i] += a[i] * b[k - i];
  }

  renormalization4(r_ext, k + 1, r, k);

  return;
}
