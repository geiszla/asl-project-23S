
/* multiplies 4 doubles at a time 
it's important that the input sizes of all doubles and the outputsizes are the same 
no cehck for this consition is performed*/
void fourtimesmultiplicationversion0(double *a0, double *b0, double *a1, double *b1, double *a2, double *b2,double *a3, double *b3,double *r0, double *r1, double *r2, double *r3, const int sizea, const int sizeb, const int sizer)
{
    multiplication(a0,b0,r0,sizea,sizeb,sizer);
    multiplication(a1,b1,r1,sizea,sizeb,sizer);
    multiplication(a2,b2,r2,sizea,sizeb,sizer);
    multiplication(a3,b3,r3,sizea,sizeb,sizer);
    return;
}
void fourtimesmultiplicationversion1(double *a0, double *b0, double *a1, double *b1, double *a2, double *b2,double *a3, double *b3,double *r0, double *r1, double *r2, double *r3, const int sizea, const int sizeb, const int sizer)
{
    int k = sizea;
    double *err = (double *)alloca(2*(sizea * sizea +1*sizea +1) * sizeof(double));
    for(int i = 0; i< 2*(sizea * sizea +1*sizea +1);i++){
        err[i] = 0;
    }
    double *r_ext = (double *)alloca((sizea  + 1) * sizeof(double));
    for (int i = 0; i <(sizea  + 1) ; i++)
    {
        r_ext[i] = 0.0;
    }
    multiplication(a0,b0,r0,sizea,sizeb,sizer);
    multiplication(a1,b1,r1,sizea,sizeb,sizer);
    multiplication(a2,b2,r2,sizea,sizeb,sizer);
    multiplication(a3,b3,r3,sizea,sizeb,sizer);
    return;
}