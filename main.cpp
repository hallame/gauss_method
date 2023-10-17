#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>

void ShowVector(int n, double * vec);
void PryamoiHod(int n, double **a, double *b);
void ObratniHod(int n, double **a, double *b, double *x);
double* SumVect(int n, double *a, double *b);
double* RaznVect(int n, double *a, double *b);
double* Mat_X_Vect(int n, double **a, double *b);
double Norm(int n, double *b);


int main()
{
	const int n = 4;
    double **a, *b, *x;
;
        a = (double **)malloc(n*sizeof(double));
        b = (double *)malloc(n*sizeof(double));
        x = (double *)malloc(n*sizeof(double));
        for(int i = 0; i < n; i++)
        {
            a[i] = (double *)malloc(n*sizeof(double));
        }
        a[0][0]=7*M_PI;
        a[0][1]=4;
        a[0][2]=3;
        a[0][3]=3;

        a[1][0]=4;
        a[1][1]=8*M_PI;
        a[1][2]=3;
        a[1][3]=2;

        a[2][0]=4;
        a[2][1]=4;
        a[2][2]=7*M_PI;
        a[2][3]=1;

        a[3][0]=1;
        a[3][1]=2;
        a[3][2]=2;
        a[3][3]=5*M_PI;

        b[0]=2;
        b[1]=5;
        b[2]=2;
        b[3]=5;


        printf("Matrix A:\r\n");
        for(int i = 0; i < n; i++)
            ShowVector(n, a[i]);
        printf("Vector B:\r\n");
        ShowVector(n, b);

        printf("\n\nSolving on Gauss method\r\n");
        PryamoiHod(n, a, b);
        printf("Forvard Gauss course\r\n");                        
        printf("Matrix A:\r\n");
        for(int i = 0; i < n; i++)
            ShowVector(n, a[i]);
        printf("Vector B:\r\n");
        ShowVector(n, b);

        ObratniHod(n, a, b, x);

        printf("\n\nResults :\r\n");
        ShowVector(n, x);


 		double d = Norm(n, RaznVect(n, Mat_X_Vect(n, a, x), b));
 		printf("DELTA :%.5f\r\n", d);


        free((void *)a);
        free((void *)b);
        free((void *)x);


      return 0;
}

void ShowVector(int n, double * vec)
{
    for(int i = 0; i < n; i++)
        printf("%.5f ",vec[i]);
    printf("\r\n");
}

void PryamoiHod(int n, double **a, double *b)
{
        double v;
        for(int k = 0,i,j,im; k < n - 1; k++)
        {
                im = k;
                for(i = k + 1; i < n; i++)
                {
                        if(fabs(a[im][k]) < fabs(a[i][k]))
                        {
                                im = i;
                        }
                }
                if(im != k)
                {
                        for(j = 0; j < n; j++)
                        {
                                v                = a[im][j];
                                a[im][j] = a[k][j];
                                a[k][j]  = v;
                        }
                        v     = b[im];
                        b[im] = b[k];
                        b[k]  = v;
                }
                for(i = k + 1; i < n; i++)
                {
                        v               = 1.0*a[i][k]/a[k][k];
                        a[i][k] = 0;
                        b[i]    = b[i] - v*b[k];
                        if(v != 0)
                        for(j = k + 1; j < n; j++)
                        {
                                a[i][j] = a[i][j] - v*a[k][j];
                        }
                }
        }
}

void ObratniHod(int n, double **a, double *b, double *x)
{
        double s = 0;
        x[n - 1] = 1.0*b[n - 1]/a[n - 1][n - 1];
        for(int i = n - 2, j; 0 <= i; i--)
        {
                s = 0;
                for(j = i + 1; j < n; j++)
                {
                        s = s+a[i][j]*x[j];
                }
                x[i] = 1.0*(b[i] - s)/a[i][i];
        }
}

double* SumVect(int n, double *a, double *b)
{
	double *x = (double *)malloc(n*sizeof(double));
	for(int i=0; i<n; i++)
	{
		x[i]=a[i]+b[i];
	}
	return x;
}

double* RaznVect(int n, double *a, double *b)
{
	double *x = (double *)malloc(n*sizeof(double));
	for(int i=0; i<n; i++)
	{
		x[i]=a[i]-b[i];
	}
	return x;
}

double* Mat_X_Vect(int n, double **a, double *b)
{
	double *out = (double *)malloc(n*sizeof(double));
    for (int i=0; i<n; i++)
    {
        out[i]=0;
        for (int j=0; j<n; j++)
            *(out+i) += *((*(a+i))+j) * *(b+j);
    }
    return out;
}

double Norm(int n, double *b)
{
	double res=0;
	for(int i=0; i<n; i++)
	{
		res=res+b[i]*b[i];
	}
	return sqrt(res);
}
