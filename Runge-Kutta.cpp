#include <iostream>
#include <cmath>
#include <random>
#include <ctime>
using namespace std;


// function to display array
void printArray(double** arr, int row, int col)
{
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            cout << arr[i][j]<<"\t";
        }
        cout << endl;
    }
    cout << endl;
}

double** Multiply_Matrices(double** M1, int row_1, int col_1,
                           double** M2, int row_2, int col_2)
{
    double** arr = new double*[row_1];

    if (col_1 != row_2)
    {
        cout << "Error: Matrices not compatibles" << endl;
    }
    else
    {
        for (int i = 0; i < row_1; ++i)
        {
            arr[i] = new double[col_2];
            for (int j = 0; j < col_2; ++j)
            {
                for (int k = 0; k < col_1 ; ++k)
                arr[i][j] += M1[i][k] * M2[k][j];
            }
        }
    }
    return arr;
}

// function to initialize and returning array in a Guassian distribution
double** Gaussian(double c, double s_c, int row, int col)
{
    double** arr = new double*[row];

    std::default_random_engine gen{static_cast<long unsigned int>(time(0))};
    std::normal_distribution<double> d(c,s_c);

    for (int i = 0; i < row; ++i)
    {
        arr[i] = new double[col];
        for (int j = 0; j < col; ++j)
        {
            arr[i][j] = d(gen);
        }
    }
    return arr;
}

// function to initialize and returning array in a uniform distribution
double** Random(double min, double max, int row, int col)
{
    double** arr = new double*[row];

    std::default_random_engine gen{static_cast<long unsigned int>(time(0))};
    std::uniform_real_distribution<double> dis(min, max);

    for (int i = 0; i < row; ++i)
    {
        arr[i] = new double[col];
        for (int j = 0; j < col; ++j)
        {
            arr[i][j] = dis(gen);
        }
    }
    return arr;
}

double** Transpose(double** M1, int row, int col)
{
    double** arr = new double*[col];

    for (int i = 0; i < col; ++i)
    {
        arr[i] = new double[row];
        for (int j = 0; j < row; ++j)
        {
            arr[i][j] = M1[j][i];
        }
    }
    return arr;
}

// function to initialize and returning array in a uniform distribution
double** fi(double** ni, double** Ki, double** ru, double** mi, double** cu,
            double t, double fac, int N, int M)
{
    double** arr = new double*[N];
    double** cu_v; double** ni_v;

    ni_v = ni;

    for(int i = 0; i < N; ++i)
    {
        for(int j = 0; j < 1; ++j)
        {
            ni_v[i][j] =  ni_v[i][j] + (fac)*Ki[i][j];
        }
    }

    cu_v   = cu;
    cu_v   = Multiply_Matrices(cu_v, N, M, ru, M, 1);

    for(int i = 0; i < N; ++i)
    {
        arr[i] = new double[1];
        for(int j = 0; j < 1; ++j)
        {
            arr[i][j] = cu_v[i][j] -  mi[i][j];
            arr[i][j] = arr[i][j]*ni_v[i][j];
        }
    }
    return arr;
}

// function to initialize and returning array in a uniform distribution
double** fu(double** ni, double** ru, double** Ku, double** ku, double** ci,
            double au, double t, double fac, int N, int M)
{
    double** arr = new double*[M];
    double** ci_v; double** ru_v;

    ru_v = ru;

    for(int i = 0; i < M; ++i)
    {
        for(int j = 0; j < 1; ++j)
        {
            ru_v[i][j] =  ru_v[i][j] + (fac)*Ku[i][j];
        }
    }

    ci_v   = ci;
    ci_v   = Multiply_Matrices(ci_v, M, N, ni, N, 1);

    for(int i = 0; i < M; ++i)
    {
        arr[i] = new double[1];
        for(int j = 0; j < 1; ++j)
        {
            arr[i][j] = au*(ku[i][j]-ru[i][j]) - ci_v[i][j];
            arr[i][j] = arr[i][j]*ru_v[i][j];
        }
    }
    return arr;
}

// Driver Code
int main()
{
    int ite;
    double v,c,s_c,k,s_k,m,t, tsim, h, au;



    //Establish parameters for the simulation
    const int N = 500;
    v    = 10;
    const int M = N/v ;
    c    = 1.0;
    s_c  = sqrt(1.0);
    k    = 5.0;
    s_k  = 0.0;
    m    = 1.0;
    au   = 1.0;

    //Establish  parameters for time
    t    = 0.0;
    tsim = 500;
    h    = 0.001;
    ite  = tsim/h;


    double** ru;  double** ku;
    double** ni;  double** mi;
    double** cu;  double** ci;
    double** k1i; double** k2i;   double** k3i; double** k4i;
    double** k1u; double** k2u;   double** k3u; double** k4u;


    ru   = Random(0, k, M, 1);
    ku   = Gaussian(k,s_k,M,1);

    ni   = Random(0,10, N, 1);
    mi   = Gaussian(m,0,N,1);

    cu   = Gaussian(c/N,s_c/sqrt(N),N,M);
    ci   = Transpose(cu,N,M);


    for(int i=0; i<ite; ++i)
    {
        k1i = fi(ni, ni, ru, mi, cu, t, 0, N, M);
        k2i = fi(ni, k1i, ru, mi, cu, t+h/2, h/2, N, M);
        k3i = fi(ni, k2i, ru, mi, cu, t+h/2, h/2, N, M);
        k4i = fi(ni, k3i, ru, mi, cu, t+h, h, N, M);

        k1u = fu(ni, ru, ru, ku, ci, au, t, 0, N, M);
        k2u = fu(ni, ru, k1u, ku, ci, au, t+h/2, h/2, N, M);
        k3u = fu(ni, ru, k2u, ku, ci, au, t+h/2, h/2, N, M);
        k4u = fu(ni, ru, k1u, ku, ci, au, t+h, h, N, M);

        for(int j=0; j<N; ++j)
        {
            for(int k=0; k<1; ++k)
            {
                ni[j][k] = ni[j][k] + h*(k1i[j][k]+2*k2i[j][k]+2*k3i[j][k]+k4i[j][k])/6;
            }
        }

        for(int j=0; j<M; ++j)
        {
            for(int k=0; k<1; ++k)
            {
                ru[j][k] = ru[j][k] + h*(k1u[j][k]+2*k2u[j][k]+2*k3u[j][k]+k4u[j][k])/6;
            }
        }

        t = t + h;
    }

    printArray(ni,N,1);


    return 0;
}
