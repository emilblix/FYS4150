#include <iostream>
#include <armadillo>
#include <string>
#include <sstream>
using namespace std;
using namespace arma;
#include "time.h"

int main()
{
    // Time measurement
    clock_t start, finish;
    start = clock();

    double h,h_squared;
    int n;
    cout << "Enter value for n: ";
    cin >> n;
    h=1.0/(n+1);
    vec x= linspace(0,1,n+2);

    mat A=zeros<mat>(n,n);

    // Vector a, b og c in tridiagonal matrix
    vec tridiag_a= vec(n); tridiag_a.fill(-1);
    vec tridiag_b= vec(n); tridiag_b.fill(2);
    vec tridiag_c= tridiag_a;
    tridiag_a(0)=0;
    tridiag_c(n-1)=0;

    // Vector with initial function
    h_squared=h*h;
    vec init_b= h_squared*100*exp(-10*x);
    init_b(0)=0;init_b(n+1)=0; //Start and end point are unaffected by calculations
    vec b_marked=init_b.subvec(1,n);       // Copy of inital function to be used in later calculations

    // Forward sweep
    tridiag_c(0) = tridiag_c(0) / tridiag_b(0);
    init_b(1) = init_b(1) / tridiag_b(0);
    // The first and last element of init_func are outside the scope of the tridiagonal matrix, and are omitted from calculations
    for (int i=1;i<n;i++)
    {
        tridiag_c(i) = tridiag_c(i)/(tridiag_b(i) - tridiag_a(i)*tridiag_c(i-1));
        init_b(i+1) = (init_b(i+1)-tridiag_a(i) * init_b(i)) / (tridiag_b(i)-tridiag_a(i) * tridiag_c(i-1));
    }

    // Back substitution
    vec ans= init_b;
    for (int i=n-2;i>=0;i--)
    {
        ans(i+1)=init_b(i+1)-tridiag_c(i)*ans(i+2);
    }

    // End timing for first calculation method
    finish = clock();
    double time_normal = ((finish-start)/(double) CLOCKS_PER_SEC);
    cout << "Algorithm calc finished"<<endl;

    //Answer matrix for plotting of closed-form and calculated function
    vec closedform= 1.0-(1.0-exp(-10.0))*x-exp(-10.0*x);
    mat svar= mat(n+2,3);
    for (int i=0;i<n+2;++i)
    {
        svar(i,0)=h*i;
        svar(i,1)=closedform(i);
        svar(i,2)=ans(i);
    }

    // Saving answer matrix with file name including value of n
    std::ostringstream namein;
    namein << "Solution_" << n << ".dat";
    std::string filename = namein.str();
    svar.save(filename, raw_ascii);

    // Relative error
    vec error=log10(abs(ans-closedform)/closedform);
    error=error.subvec(1,n);

//-----------------------------------------------------------------------------------------
    // Gaussian elimination
    // The if loop is present to avoid having conflicts from missing variables furter down which
    // would have risen from commenting out this part at n>1000

    vec ans_Gaussian= zeros<vec>(n);
    if(n<1001)
    {
    // Creating matrix A
    mat A=zeros<mat>(n,n);
    A.diag(0).fill(2);
    A.diag(1).fill(-1);  // Setting upper and lower tridiagonal as -1
    A.diag(-1).fill(-1);

    start = clock () ;
    // Solving Ax = b

    ans_Gaussian= solve(A,b_marked);

    finish = clock();
    }
    else{start=0;finish=0;}

    double time_Gaussian = ((finish-start)/(double) CLOCKS_PER_SEC);

    cout<<"Gauss calc finished"<< endl;

    // Relative error for Gaussian elimination
    vec closedform_short=closedform.subvec(1,n);
    vec error_Gaussian=log10(abs(ans_Gaussian-closedform_short)/closedform_short);

//------------------------------------------------------------------------------------------
    // LU decomposition
    // Creating matrices for LU decomp
    mat P,L,U;
    vec ans_LU= zeros<vec>(n);

    if(n<1001)
    {
    mat A=zeros<mat>(n,n);
    A.diag(0).fill(2);
    A.diag(1).fill(-1);  // Setting upper and lower tridiagonal as -1
    A.diag(-1).fill(-1);

    start = clock () ;

    // Solving for LU decomp
    lu(L, U, P, A);
    vec y   = solve(L,b_marked);
    ans_LU  = solve(U,y);

    finish = clock();
    }
    else{start=0;finish=0;}

    double time_LU = ((finish-start)/(double) CLOCKS_PER_SEC);

    cout<<"LU calc finished"<< endl;

    // Relative error for LU decomp
    vec error_LU=log10(abs(ans_LU-closedform_short)/closedform_short);

//------------------------------------------------------------------------------------------
    //Answer matrix for comparison of closed-form and calculated functions
    mat svar_all= mat(n,5);
    for (int i=0;i<n;++i)
    {
        svar_all(i,0)=h*(i+1);
        svar_all(i,1)=closedform(i+1);
        svar_all(i,2)=ans(i+1);
        svar_all(i,3)=ans_Gaussian(i);
        svar_all(i,4)=ans_LU(i);
    }

    // Saving comparison matrix
    namein.str(""); namein.clear(); // Clearing previous value of the string
    namein << "Solution_all_" << n << ".dat";
    filename = namein.str();
    svar_all.save(filename, raw_ascii);


    // Saving relative error as a matrix for comparison
    mat errormat=randn<mat>(n,3);
    errormat.col(0)=error;
    errormat.col(1)=error_Gaussian;
    errormat.col(2)=error_LU;

    namein.str(""); namein.clear(); // Clearing previous value of the string
    namein << "Error_" << n << ".dat";
    filename = namein.str();
    errormat.save(filename, raw_ascii);

//------------------------------------------------------------------------------------------
    // Printing answers
    cout << "For n= " << n << ":" << endl;
    cout << "Max error algorithm = "<< error.max() << endl;
    cout << "Max error Gaussian = " << error_Gaussian.max() << endl;
    cout << "Max error LU = "       << error_LU.max() << endl;
    cout << "Computational time using algorithm: " << time_normal << endl;
    cout << "Computational time using Gaussian: " << time_Gaussian << endl;
    cout << "Computational time using LU decomposition: " << time_LU << endl;
    return 0;
}


