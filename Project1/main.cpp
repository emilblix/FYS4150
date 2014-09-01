#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main()
{
    double n,gausstall,h;
    //n=5;                                    // Testverdi for n
    cout << "Enter value for n: ";
    cin >> n;
    h=1/(n+1);
    vec A= vec(n); A.fill(-1); A(0)=0;      // Lager en vektor av lengde n der alle verdier bortsett fra den første er -1 (nedre tridiagonal)
    vec B= vec(n); B.fill(2);               // Lager en vektor av lengde n der alle verdier er 2
    vec C= vec(n); C.fill(-1); C(n-1)=0;    // Lager en vektor av lengde n der alle verdier bortsett fra den siste er -1 (øvre tridiagonal)
    vec X= vec(n); //X(0)=5;X(1)=3;X(2)=2;X(3)=5;X(4)=3; // testverdier for n=5

    for (int i=0;i<n;++i)           // Lager funksjonsvektoren
    {
        X(i)= 100*exp(-10*h*i);
    }
    //cout << "X:" << endl << X << endl;

    for (int i=1;i<n;++i)           // Redusere nedre trekant av "matrisen" til 0
    {
        gausstall=A(i)/B(i-1);
        B(i)=B(i)-C(i-1)*gausstall;
        //A(i)=0;//A(i)=A(i)-B(i-1)*gausstall;     // Skal være =0
        X(i)=X(i)-X(i-1)*gausstall;
    }

    for (int i=n-2;i>=0;--i)        // Redusere øvre trekant av "matrisen" til 0
    {
        gausstall=C(i)/B(i+1);
        //B(i)=B(i)-A(i+1)*gausstall;     // Skal være =B(i)
        //C(i)=0;//C(i)=C(i)-B(i+1)*gausstall;     // Skal være =0
        X(i)=X(i)-X(i+1)*gausstall;
    }

    for (int i=0;i<n;++i)           // Normere diagonalen
    {
        X(i)=X(i)/B(i);
        //B(i)=B(i)/B(i);
    }
    //cout << "A:" << endl << A << endl;
    //cout << "B:" << endl << B << endl;
    //cout << "C:" << endl << C << endl;
    cout << "X.n_cols:" << X.n_cols << endl;
    cout << "X.n_rows:" << X.n_rows << endl;

    X.reshape(2,X.n_rows);
    X.save("Solution.dat",raw_ascii);
    //cout << "X:" << endl << X << endl;
    //cout << A(2) << endl;
    return 0;
}

