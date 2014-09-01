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
    vec A= vec(n); A.fill(-1); A(0)=0;      // Vektor av lengde n der alle verdier bortsett fra den første er -1 (nedre tridiagonal)
    vec B= vec(n); B.fill(2);               // Vektor av lengde n der alle verdier er 2
    vec C= vec(n); C.fill(-1); C(n-1)=0;    // Vektor av lengde n der alle verdier bortsett fra den siste er -1 (øvre tridiagonal)
    vec U= vec(n);                          // Vektor for closed-form (analytiske) verdier
    vec V= vec(n); //V(0)=5;V(1)=3;V(2)=2;V(3)=5;V(4)=3; // testverdier for n=5

    for (int i=0;i<n;++i)           // Lager funksjonsvektoren
    {
        V(i)= 100  *exp(-10*h*i);
        U(i)= 10000*exp(-10*h*i);
    }
    //cout << "V:" << endl << V << endl;

    for (int i=1;i<n;++i)           // Redusere nedre trekant av "matrisen" til 0
    {
        gausstall=A(i)/B(i-1);
        B(i)=B(i)-C(i-1)*gausstall;
        //A(i)=0;//A(i)=A(i)-B(i-1)*gausstall;     // Skal være =0
        V(i)=V(i)-V(i-1)*gausstall;
    }

    for (int i=n-2;i>=0;--i)        // Redusere øvre trekant av "matrisen" til 0
    {
        gausstall=C(i)/B(i+1);
        //B(i)=B(i)-A(i+1)*gausstall;     // Skal være =B(i)
        //C(i)=0;//C(i)=C(i)-B(i+1)*gausstall;     // Skal være =0
        V(i)=V(i)-V(i+1)*gausstall;
    }

    for (int i=0;i<n;++i)           // Normere diagonalen
    {
        V(i)=V(i)/B(i);
        //B(i)=B(i)/B(i);
    }
    //cout << "A:" << endl << A << endl;
    //cout << "B:" << endl << B << endl;
    //cout << "C:" << endl << C << endl;
    //cout << "V.n_cols:" << V.n_cols << endl; // =1
    //cout << "V.n_rows:" << V.n_rows << endl; // = n
    mat svar= mat(V.n_rows,3);
    for (int i=0;i<n;++i)
    {
        svar(i,0)=h*i;
        svar(i,1)=U(i);
        svar(i,2)=V(i);
    }
    //cout << svar;
    //V.reshape(2,V.n_rows);
    svar.save("Solution.dat",raw_ascii);
    //cout << "V:" << endl << V << endl;
    //cout << A(2) << endl;
    return 0;
}

