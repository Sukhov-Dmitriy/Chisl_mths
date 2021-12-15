#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
using namespace std;
const double eps = 0.0000001;
const double AA = 0, BB = 1;// пределы интегрирования


void print(const vector <double> &x)
{
    for (int i = 0; i < x.size(); i++)
    {
        cout<< x[i] << "   ";
    }
    cout<<endl;
}
void print(const vector <vector<double>> &x)
{
    for (int i = 0; i < x.size(); i++)
    {
        for (int j = 0; j < x[i].size(); j++)
        {
        cout<< x[i][j] << "   ";
        }
        cout<<endl;
    }
    cout<<endl;
}


vector<double> smult(double a, const vector <double> &b)//умножение вектора на скаляр
{
    vector <double> ab(b.size());
    for(int i = 0; i < b.size(); i++)
        ab[i] = a*b[i];
    return ab;

}

vector<double> sum(const vector <double> &a, const vector <double> &b)//сумма векторов
{
    vector <double> ab(b.size());
    for(int i = 0; i < a.size(); i++)
        ab[i] = a[i] + b[i];
    return ab;

}

vector<vector<double>> invers(vector<vector<double>> A, int N)// поиск обратной матрицы методом Гаусса
{
    int k = 0, p = 0;
    int mu = 0;
    for(int y = 0; y < N; y++)//приписывание единичной матрицы
    {
        for(int z = 0; z <N; z++)
        {
            A[y].push_back(0.0);
            if(y == z)
            A[y][z+N] = 1.0;
        }
    }
    vector <double> temp(2*N);
    for(int i = 0; i < N; i++)//вверхтреугольность
    {
        if(abs(A[i][i]) < pow(eps,3))
        {
            k = i + 1;
            while(A[k][i] < pow(eps,3) && k < N)
            {
                k++;
                if(k == N)//проверяю на равенство определителя нулю
                {
                 cout<<"Err. Nulevoi determinant u matrici vtor proizvodnih"<<endl;
                 exit(-1);
                }

            }
            temp = A[i];
            A[i] = A[k];
            A[k] = temp;
        }
        p = i + 1;
        while(p < N)
        {
            A[p] = sum(A[p],smult(-(A[p][i]/A[i][i]),A[i]));
            p++;
        }

    }
    for(int i = N-1; i >= 0; i--)//диагонализация
    {
        p = i - 1;
        while(p >= 0)
        {
            A[p] = sum(A[p],smult(-(A[p][i]/A[i][i]),A[i]));
            p--;
        }
    }
    for(int i = 0; i < N; i++)//нормизация
    {
        A[i] = smult(1/A[i][i],A[i]);
    }
    for(int nm = 0; nm < N; nm++)//отделяем обратную матрицу
     {
         A[nm].erase(A[nm].begin(), A[nm].end()-N);

     }
    return A;
}

double math_func (double x)
{
    //здесь вбивается функция
    return exp(-x);
}
double math_Kerr (double x, double t)
{   //ядро
    return 0.5*x*exp(t);
}

double math_real_u(double x)
{
    return x+exp(-x);
}

vector <double> u_n_t (int h)
{
    double step = (BB - AA)/h;
    vector <double> b(h);
    vector <double> u(h);
    vector < vector <double> > A (h, vector <double> (h,0));
    vector < vector <double> > invA (h, vector <double> (h,0));
    for(int i = 0; i < h; i++)
    {
        b[i] = math_func(AA + i*step + step/2);
    }
    for(int i = 0; i < h; i++)
    {
        for(int j = 0; j < h; j++)
        {
            if(i == j)
            {
               A[i][j] = (1 - step*math_Kerr(AA + i*step + step/2,AA + j*step + step/2));
            }
            else
            {
                A[i][j] = - step*math_Kerr(AA + i*step + step/2,AA + j*step + step/2);
            }
        }
    }
    invA = invers(A,h);
    for(int i = 0; i < h; i++)
    {   u[i] = 0;
        for(int j = 0; j < h; j++)
        {
            u[i] += invA[i][j]* b[j];
        }
    }
    return u;
}
double dif_u(double x, vector <double> u_prs, vector <double> u_nxt)
{
    double zun = 0, zup = 0;

    double stepn = (BB - AA)/u_nxt.size();
    double stepp = (BB - AA)/u_prs.size();
    for(int i = 0; i < u_nxt.size(); i++)
        {zun += u_nxt[i]*stepn*math_Kerr(x,AA + i*stepn + stepn/2);}
    for(int i = 0; i < u_prs.size(); i++)
        {zup += u_prs[i]*stepp*math_Kerr(x,AA + i*stepp + stepp/2);}
    return (zun - zup)*(zun - zup);
}
double gauss3(double a,double b, vector <double> u_prs, vector <double> u_nxt)//3 узла
{
    double A1, A2, A3, sum;
    A1 = dif_u((a+b)/2.0 - 0.7745967*(b-a)/2,u_prs,u_nxt);
    A2 = dif_u((a+b)/2.0,u_prs,u_nxt);
    A3 = dif_u((a+b)/2.0 + 0.7745967*(b-a)/2,u_prs,u_nxt);
    //cout << " A1 A2 A3 "<< A1<<" " << A2<< " " << A3<< endl;
    sum = 0.5*(b-a)*((5.0/9.0)*A1 + (8.0/9.0)*A2 + (5.0/9.0)*A3);
    //cout<<"sum "<< sum<< endl;
    return sum;
}
double gauss4(double a,double b, vector <double> u_prs, vector <double> u_nxt)//4 узла
{
    double A1, A2, A3, A4, sum;

    A1 = dif_u((a+b)/2.0 - 0.8611363*(b-a)/2,u_prs,u_nxt);
    A1 = dif_u((a+b)/2.0 - 0.3399810*(b-a)/2,u_prs,u_nxt);
    A3 = dif_u((a+b)/2.0 + 0.3399810*(b-a)/2,u_prs,u_nxt);
    A1 = dif_u((a+b)/2.0 + 0.8611363*(b-a)/2,u_prs,u_nxt);
   // cout << " A1 A2 A3 "<< A1<<" " << A2<< " " << A3<< endl;
    sum = 0.5*(b-a)*(0.3478548*A1 + 0.6521451*A2 + 0.6521451*A3 + 0.3478548*A4);
    //cout<<"sum "<< sum<< endl;
    return sum;
}
double norm(vector <double> u_prs, vector <double> u_nxt, double pet)
{
    double Ih, Ih2, al, be, delta, I = 0, k = 0;
    double h = BB - AA;
    h = h/16;
    al = AA;
    be = BB;
    while(al < BB)
    {
    Ih = gauss3(al, al + h, u_prs, u_nxt);
    //cout << " Ih  " << Ih << endl;
    Ih2 = gauss4(al, al + h, u_prs, u_nxt);
    //cout << " Ih2  " << Ih2 << endl;
    delta = (Ih2 - Ih)/(63.0);
    if(abs(delta) > eps)
    {
        h /= 2;
        k++;
        if(k>50)
        {
            cout<< "UNINTEGRATED";
            exit(-1);
        }
    }
    else
    {
        I += Ih2;
        if(delta > (eps/64))
        {
            al = al + h;
            if(al + h > BB)
                h = BB - al;
        }
        else
        {
            al = al + h;
            h = 2*h;
            if(al + h > BB)
                h = BB - al;
        }

    }

    }
    if(pet == 0)
    cout << "Norma L2 ||u_next - u_present|| equals " << I <<endl;
    return I;
}

int main()
{
    setlocale(LC_ALL,"ru");
    int n = 2;
    vector <double> up, un;
    up = u_n_t(n);
    un = u_n_t(2*n);
    while(norm(up,un,0)>eps)
    {
        up = un;
        n *= 2;
        un = u_n_t(2*n);
    }
    vector <double> ur(2*n);
    double st;
    st = (BB-AA)/(2.0*n);
    for(int i = 0; i < 2*n; i++)
    {
        ur[i] = math_real_u(AA + i*st+(st/2.0));
    }
    cout << "Norma L2 ||u_next - u_real|| equals "<< norm(ur,un,1);
    return 0;
}
