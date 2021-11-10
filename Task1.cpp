#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
using namespace std;
const double eps = 0.001;
const int N = 2;//размерность
int ind = 0;//номер последовательности

void print(vector <double> x)
{
    for (int i = 0; i < x.size(); i++)
    {
        cout<< x[i] << "   ";
    }
    cout<<endl;
}
void print(vector <vector<double>> x)
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
double sprod(vector <double> a, vector <double> b)//скалярное произведение
{
    double ab=0;
    for(int i = 0; i < a.size(); i++)
        ab += a[i] * b[i];
    return ab;

}

vector<double> smult(double a, vector <double> b)//умножение вектора на скаляр
{
    vector <double> ab(b.size());
    for(int i = 0; i < b.size(); i++)
        ab[i] = a*b[i];
    return ab;

}

vector<double> sum(vector <double> a, vector <double> b)//сумма векторов
{
    vector <double> ab(b.size());
    for(int i = 0; i < a.size(); i++)
        ab[i] = a[i] + b[i];
    return ab;

}

vector<vector<double>> invers(vector<vector<double>> A)// поиск обратной матрицы методом Гаусса
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
            while(A[k][i] < pow(eps,3)&& k < N)
            {
                k++;
                if(k == N)//проверяю на равенство определителя нулю
                {
                 cout<<"Err. Nulevoi determinant"
                 exit(0);
                }

            }
            temp = A[i];
            A[i] = A[k];
            A[k] = temp;
        }
        p = i + 1;
        while(p < N)
        {//print(A[1]);
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

double math_func (vector <double> x)
{
    //здесь вбивается функция
    return pow(x[0],2)+pow(x[1],4);
}

vector <double> gradient(vector <double> pnt){//градиент
    vector <double> pnt_tmp1 = pnt;
    vector <double> pnt_tmp2 = pnt;
    vector <double> grad(N);
    double tau = 0.1*sqrt(eps);
    for (int i = 0; i < N; i++)
    {
        pnt_tmp1[i] += tau;
        pnt_tmp2[i] -= tau;
        grad[i] = (math_func(pnt_tmp1)-math_func(pnt_tmp2))/(2*tau);
        pnt_tmp1 = pnt;
        pnt_tmp2 = pnt;
    }
    return grad;
}

vector < vector <double>> hessian_matr(vector <double> pnt){// матрица вторых производных
    vector <double> pnt_tmp11 = pnt;
    vector <double> pnt_tmp12 = pnt;
    vector <double> pnt_tmp21 = pnt;
    vector <double> pnt_tmp22 = pnt;
    vector < vector <double>> hessian(N,vector <double>(N));

    double tau = 0.1*sqrt(eps);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            pnt_tmp11[i]+=tau;
            pnt_tmp11[j]+=tau;
            pnt_tmp12[i]+=tau;
            pnt_tmp12[j]-=tau;
            pnt_tmp21[i]-=tau;
            pnt_tmp21[j]+=tau;
            pnt_tmp22[i]-=tau;
            pnt_tmp22[j]-=tau;
            hessian[i][j]=((math_func(pnt_tmp11)-math_func(pnt_tmp21)
                        -math_func(pnt_tmp12)+math_func(pnt_tmp22))/(4*tau*tau));
            pnt_tmp11=pnt_tmp12=pnt_tmp21=pnt_tmp22=pnt;
        }
    }
    print(hessian);
    return hessian;
}

bool criteria1(vector <double> Xcur, vector <double> Xprv)
{
   ::ind++;
   vector <double> tmpgr = gradient(Xcur);
   if(
      sprod(sum(Xcur,smult(-1,Xprv)),sum(Xcur,smult(-1,Xprv))) < sqrt(eps) &&
      abs(math_func(Xcur)-math_func(Xprv) < sqrt(eps)) &&
      sprod(tmpgr,tmpgr) < sqrt(eps)
      )
        return true;
   else
        return false;
}

bool criteria2(vector <double> Xcur, vector <double> Xprv)
{
   ::ind++;
   vector <double> tmpgr = gradient(Xcur);
   if(
      sprod(sum(Xcur,smult(-1,Xprv)),sum(Xcur,smult(-1,Xprv))) < eps &&
      abs(math_func(Xcur)-math_func(Xprv)) < eps &&
      sprod(tmpgr,tmpgr) < eps
      )
        return true;
   else
        return false;
}

vector <double> fst_method (vector <double> Xcur)//выбор альфа_к методом золотого сечения
{
    double l, a = 0 , b = 1 , r , c , d;
    r = (3 - sqrt(5))/2;
    l = a - b;
    vector <double> tgrad = gradient(Xcur);
    double tempeps = sqrt(eps);
    while(l > tempeps)
    {
        c = a + r*(b - a);
        d = b - r*(b - a);
        if( math_func(sum(Xcur,smult(-c,tgrad))) < math_func(sum(Xcur,smult(-d,tgrad))))
            b = d;
        else
            a = c;
        l = a - b;
    }
    return sum(Xcur,smult(-(b-a)/2,tgrad));
}

vector <double> sec_method (vector <double> Xcur)//выбор альфа_к методом дробления шага
{
    double lam = 1/2, mu = 2, a = 1;
    vector <double> tgrad = gradient(Xcur);
    vector <double> Hk(N);
    for(int i  = 0; i<N; i++)
    {
        Hk[i] = -sprod(invers(hessian_matr(Xcur))[i],tgrad);
    }
    while(math_func(sum(Xcur,smult(a,Hk))) > math_func(Xcur))
    {
        a *= lam;
    }
    while(math_func(sum(Xcur,smult (a,Hk))) > math_func(sum(Xcur,smult(a*mu,Hk))))
    {
        a *= mu;
    }
    return sum(Xcur,smult(a,Hk));
}

void finding_min(vector <double> Xapr)
{
    vector <double> Xnxt;
    Xnxt = fst_method(Xapr);
    while(!criteria1(Xnxt, Xapr))
    {
        Xapr = Xnxt;
        Xnxt = fst_method(Xnxt);
    }
    cout<<"Pervii shag zavershilsya za "<< ind<<" iteracii"<<endl;
    cout<<"V tochke "; print(Xnxt); cout<<endl;
    cout<<"Znachenie funkcii "<< math_func(Xnxt)<<endl;
    double perv = ind;
    Xapr = Xnxt;
    Xnxt = sec_method(Xapr);
    while(!criteria2(Xnxt, Xapr))
    {
        Xapr = Xnxt;
        Xnxt = sec_method(Xnxt);
    }
    cout<<"Vtoroi shag zavershilsya za "<< ind - perv<<" iteracii"<<endl;
    cout<<"V tochke "; print(Xnxt); cout<<endl;
    cout<<"Znachenie funkcii "<< math_func(Xnxt)<<endl;
}

int main()
{
    setlocale(LC_ALL,"ru");
    vector<double> x = {1,2};//начальное приблежение
    finding_min(x);
    return 0;
}
