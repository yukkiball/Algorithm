#include <iostream>
#include <vector>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;
void Gauss(int n);  //˳��Gauss��Ԫ��
void main_Gauss(int n);  //����ԪGauss��Ԫ��
void Jacobi(int n);  //Jacobi������
void Gauss_Seidel(int n); //Gauss_Seidel������
const int l = 10;

int main()
{
    Gauss(l);
    main_Gauss(l);
    Jacobi(l);
    Gauss_Seidel(l);
    return 0;
}

void Gauss(int n)
{
    //��������[A b]
    double a[n][n + 1] = {0};
    for (int k = 0; k < n + 1; ++k)
        a[0][k] = 0;
    a[0][0] = 6;
    a[0][1] = 1;
    a[0][n] = 7;
    a[n - 1][n - 1] = 6;
    a[n - 1][n - 2] = 8;
    a[n - 1][n] = 14;
    for (int k = 1; k < n - 1; ++k)
    {
        a[k][k] = 6;
        a[k][k - 1] = 8;
        a[k][k + 1] = 1;
        a[k][n] = 15;
    }
    //��Ԫ
    for (int k = 1; k < n ; ++k)
    {
        a[k][k] = a[k][k] - a[k][k - 1] * a[k - 1][k]/ a[k - 1][k - 1];
        a[k][n] = a[k][n] - a[k][k - 1] * a[k - 1][n]/ a[k - 1][k - 1];
        a[k][k - 1] = 0;
    }
    //��ʾ����
    cout << "˳��Gauss��Ԫ������Ϊ��" << endl;
    for (int k = 0; k < n; ++k)
    {
        for (int j = 0; j < n + 1; ++j)
            cout << a[k][j] << " ";
        cout << endl;
    }
    //���
    double x[n] = {0};
    x[n - 1] = a[n - 1][n]/a[n - 1][n - 1];
    for (int k = n - 2; k > -1; --k)
        x[k] = (a[k][n] - x[k + 1] * a[k][k + 1])/a[k][k];
    //��ʾ��
    cout << endl;
    cout << "˳��Gauss��Ԫ����Ϊ��" << endl;
    for (int k = 0; k < n; ++ k)
        cout << x[k] << endl;
    cout << endl;
}

void main_Gauss(int n)
{
    //��������[A b]���������б任����ÿһ������һ�У���һ���ƶ������һ�У�
    double a[n][n + 1] = {0};
    for (int k = 0; k < n + 1; ++k)
        a[0][k] = 0;
    for (int k = 0; k < n - 2 ; ++k)
    {
        a[k][k] = 8;
        a[k][k + 1] = 6;
        a[k][k + 2] = 1;
        a[k][n] = 15;
    }
    a[n - 2][n - 2] = 8;
    a[n - 2][n - 1] = 6;
    a[n - 2][n] = 14;
    a[n - 1][0] = 6;
    a[n - 1][1] = 1;
    a[n - 1][n] = 7;
    double m = 6;
    double z = 1;
    double p = 7;
    //ÿһ����Ԫ��Ϊ8
    for (int k = 0; k < n - 2; ++k)
    {
        double tempm = m;
        double tempz = z;
        m = tempz - 0.75 * tempm;
        z = - 0.125 * tempm;
        p -= 1.875 * tempm;
    }
    double tempm = m;
    double tempz = z;
    m = tempz - 0.75 * tempm;
    z = - 0.125 * tempm;
    p -= 1.75 * tempm;
    a[n - 1][0] = 0;
    a[n - 1][1] = 0;
    a[n - 1][n - 1] = m;
    a[n - 1][n] = p;
    //��ʾ����
    cout << "����ԪGauss��Ԫ������Ϊ��" << endl;
    for (int k = 0; k < n; ++k)
    {
        for (int j = 0; j < n + 1; ++j)
            cout << a[k][j] << " ";
        cout << endl;
    }

    //���
    double x[n] = {0};
    x[n - 1] = a[n - 1][n]/a[n - 1][n - 1];
    x[n - 2] = (a[n - 2][n] - 6 * x[n - 1])/a[n - 2][n - 2];
    for (int k = n - 3; k > -1; --k)
        x[k] = (a[k][n] - x[k + 2] - 6 * x[k + 1])/a[k][k];
    //��ʾ��
    cout << endl;
    cout << "����ԪGauss��Ԫ����Ϊ��" << endl;
    for (int k = 0; k < n; ++ k)
        cout << x[k] << endl;
    cout << endl;
}

void Jacobi(int n)
{
    vector<double> x(n, 0);
    vector<double> y(n, 0);
    vector<double> z(n, 0);
    int i = 0;
    double max_e = 0;
    do
    {
        y.assign(x.begin(), x.end());
        x[0] = (7.0 - y[1])/6.0;
        for (int k = 1; k < n - 1; ++k)
            x[k] = (15 - 8.0 * y[k - 1] - y[k + 1])/6.0;
        x[n - 1] = (14.0 - 8.0 * y[n - 2])/6.0;
        for (int k = 0; k < n; ++k)
            z[k] = abs(x[k] - y[k]);
        auto max_ele = max_element(z.begin(), z.end());
        max_e = *max_ele;
        ++i;
    }while(max_e >= 0.0001); //ʹ��һ��������
    cout << "Jacobi����������������" << i << endl;
    cout << "��Ϊ��" << endl;
    for (int k = 0 ; k < n; ++k)
        cout << x[k] << endl;
    cout << endl;
}

void Gauss_Seidel(int n)
{
    vector<double> x(n, 0);
    vector<double> y(n, 0);
    vector<double> z(n, 0);
    int i = 0;
    double max_e = 0;
    do
    {
        y.assign(x.begin(), x.end());
        x[0] = (7.0 - x[1])/6.0;
        for (int k = 1; k < n - 1; ++k)
            x[k] = (15 - 8.0 * x[k - 1] - x[k + 1])/6.0;
        x[n - 1] = (14.0 - 8.0 * x[n - 2])/6.0;
        for (int k = 0; k < n; ++k)
            z[k] = abs(x[k] - y[k]);
        auto max_ele = max_element(z.begin(), z.end());
        max_e = *max_ele;
        ++i;
    }while(max_e >= 0.0001); //ʹ��һ��������
    cout << "Gauss-Seidel����������������" << i << endl;
    cout << "��Ϊ��" << endl;
    for (int k = 0 ; k < n; ++k)
        cout << x[k] << endl;
    cout << endl;
}
