//
//  main.cpp
//  main
//
//  Created by Anastasiya Moskovchenko on 01.11.16.
//  Copyright © 2016 Anastasiya Moskovchenko. All rights reserved.
//

#include <iostream>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include <cstdio>
#include <cmath>


#define eps 0.0001
#define x0 10.0
#define x1 50.0

using namespace std;

double my_function(vector<double> &x){
    return x[0] + x[1];
}

double barrier_function(vector<double> &x){
    
    return 1.0 /(2.0 * x[0]) + 1.0 /(2.0 * x[1]) ;
    
}

double new_function(vector<double> &x, int k){
    return my_function(x) + (1.0/(k * k)) * barrier_function(x);
}

vector<double> func_grad(vector<double> &x, int k){
    
    vector<double> grad(x.size());
    
    for(int i = 0; i < x.size(); ++i) {
        grad[i] = 1 - (1.0/(k * k)) * (1.0/pow(x[i], 2.0));
        
    }
    
    return grad;
}

vector<double> vector_alpha(vector<double> &x, double alpha){
    vector<double> result(x.size());
    for (int i = 0; i < x.size(); ++i)
        result[i] = x[i] * alpha;
    return result;
}

vector<double> vector_minus(vector<double> &x, vector<double> &y){
    vector<double> result(x.size());
    for (int i = 0; i < x.size(); ++i)
        result[i] = x[i] - y[i];
    return result;
}

double vector_scalar(vector<double> &x, vector<double> &y){
    double pr = 0;
    for (int i = 0; i < x.size(); ++i){
        pr += x[i] * y[i];
    }
    return pr;
}

vector<vector<double> > vector_multipl(vector<double> &x, vector<double> &y){
    vector<vector<double> > M(x.size(), vector<double>(x.size()));
    for (int i = 0; i < x.size(); ++i){
        for (int j = 0; j < x.size(); ++j){
            M[i][j] = x[i] * y[j];
        }
    }
    return M;
}

vector<vector<double> > divide_matrix(vector<vector<double> > &M, double alpha){
    vector<vector<double> > result(M.size(), vector<double>(M.size()));
    for (int i  = 0; i < M.size(); ++i){
        for (int j = 0; j < M.size(); ++j){
            result[i][j] = M[i][j] / alpha;
        }
    }
    return result;
}

vector<vector<double> > sub_matrix(vector<vector<double> > &M1, vector<vector<double> > &M2){
    vector<vector<double> > result(M1.size(), vector<double>(M1.size()));
    for (int i = 0; i < M1.size(); ++i)
        for (int j = 0; j < M1.size(); ++j)
            result[i][j] = M1[i][j] - M2[i][j];
    return result;
}

vector<vector<double> > sum_matrix(vector<vector<double> > &M1, vector<vector<double> > &M2){
    vector<vector<double> > result(M1.size(), vector<double>(M1.size()));
    for (int i = 0; i < M1.size(); ++i)
        for (int j = 0; j < M1.size(); ++j)
            result[i][j]  = M1[i][j] + M2[i][j];
    return result;
}
vector<double> matrix_vector(vector<vector<double> > &M, vector<double> &x){
    vector<double> result(x.size(),0);
    for (int i = 0; i < x.size(); ++i){
        for (int j = 0; j< x.size(); ++j){
            result[i] += M[i][j] * x[i];
        }
    }
    return result;
}


vector<vector<double> > copy_matrix(vector<vector<double> > &M){
    vector<vector<double> > result(M.size(), vector<double>(M.size()));
    
    for (int i = 0; i < M.size(); ++i){
        
        result[i] = M[i];
    }
    return result;
}

double find_min(vector<double> &x, vector<double> &d, double a, double b, int count){
    
    double e = 0;
    double x_1 = 0;
    double x_2 = 0;
    vector<double> V_a1;
    vector<double> V_a2;
    
    while(fabs(b - a) >= eps) {
        e = (b - a) * 0.001; // 1E-2 * 1E-3;
        x_1 = (b + a) / 2.0 - e;
        x_2 = (b + a) / 2.0 + e;
        
        V_a1 = vector_alpha(d, x_1);
        V_a2 = vector_alpha(d, x_2);
        
        V_a1 = vector_minus(x, V_a1);
        V_a2 = vector_minus(x, V_a2);
        
        double f1 = new_function(V_a1, count);
        double f2 = new_function(V_a2, count);
        
        if (f1 > f2){
            a = x_1;
        }else{
            b = x_2;
        }
    }
    
    return (a + b) / 2.0;
    
}

int main(){
    
    
    double alpha;

    int count = 0;
    double a = 0;
    double b = 1;
    double vec_scalar = 0;
    
    double e = 0;
    double x_1 = 0;
    double x_2 = 0;
    vector<double> V_a1;
    vector<double> V_a2;
    
    vector<double> gradient;
    vector<vector<double> > matrix;
    //vector<vector<double> > tmp;
    const double arr0[] = {0,0};
    vector<double> x;
    x.assign(arr0, arr0 + sizeof arr0 / sizeof * arr0);
    
    vector<double> y;
    y.assign(arr0, arr0 + sizeof arr0 / sizeof * arr0);
    
    vector<double> d;
    d.assign(arr0, arr0 + sizeof arr0  / sizeof * arr0);
    
    double norm;
   
    
    vector<double> old_x = x;
    vector<double> old_grad = gradient;
    
    vector<double> delta_x(x.size());
    vector<double> delta_d(x.size());
    vector<double> m_v(x.size());
    
    // initial estimate x0
    const double arr1[] = {x0,x1};
    vector<double> values;
    values.assign(arr1, arr1 + sizeof arr1 / sizeof * arr1);
    

    // select a positive definite matrix S0
    vector<vector<double> > vec_mult(values.size(), vector<double>(values.size()));
    vector<vector<double> > S(values.size(), vector<double>(values.size()));
    vector<vector<double> > dev_S(values.size(), vector<double>(values.size()));
    vector<vector<double> > tmp(values.size(), vector<double>(values.size()));
    

    
    vector<double> X0;
    
    do{
        ++count;
        x = values; /// почему x затирается при этом был size = 2 стал size = 1
        // а потом идут операции с не корректной размерностью
        X0 = values;

        
        //
        for (int i = 0; i < values.size(); ++i){
            for (int j = 0; j < values.size(); ++j){
                if (i == j){
                    S[i][j] = 1.0;
                }else{
                    S[i][j] = 0.0;
                }
            }
        }
        
        
        do{
            
            
            
            a = 0;
            b = 1;
            
            tmp = copy_matrix(S);
     
            // compulate gradient
            gradient = func_grad(x, count);
            
            old_x = x;
            old_grad = gradient;
            
            // compulate d
            d = matrix_vector(S, gradient);
            // looking for the minimum of F(x - alpha * d) using the method of one-dimensional optimization
            
            
            do
                {
                    
                    e = (b - a) * 0.000001;
                    x_1 = (a + b) / 2.0 - e;
                    x_2 = (a + b) / 2.0 + e;
                    V_a1 = vector_alpha(d, x_1);
                    V_a2 = vector_alpha(d, x_2);
                    
                    V_a1 = vector_minus(x, V_a1);
                    V_a2 = vector_minus(x, V_a2);
                    double f1 = new_function(V_a1, count);
                    double f2 = new_function(V_a2, count);
                    
                    if (f1 > f2){
                        a = x_1;
                    }else{
                        b = x_2;
                    }
                    
                    
                } while ((b - a) > eps);
            
            if (new_function(V_a1, count) > new_function(V_a2, count)){
                alpha = x_1;
            } else {
                alpha = x_2;
            }
            
            
            y = vector_alpha(d, alpha);
            x = vector_minus(x, y);
            
            gradient = func_grad(x, count);
            
            delta_x = vector_minus(x, old_x);
            delta_d = vector_minus(gradient, old_grad);
            
            m_v = matrix_vector(tmp, delta_d);
            
            vec_mult = vector_multipl(delta_x, delta_x);
            vec_scalar = vector_scalar(delta_x, delta_d);
            dev_S = divide_matrix(vec_mult, vec_scalar);
            
            S = sum_matrix(tmp, dev_S);
            
            vec_mult = vector_multipl(m_v, m_v);
            vec_scalar = vector_scalar(m_v, delta_d);
            dev_S = divide_matrix(vec_mult, vec_scalar);
            
            S = sub_matrix(S, dev_S);
            
            norm = sqrt(gradient[0] * gradient[0] + gradient[1] * gradient[1]);
            
            
        }while(fabs(new_function(x, count) - new_function(old_x, count)) > eps && count < 1E6 && norm >= eps);
        
        
        values = x;
        x.clear();
        
        
    } while (fabs(my_function(X0) - my_function(values)) > eps);
    
    
    cout<< "Значение точки минимума:\n";
    cout<< "x_0 = "<< values[0]<<endl;
    cout<< "x_1 = "<< values[1]<<endl;
    cout<< "Значение функции в точке минимума:\n";
    cout<< "F(x) = "<< my_function(values)<<endl;
    return 0;
}
