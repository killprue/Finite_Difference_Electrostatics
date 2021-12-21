#include <iostream>
#include <fstream>

using namespace std;

extern "C" void  sgesv_(int *n, int *nrhs, float *a, int *lda, int *ipiv,
float *b, int *ldb, int *info );

int main() {
    int n = 30;

    float** matrix = new float*[n*n*n];
    float* RHS = new float[n*n*n];
    float* X = new float[(n*n*n)*(n*n*n)];

    for (int i=0;i<n*n*n;i++) {
        matrix[i]=&X[i*n*n*n];
    }

    for(int i = 0; i < n*n*n; i++){
        for(int j = 0; j < n*n*n; j++){
            matrix[i][j] = 0;
        }
    }

    for(int i = 0; i < n*n*n; i++){
        RHS[i] = 0;
    }

    int q = 0;
    int j = 0;
    int k = 0;
    int i = 0;

    j = 0;
    for(i = 0; i < n; i++){
        for(k = 0; k < n; k++){
            q = (i)*(n*n)+j*(n)+k;
            matrix[q][q] = 1.0;
        }
    }

    j = n-1;
    for(i = 0; i < n; i++){
        for(k = 0; k < n; k++){
            q = (i)*(n*n)+j*(n)+k;
            matrix[q][q] = 1.0;
        }
    }

    i = n-1;
    for(j = 0; j < n; j++){
        for(k = 0; k < n; k++){
            q = (i)*(n*n)+j*(n)+k;
            matrix[q][q] = 1.0;
        }
    }

    i = 0;
    for(j = 0; j < n; j++){
        for(k = 0; k < n; k++){
            q = (i)*(n*n)+j*(n)+k;
            matrix[q][q] = 1.0;
        }
    }

    k = n-1;
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            q = (i)*(n*n)+j*(n)+k;
            matrix[q][q] = 1.0;
        }
    }

    k = 0;
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            q = (i)*(n*n)+j*(n)+k;
            matrix[q][q] = 1.0;
        }
    }

    for(i = 1; i < n-1; i++){
        for(j = 1; j < n-1; j++){
            for(k = 1; k < n-1; k++){
                q = (i)*(n*n)+(j*n)+k;
                matrix[q][q-(n*n)] = 1;
                matrix[q][q-(n)] = 1.0;
                matrix[q][q-1] = 1.0;
                matrix[q][q] = -6.0;
                matrix[q][q+1] = 1.0;
                matrix[q][q+(n)] = 1.0;
                matrix[q][q+(n*n)] = 1;
            }
        }
    }
    
    float length = 10;
    float n_float = float(n);
    float h = length/n_float;
    float rho = 1;
    float ep = 1;
    float sphere_radius = 8.0;
    
    float x = 0;
    float y = 0;
    float z = 0;
    float dist = 0;
    float sphere_distance = 0;
    
    float x_center = 0;
    float y_center = 0;
    float z_center = 0;
    
    for(i = 1; i < n-1; i++){
        for(j = 1; j < n-1; j++){
            for(k = 1; k < n-1; k++){
                y = h*(i+1 - n_float/2);
                x = h*(j+1 - n_float/2);
                z = h*(k+1 - n_float/2);
                dist = sqrt(x*x + y*y + z*z);
                if(dist < sphere_radius){
                    x_center = -sphere_radius/2;
                    
                    sphere_distance = sqrt((x-x_center)*(x-x_center) + (y-y_center)*(y-y_center) + (z-z_center)*z-(z_center));
                    
                    if(sphere_distance > sphere_radius/4){
                        RHS[(i)*(n*n)+(j*n)+k] = h*h*rho/ep;
                    }
                }
            }
        }
    }

    int N = n*n*n;
    int nrhs = 1;
    int LDA = n*n*n;
    int LDB = n*n*n;
    int *IPIV;
    IPIV = new int[n*n*n];
    int INFO;

    ofstream file;
    file.open("");
    file << "RHS" << endl;
    for(i = 0; i < n*n*n; i++){
        file << RHS[i] << endl;
    }
    file.close();
    
    sgesv_(&N,&nrhs,X,&LDA,IPIV,RHS,&LDB,&INFO);
    
    file.open("");
    file << "V" << endl;
    for(i = 0; i < n*n*n; i++){
        if(RHS[i] > 0){
            RHS[i] = 0;
        }
        file << RHS[i] << endl;
    }
    file.close();
    
    delete[] matrix;
    delete[] X;
    delete[] RHS;
    
}











