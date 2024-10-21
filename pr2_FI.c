#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define N 512

// Matrius i vectors
float Mat[N][N];        // Matriu de N x N
float MatDD[N][N];     // Matriu diagonal dominant
float V1[N];           // Vector 1
float V2[N];           // Vector 2
float V3[N];           // Vector 3
float V4[N];           // Vector 4

void InitData() {
    int i, j;
    srand(334411);  // Inicialització de la llavor
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++) {
            Mat[i][j] = (((i * j) % 3) ? -1 : 1) * (100.0 * (rand() / (1.0 * RAND_MAX)));
            if ((abs(i - j) <= 3) && (i != j))
                MatDD[i][j] = (((i * j) % 3) ? -1 : 1) * (rand() / (1.0 * RAND_MAX));
            else if (i == j)
                MatDD[i][j] = (((i * j) % 3) ? -1 : 1) * (10000.0 * (rand() / (1.0 * RAND_MAX)));
            else 
                MatDD[i][j] = 0.0;
        }
    for (i = 0; i < N; i++) {
        V1[i] = (i < N / 2) ? (((i * j) % 3) ? -1 : 1) * (100.0 * (rand() / (1.0 * RAND_MAX))) : 0.0;
        V2[i] = (i >= N / 2) ? (((i * j) % 3) ? -1 : 1) * (100.0 * (rand() / (1.0 * RAND_MAX))) : 0.0;
        V3[i] = (((i * j) % 5) ? -1 : 1) * (100.0 * (rand() / (1.0 * RAND_MAX)));
    }
}

void PrintVect( float vect[N], int from, int numel ){
// Validem que 'from + numel' no excedeixi la mida del vector
    if (from + numel > N) {
        printf("Error: El rang excedeix la mida del vector.\n");
        return;
    }

    // Imprimim els elements del vector
    for (int i = from; i < from + numel; i++) {
        printf("%f ", vect[i]);
    }
}

void PrintRow( float mat[N][N], int row, int from, int numel ){
    if (from + numel > N) {
        printf("Error: El rang excedeix la mida del vector.\n");
        return;
    }

    for(int i=from;i<from+numel;i++){
        print("%f ", mat[row][i]);
    }
}

void MultEscalar( float vect[N], float vectres[N], float alfa ){
    for (int i=0;i<N;i++){
        vectres[i]=vect[i]*alfa;
    }
}

float Scalar( float vect1[N], float vect2[N] ){
    float scalar=0.0;
    for (int i=0;i<N;i++){
        scalar+=(vect1[i]*vect2[i]);
    }
    return scalar;
}

float Magnitude( float vect[N] ){
    return sqrt(Scalar(vect, vect));
}

int Ortogonal( float vect1[N], float vect2[N] ){
    for(int i=0; i<N;i++){
        if(Scalar(vect1, vect2)==0){
            return 1;
        }
        else{
            return 0;
        }
    }
}

void Projection( float vect1[N], float vect2[N], float vectres[N] ){
    float numerador=Scalar(vect1, vect2);
    float denominador=Magnitude(vect2);
    for (int i=0;i<N;i++){
        vectres[i]=(numerador/denominador)*vect2[i];
    }
}

float Infininorm(float M[N][N]) {
    float maxSum = 0.0;  
    for (int i = 0; i < N; i++) {
        float rowSum = 0.0;  // Suma dels valors absoluts de la fila
        for (int j = 0; j < N; j++) {
            rowSum += fabs(M[i][j]);
        }
        if (rowSum > maxSum) {
            maxSum = rowSum;
        }
    }
    return maxSum;  // Retorna la infini-norma
}

float Onenorm( float M[N][N] ){
    float maxSum = 0.0;  
    for (int i = 0; i < N; i++) { 
        float colSum = 0.0; 
        for (int j = 0; j < N; j++) {
            colSum += fabs(M[j][i]);
        }
        if (colSum > maxSum) {
            maxSum = colSum;
        }
    }
    return maxSum;  
}

float NormFrobenius( float M[N][N] ){
    float sum=0.0;
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            sum+=M[i][j]*M[i][j];
        }
    }
    return sqrt(sum);
}

int DiagonalDom( float M[N][N] ){
    float sumDiagonal=0.0;
    float sumValors=0.0;

    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            if(i==j){
                sumDiagonal+=M[i][j];
            }
            else{
                sumValors+=fabs(M[i][j]);
            }
        }
    }
    if(fabs(sumDiagonal)>sumValors){
        return 1;
    }
    else{
        return 0;
    }
}

void Matriu_x_Vector( float M[N][N], float vect[N], float vectres[N] ){
    for (int i = 0; i < N; i++) {
        vectres[i] = 0;  // Inicialitzar a 0 abans de sumar
    }
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            vectres[i]+=M[j][i]*vect[i];
        }
    }
}

int Jacobi( float M[N][N] , float vect[N], float vectres[N], unsigned iter ){

}

int main(){
    InitData();  // Inicialitzar dades

}