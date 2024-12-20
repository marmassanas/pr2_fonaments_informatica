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

#1
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
    printf("\n");
}

#2
void PrintRow( float mat[N][N], int row, int from, int numel ){
    // Validem que 'from + numel' no excedeixi la mida del vector
    if (from + numel > N) {
        printf("Error: El rang excedeix la mida del vector.\n");
        return;
    }
    //imprimeix els elements seleccionats d’una fila específica d'una matriu.
    for(int i=from;i<from+numel;i++){
        printf("%f ", mat[row][i]);
    }
    printf("\n");
}

#3
void MultEscalar( float vect[N], float vectres[N], float alfa ){
    //multiplica un vector per un escalar (alfa) i ho guarda a un vector nou (vectres)
    for (int i=0;i<N;i++){
        vectres[i]=vect[i]*alfa;
    }
}

#4
float Scalar( float vect1[N], float vect2[N] ){
    //calcula el producte escalar entre dos vectors
    float scalar=0.0;
    for (int i=0;i<N;i++){
        scalar+=(vect1[i]*vect2[i]);
    }
    return scalar;
}

#5
float Magnitude( float vect[N] ){
    //calcula la magnitud d'un vector utilitzant la funció Scalar per obtenir el quadrat de la magnitud i després retorna l'arrel quadrada
    return sqrt(Scalar(vect, vect));
}

#6
int Ortogonal( float vect1[N], float vect2[N] ){
    //determina si dos vectors són ortogonals
    //Si el producte escalar és zero i per tant són ortogonals, retorna 1
    if(Scalar(vect1, vect2)==0){
        return 1;
    }
    //sino retorna 0
    else{
        return 0;
    }
}

#7
void Projection( float vect1[N], float vect2[N], float vectres[N] ){
    //calcula la projecció d'un vector (vect1) sobre un altre vector (vect2)
    //obtenir el producte escalar dels dos vectors
    float numerador=Scalar(vect1, vect2);
    //calcula la magnitud de vect2
    float denominador=Magnitude(vect2);
    MultEscalar(vect2, vectres, numerador/denominador);
}

#8
float Infininorm(float M[N][N]) {
    //calcula la infini-norma d'una matriu, que és el màxim de les sumes dels valors absoluts de cada fila
    float maxSum = 0.0;  
    for (int i = 0; i < N; i++) {
        float rowSum = 0.0;  // Suma dels valors absoluts de la fila
        for (int j = 0; j < N; j++) { //recorrem la matriu
            rowSum += fabs(M[i][j]); //sumem el valor absolut de cada casella de la fila a rowSum
        }
        if (rowSum > maxSum) { //comparem valors per obtenir el més gran
            maxSum = rowSum;
        }
    }
    return maxSum;  // Retorna la infini-norma
}

#9
float Onenorm( float M[N][N] ){
    // calcula la norma 1 d'una matriu, que és el màxim de les sumes de les columnes.
    float maxSum = 0.0;  
    for (int i = 0; i < N; i++) { 
        float colSum = 0.0; // Suma dels valors absoluts de la columna
        for (int j = 0; j < N; j++) { //recorrem la matriu verticalment
            colSum += fabs(M[j][i]); //sumem el valor absolut de cada casella de la columna a colSum
        }
        if (colSum > maxSum) { //comparem valors per obtenir el més gran
            maxSum = colSum;
        }
    }
    return maxSum;  
}

#10
float NormFrobenius( float M[N][N] ){
    //calcula la norma de Frobenius d'una matriu
    float sum=0.0;
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            sum+=pow(M[i][j],2);//sumant els quadrats de tots els elements de la matriu
        }
    }
    return sqrt(sum); //retornant de l'arrel quadrada del resultat
}

#11
int DiagonalDom( float M[N][N] ){
    //comprova si una matriu és diagonal dominant
    float sumDiagonal=0.0;
    float sumValors=0.0;
    //compara la suma dels valors absoluts fora de la diagonal amb el valor absolut de la suma dels elements diagonals
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){//recorrem la matriu
            if(i==j){
                sumDiagonal+=M[i][j];//suma dels valors absoluts dels elements diagonals
            }
            else{
                sumValors+=fabs(M[i][j]); //suma dels valors absoluts fora de la diagonal
            }
        }
    }
    if(fabs(sumDiagonal)>sumValors){
        //si la suma dels valors de la diagonal és major a la suma de la resta de valors, significa que és dominant i per tant retorna 1
        return 1;
    }
    else{
        //sino no és dominant i retorna 0
        return 0;
    }
}

#12
void Matriu_x_Vector( float M[N][N], float vect[N], float vectres[N] ){
    //realitza la multiplicació d'una matriu per un vector
    for (int i = 0; i < N; i++) {
        vectres[i] = 0;  // Inicialitzar a 0 abans de sumar del vector on es guarda el resultat
    }
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){//recorrem la matriu i el vector
            vectres[i]+=M[i][j]*vect[j]; //sumem els productes de cada element corresponent i ho guardem a vectres
        }
    }
}

#13
int Jacobi(float M[N][N], float vect[N], float vectres[N], unsigned iter) {
    // Comprovem si la matriu és diagonalment dominant
    if (!DiagonalDom(M)) {
        return 0;  // No es pot aplicar el mètode de Jacobi
    }

    float xnew[N];  // Vector temporal per a les noves solucions

    // Iteracions del mètode de Jacobi
    for (unsigned k = 0; k < iter; k++) {
        for (int i = 0; i < N; i++) {
            float sum = 0;
            for (int j = 0; j < N; j++) {
                if (i != j) {
                    sum += M[i][j] * vectres[j];
                }
            }
            xnew[i] = (vect[i] - sum) / M[i][i];
        }

        // Actualitzem el vector de resultats
        for (int i = 0; i < N; i++) {
            vectres[i] = xnew[i];
        }
    }

    return 1;  // El mètode de Jacobi s'ha aplicat correctament
}



int main(){
    InitData();  // Inicialitzar dades

    //a
    printf("V1 del 0 al 9 i del 256 al 265: \n");
    PrintVect(V1, 0, 10);
    PrintVect(V1, 256, 10);
    printf("\n");

    printf("V2 del 0 al 9 i del 256 al 265: \n");
    PrintVect(V2, 0, 10);
    PrintVect(V2, 256, 10);
    printf("\n");

    printf("V3 del 0 al 9 i del 256 al 265: \n");
    PrintVect(V3, 0, 10);
    PrintVect(V3, 256, 10);
    printf("\n");

    //b
    printf("Mat fila 0 i fila 100 del 0 al 9: \n");
    PrintRow(Mat, 0, 0, 10);
    PrintRow(Mat, 100, 0, 10);
    printf("\n");

    //c
    printf("MatDD fila 0 del 0 al 9 i fila 100 del 95 al 104: \n");
    PrintRow(MatDD, 0, 0, 10);
    PrintRow(MatDD, 100, 95, 10);
    printf("\n");
    
    //d
    printf("Infininorma de Mat = %f\n", Infininorm(Mat));
    printf("Norma ú de Mat = %f\n", Onenorm(Mat));
    printf("Norma de Frobenius de Mat = %f\n", NormFrobenius(Mat));  
    if(DiagonalDom(Mat)==0){
        printf("La matriu Mat no és diagonal dominant. \n");
    }
    else{
        printf("La matriu Mat és diagonal dominant. \n");
    }

    printf("\n");

    printf("Infininorma de MatDD = %f\n", Infininorm(MatDD));
    printf("Norma ú de MatDD = %f\n", Onenorm(MatDD));
    printf("Norma de Frobenius de MatDD = %f\n", NormFrobenius(MatDD));  
    if(DiagonalDom(MatDD)==0){
        printf("La matriu MatDD no és diagonal dominant. \n");
    }
    else{
        printf("La matriu MatDD és diagonal dominant. \n");
    }
    printf("\n");

    //e
    printf("Escalar <V1,V2> = %f Escalar <V1,V3> = %f Escalar <V2,V3> = %f\n", Scalar(V1, V2), Scalar(V1, V3), Scalar(V2, V3));
    printf("\n");

    //f
    printf("Magnitud V1,V2 i V3 = %f, %f, %f\n", Magnitude(V1), Magnitude(V2), Magnitude(V3));
    printf("\n");

    //g
    if(Ortogonal(V1, V2)==1){
        printf("V1 i V2 són ortogonals. \n");
    }
    else{
        printf("V1 i V2 no són ortogonals. \n");
    }

    if(Ortogonal(V1, V3)==1){
        printf("V1 i V3 són ortogonals. \n");
    }
    else{
        printf("V1 i V3 no són ortogonals. \n");
    }

    if(Ortogonal(V2, V3)==1){
        printf("V2 i V3 són ortogonals. \n");
    }
    else{
        printf("V2 i V3 no són ortogonals. \n");
    }
    printf("\n");

    //h
    float VH[N];
    // Multiplicació de V3 per l'escalar 2.0
    MultEscalar(V3, VH, 2.0);

    // Visualitzar els elements 0 al 9 i 256 al 265 de V3x2
    printf("Els elements 0 al 9 i 256 al 265 del resultat de multiplicar V3x2.0 són:\n");
    PrintVect(VH, 0, 10);
    PrintVect(VH, 256, 10);
    printf("\n");

    //i
    float VI1[N];
    float VI2[N];
    Projection(V2, V3, VI2);
    Projection(V1, V2, VI1);
    
    printf("Els elements 0 a 9 del resultat de la projecció de V2 sobre V3 són:\n");
    PrintVect(VI2, 0, 10);
    printf("Els elements 0 a 9 del resultat de la projecció de V1 sobre V2 són:\n");
    PrintVect(VI1, 0, 10);
    printf("\n");

    //j
    float VJ[N];
    Matriu_x_Vector(Mat, V2, VJ);
    printf("Els elements 0 a 9 del resultat de la multiplicació de Mat per v2 són:\n");
    PrintVect(VJ, 0, 10);
    printf("\n");

    //k
    float VK[N];
    Jacobi(MatDD, V3, VK, 1);
    printf("Els elements 0 a 9 de la solució (1 iter) del sistema d'equacions són:\n");
    PrintVect(VK, 0, 10);

    // Tornar a inicialitzar VK amb zeros
    for (int i = 0; i < N; i++) {
        VK[i] = 0.0;
    }

    // Calcular la solució amb 1000 iteracions
    Jacobi(MatDD, V3, VK, 1000);
    printf("Els elements 0 a 9 de la solució (1000 iters) del sistema d'equacions són:\n");
    PrintVect(VK, 0, 10);

    // Resoldre Mat * X = V3
    // Tornar a inicialitzar VK amb zeros
    for (int i = 0; i < N; i++) {
        VK[i] = 0.0;
    }

    // Comprovar si Mat és diagonal dominant
    if (DiagonalDom(Mat) == 0) {
        printf("La matriu Mat no és diagonal dominant, no es pot aplicar Jacobi.\n");
    } else {
        Jacobi(Mat, V3, VK, 1000);
        printf("Els elements 0 a 9 de la solució (1000 iters) del sistema d'equacions Mat * X = V3 són:\n");
        PrintVect(VK, 0, 10);
    }
    printf("\n");

}