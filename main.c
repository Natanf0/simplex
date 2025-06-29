#include<stdio.h>
#include<stdlib.h>
#include <string.h>
#include<math.h>
int max=-1; // se o PPL for de maximização, max = 1. Caso contrário, 0;
float *vetC, *vetX, *vetB, *matA;
double* output;
int nVar, nRestricoes;
float* A;


void printTableux(){
    char *str = malloc(nVar*10 + 1); memset(str, '-', nVar*10 + 6); str[nVar*10+6] = '\0';
    printf("\nTableux:\n\n    ");

    for(int i = 0; i<nVar; i++) printf("VNB%d  ", i);
    printf("   SBV\n");
    for(int i = 0; i<nRestricoes; i++){
        printf(" VB%i ", i);
        for(int j = 0; j<nVar; j++){
            printf("%.2lf  ", matA[i*nVar + j]);
        }
        printf(" | %.2lf \n", vetB[i]);
    }
    if(max) printf("%s\n Z = ", str);
    else printf("%s\n-Z = ", str);
    for(int i=0; i < nVar ; i++){
       printf("%.2lf  ", vetC[i]);
    }
    printf(" | %.2lf\n", vetC[nVar]);
    printf("\n\n");
}

int escolherVNBEntraNaBase(){
    // método que retorna o índice da coluna pivô   
    int indice = 0;
    float menorValor = 0; // o pivô deve ser um valor negativo

    for(int i = 0; i < nVar; i++){
        if(vetC[i] < menorValor){
            menorValor = vetC[i];
            indice = i;         
        }
    }
    return indice;
}

int escolheVBSaiDaBase(int colunaPivo){
    double r = 0.f, R_Min = pow(10,4);
    int indice = 0;
    for(int i = 0; i<nRestricoes; i++){
        if(matA[colunaPivo + i*nVar] > 0.f){
            r = vetB[i]/matA[colunaPivo + i*nVar];
            if(r<R_Min){
                R_Min = r;
                indice = i;
            } 
        }
    }
    return indice;
}


void simplex(){
    // escolhe VNB que deve entrar na base
    // escolher VB que vai deixar a base
    // selecionada pivo e sua linha, deve fazer o pivoteamento

    int colunaPivo = escolherVNBEntraNaBase();  
    int linhaPivo = escolheVBSaiDaBase(colunaPivo);

    vetX[linhaPivo] = 0;
    
    float pivo = matA[colunaPivo + linhaPivo*nVar];
    if(pivo != 1.f){
        // quando o pivo não é 1, é necessaŕio normalizar a linha pivô. Ou seja, dividir a linha do tableux pelo pivô
        for(int i = 0; i<nVar; i++){
            matA[linhaPivo*nVar + i] = matA[linhaPivo*nVar + i]/pivo;
        }
        vetB[linhaPivo] = vetB[linhaPivo]/pivo; 
    }

    // pivoteamento na Z-linha:
    float fator = vetC[colunaPivo]*(-1);
    for(int i = 0; i<nVar; i++){
        vetC[i] = matA[linhaPivo*nVar + i]*fator + vetC[i];
    }
    vetC[nVar] = fator*vetB[linhaPivo] + vetC[nVar]; // SBV ( decidi armazenar aqui ).

    // pivoteamento na matriz das restições:
    for(int i = 0; i<nRestricoes; i++){
        if(i!=linhaPivo){
            float fator = matA[colunaPivo + i*nVar]*(-1);
            for(int j = 0; j<nVar; j++){
                matA[i*nVar + j] = matA[linhaPivo*nVar + j]*fator + matA[i*nVar + j];
            }
            vetB[i] = vetB[linhaPivo]*fator + vetB[i];
        }  
    }
    printTableux();
}

int verificaOtimalidade(){

        for(int i = 0; i< nVar; i++){
            if(vetC[i]<0) return 0; // se há algum coef. negativo na função objetivo, então não estamos no ótimo
        }
        return 1; // se todos são positivos, então estamos no ótimo
    
}

int main(){
    // objetivo: maximizar / minimizar uma PPL na forma: c^t*x 
    //                                        sujeito a: A*x = b
    // vetor  C: coeficientes da funçao objetivo.
    // vetor  X: variáveis x.                                      
    // matriz A: matriz dos coeficientes das restrições
    // vetor  B: constantes das restrições
    int numIteracoes = 0;
    char tipo[12];

    printf("O problema é de maximização (max) ou minimização (min) ?:\n");
    scanf("%s", tipo);
    max = (strcmp(tipo, "max") == 0) ? 1 : 0;

    printf("Digite a quantidade de variáveis distintas:\n");
    scanf("%d", &nVar);

    printf("Digite a quantidade de restrições (exceto não-negatividade)\n");
    scanf("%d", &nRestricoes);

    output = malloc(sizeof(double)*nVar);

    vetC = (float*) malloc(sizeof(float)*(nVar+1)); // resultado será armazenado aqui
    vetX = (float*) malloc(sizeof(float)*nVar);
    vetB = (float*) malloc(sizeof(float)*nRestricoes);
    matA = (float*) malloc(sizeof(float)*nRestricoes*nVar) ;

    printf("\nDigite os coeficientes da funcao objetivo:\n");
    for(int i=0; i < nVar ; i++){
        scanf("%f", &(vetC[i]));
        if(max) vetC[i] = vetC[i]*(-1); // se for um PPL de maximização, troca os sinais dos coef. da função objetivp
    }
    vetC[nVar]=0; // SBVI

    printf("\nDigite os coeficientes da matriz de restrições:\n");
    for(int i=0; i < nRestricoes ; i++){
        for(int j=0; j< nVar ; j++){
            scanf("%f", &(matA[i*nVar+j]));
        }
    }

    printf("\nDigite as constantes associadas às restrições:\n");
    for(int i=0; i < nRestricoes ; i++){
        scanf("%f", &(vetB[i]));
    }
    int offset = nVar - nRestricoes; // O número de zeros no início
    for (int i = 0; i < nVar; i++) {
        if (i < offset)  vetX[i] = 0;
        else vetX[i] = vetB[i - offset];
    }

    printf("Tableux simplex inicial: ");
    printTableux();

    while(!verificaOtimalidade()){
        simplex();
        numIteracoes++;
    }

    printf("Solução ótima: x* = ");
    for(int i=0; i<nVar; i++){
        printf("%.2lf ", vetX[i]);
    }

    printf("\nNúmero de iterações do algorimo para o PPL : %d iterações\n\n", numIteracoes);

    free(vetX);
    free(vetC);
    free(vetB);
    free(matA);
    return 0;
}