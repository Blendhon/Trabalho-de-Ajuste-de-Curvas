#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef struct Amostra {
    double t;
    double N;
} Amostra;

int contarLinhas(const char *nomeArquivo) {
    FILE *arquivo = fopen(nomeArquivo, "r");
    if (arquivo == NULL) {
        printf("Erro ao abrir o arquivo.\n");
        return -1;
    }

    int linhas = 0;
    char c;

    while ((c = fgetc(arquivo)) != EOF) {
        if (c == '\n') {
            linhas++;
        }
    }

    fclose(arquivo);

    return linhas + 1; // Contando a última linha, caso o arquivo não termine com '\n'
}

// Função para calcular os coeficientes do modelo linear
void ajuste_linear(Amostra* dados, int n, double* beta0, double* beta1) {
    double soma_t = 0, soma_N = 0, soma_tN = 0, soma_t2 = 0;

    // Calcular as somas necessárias
    for (int i = 0; i < n; i++) {
        soma_t += dados[i].t;					// Soma dos valores de t
        soma_N += dados[i].N;					// Soma dos valores de N
        soma_tN += dados[i].t * dados[i].N;		// Soma dos produtos t * N
        soma_t2 += dados[i].t * dados[i].t;		// Soma dos quadrados de t
    }

    // Calcular beta1 (coeficiente angular)
    *beta1 = (n * soma_tN - soma_t * soma_N) / (n * soma_t2 - soma_t * soma_t);

    // Calcular beta0 (coeficiente linear)
    *beta0 = (soma_N - (*beta1) * soma_t) / n;
}

// Função para calcular os coeficientes do modelo quadrático
void ajuste_quadratico(Amostra* dados, int n, double* beta0, double* beta1, double* beta2) {
    double soma_t = 0, soma_t2 = 0, soma_t3 = 0, soma_t4 = 0;
    double soma_N = 0, soma_tN = 0, soma_t2N = 0;

    // Calcular as somas necessárias
    for (int i = 0; i < n; i++) {
        double t = dados[i].t;
        double N = dados[i].N;
        double t2 = t * t;
        double t3 = t2 * t;
        double t4 = t3 * t;

        soma_t += t;			// Soma de t
        soma_t2 += t2;			// Soma de t^2
        soma_t3 += t3;			// Soma de t^3
        soma_t4 += t4;			// Soma de t^4
        soma_N += N;			// Soma de N
        soma_tN += t * N;		// Soma de t * N
        soma_t2N += t2 * N;		// Soma de t^2 * N
    }

    // Montar o sistema de equações lineares
    double A[3][4] = {
        {n, soma_t, soma_t2, soma_N},			// Primeira equação
        {soma_t, soma_t2, soma_t3, soma_tN},	// Segunda equação
        {soma_t2, soma_t3, soma_t4, soma_t2N}	// Terceira equação
    };

    // Resolver o sistema usando eliminação de Gauss
    for (int i = 0; i < 3; i++) {
        // Tornar o pivô igual a 1
        double piv = A[i][i];
        for (int j = 0; j < 4; j++) {
            A[i][j] /= piv;
        }

        // Tornar os elementos abaixo do pivô iguais a 0
        for (int k = i + 1; k < 3; k++) {
            double fator = A[k][i];
            for (int j = 0; j < 4; j++) {
                A[k][j] -= fator * A[i][j];
            }
        }
    }

    // Substituição regressiva para encontrar os valores de beta2, beta1, beta0
    for (int i = 2; i >= 0; i--) {
        A[i][3] -= (i < 2 ? A[i][2] * *beta2 : 0);
        A[i][3] -= (i < 1 ? A[i][1] * *beta1 : 0);
        if (i == 2) *beta2 = A[i][3];
        else if (i == 1) *beta1 = A[i][3];
        else *beta0 = A[i][3];
    }
}

// Função para calcular os coeficientes do modelo exponencial
void ajuste_exponencial(Amostra* dados, int n, double* beta0, double* beta1) {
    double soma_t = 0, soma_lnN = 0, soma_t_lnN = 0, soma_t2 = 0;

    // Calcular as somas necessárias
    for (int i = 0; i < n; i++) {
        double t = dados[i].t;
        double lnN = log(dados[i].N);	// Aplicar logaritmo natural em N

        soma_t += t;					// Soma de t
        soma_lnN += lnN;				// Soma de ln(N)
        soma_t_lnN += t * lnN;			// Soma de t * ln(N)
        soma_t2 += t * t;				// Soma de t^2
    }

    // Calcular beta1 (coeficiente angular) usando fórmula da regressão linear
    *beta1 = (n * soma_t_lnN - soma_t * soma_lnN) / (n * soma_t2 - soma_t * soma_t);

    // Calcular alpha0 (coeficiente linear)
    double alpha0 = (soma_lnN - (*beta1) * soma_t) / n;

    // Converter alpha0 em beta0
    *beta0 = exp(alpha0);
}

// Função para calcular o coeficiente de determinação r^2
double calcular_r2(Amostra* dados, int n, double(*modelo)(double t, double beta0, double beta1, double beta2), double* betas) {
    double soma_total = 0, soma_residual = 0, media = 0;
    
    // Calcular o valor médio de N
    for (int i = 0; i < n; i++) {
        media += dados[i].N;
    }
    media /= n;
    
    // Calcular a soma total dos quadrados
    for (int i = 0; i < n; i++) {
        soma_total += pow(dados[i].N - media, 2);
    }

    // Calcular a soma residual dos quadrados
    for (int i = 0; i < n; i++) {
        soma_residual += pow(dados[i].N - modelo(dados[i].t, betas[0], betas[1], betas[2]), 2);
    }

    return 1 - (soma_residual / soma_total);
}

// Função que representa o modelo linear
double modelo_linear(double t, double beta0, double beta1, double beta2) {
    return beta0 + beta1 * t;
}

// Função que representa o modelo quadrático
double modelo_quadratico(double t, double beta0, double beta1, double beta2) {
    return beta0 + beta1 * t + beta2 * t * t;
}

// Função que representa o modelo exponencial
double modelo_exponencial(double t, double beta0, double beta1, double beta2) {
    return beta0 * exp(beta1 * t);
}

int main() {
    const char *nomeArquivo = "dados.txt";	// Nome do arquivo com dados
    int numLinhas = contarLinhas(nomeArquivo);
    
    if (numLinhas <= 0) {
        printf("O arquivo nao contem dados validos.\n");
        return 1;
    }

    Amostra* vetor = (Amostra*)malloc(numLinhas * sizeof(Amostra));
    if (vetor == NULL) {
        printf("Erro ao alocar memória.\n");
        return 1;
    }

    FILE *arquivo = fopen(nomeArquivo, "r");
    if (arquivo == NULL) {
        printf("Erro ao abrir o arquivo.\n");
        free(vetor);
        return 1;
    }

    for (int i = 0; i < numLinhas; i++) {
        if (fscanf(arquivo, "%lf %lf\n", &vetor[i].t, &vetor[i].N) != 2) {
            printf("Erro ao ler os dados na linha %d.\n", i + 1);
            free(vetor);
            fclose(arquivo);
            return 1;
        }
    }
    fclose(arquivo);

	printf("--------------------------------\n\n"
        "     ALGORITIMOS NUMERICOS!\n"
        "        AJUSTE DE CURVAS\n"
        "        ALUNOS: BLENDHON\n"
        "\n--------------------------------\n\n");
	double beta0, beta1, beta2;

    // Ajuste do modelo linear
    ajuste_linear(vetor, numLinhas, &beta0, &beta1);
    printf("  ------ Modelo  Linear ------\n\n");
    printf("  beta0 = %.8lf\n", beta0);
    printf("  beta1 = %.8lf\n", beta1);

    // Calcular R² para o modelo linear
    double betas_linear[3] = {beta0, beta1, 0}; // Array nomeado
    double r2_linear = calcular_r2(vetor, numLinhas, modelo_linear, betas_linear);
    printf("  R^2 Linear: %.8lf\n\n", r2_linear);

    // Ajuste do modelo quadrático
    ajuste_quadratico(vetor, numLinhas, &beta0, &beta1, &beta2);
    printf("  ---- Modelo  Quadratico ----\n\n");
    printf("  beta0 = %.8lf\n", beta0);
    printf("  beta1 = %.8lf\n", beta1);
    printf("  beta2 = %.8lf\n", beta2);

    // Calcular R² para o modelo quadrático
    double betas_quadratico[3] = {beta0, beta1, beta2}; // Array nomeado
    double r2_quadratico = calcular_r2(vetor, numLinhas, modelo_quadratico, betas_quadratico);
    printf("  R^2 Quadratico: %.8lf\n\n", r2_quadratico);

    // Ajuste do modelo exponencial
    ajuste_exponencial(vetor, numLinhas, &beta0, &beta1);
    printf("  ---- Modelo Exponencial ----\n\n");
    printf("  beta0 = %.8lf\n", beta0);
    printf("  beta1 = %.8lf\n", beta1);
    
    // Calcular R² para o modelo exponencial
    double betas_exponencial[3] = {beta0, beta1, 0}; // Array nomeado
    double r2_exponencial = calcular_r2(vetor, numLinhas, modelo_exponencial, betas_exponencial);
    printf("  R^2 Exponencial: %.8lf\n\n", r2_exponencial);
    
    free(vetor);
    
    printf("--------------------------------\n\n"
    	"\t Melhor Metodo:\n");
    
    if (r2_linear > r2_quadratico && r2_linear > r2_exponencial) {
    	printf("\tRegressao Linear\n"
			"\t %.12lf\n", r2_linear);
	} else if (r2_quadratico > r2_linear && r2_quadratico > r2_exponencial) {
		printf("      Regressao Quadratica\n"
			"\t %.12lf\n", r2_quadratico);
	} else {
		printf("      Regressao Exponencial\n"
			"\t %.12lf\n", r2_exponencial);
	}
    
    return 0;
}

