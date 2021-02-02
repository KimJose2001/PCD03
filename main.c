#include "mpi.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define N 2048
#define SRAND_VALUE 1985
#define EX 2000

void MatrizVazia(int **grid, int linhas, int colunas)
{
    int i, j;
    grid = (int **)malloc(sizeof(int *) * (linhas));
    for(i = 0; i<linhas; i++)
    {
        grid[i] = (int*)malloc(sizeof(int) * colunas);
    }
}

void Populacao(int **grid)
{
    int i,j;
    int cont = 0;
    srand(SRAND_VALUE);
    for(i= 0; i<N; i++)
    {
        for (j = 0; j < N; j++)   //Percorre o Vetor de Inteiros atual.
        {
            grid[i][j] = rand() % 2;
            cont += grid[i][j];
        }
    }
    printf("Geracao Inicial: %d", cont);
}

void ImprimeMatriz(int **M, int m, int n)
{
    int i,j;
    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
            printf(" %2d ",M[i][j]);
        printf("\n");
    }
}

void LiberarMatriz(int **Matriz)
{
    int i;
    for(i = 0; i < N; i++)
        free(Matriz[i]);
    free(Matriz);
}

int numVizinhos(int **Matriz, int linha, int coluna)
{
    int vizinhos = 0;
    int x_up, y_up;
    int x_right, y_right;
    int x_left, y_left;
    int x_down, y_down;
    int x_ul, y_ul;
    int x_ur, y_ur;
    int x_dl, y_dl;
    int x_dr, y_dr;

    if(linha == 0)
    {
        x_up = N - 1;
        x_down = linha + 1;
        x_ul = N - 1;
        x_ur = N - 1;
        x_dl = linha + 1;
        x_dr = linha + 1;
    }
    else if(linha == N - 1)
    {
        x_up = linha - 1;
        x_down = 0;
        x_ul = linha - 1;
        x_ur = linha - 1;
        x_dl = 0;
        x_dr = 0;
    }
    else
    {
        x_up = linha - 1;
        x_down = linha + 1;
        x_ul = linha - 1;
        x_ur = linha - 1;
        x_dl = linha + 1;
        x_dr = linha + 1;
    }
    y_up = coluna;
    y_down = coluna;
    if(coluna == 0)
    {
        y_left = N - 1;
        y_right = coluna + 1;
        y_ul = N - 1;
        y_dl = N - 1;
        y_ur = coluna + 1;
        y_dr = coluna + 1;
    }
    else if(coluna == N - 1)
    {
        y_left = coluna - 1;
        y_right = 0;
        y_ul = coluna - 1;
        y_dl = coluna - 1;
        y_ur = 0;
        y_dr = 0;
    }
    else
    {
        y_left = coluna - 1;
        y_right = coluna + 1;
        y_ul = coluna - 1;
        y_dl = coluna - 1;
        y_ur = coluna + 1;
        y_dr = coluna + 1;
    }
    x_left = linha;
    x_right = linha;
    vizinhos = Matriz[x_up][y_up] + Matriz[x_down][y_down] + Matriz[x_left][y_left] + Matriz[x_right][y_right]
               +Matriz[x_ul][y_ul] + Matriz[x_ur][y_ur] + Matriz[x_dl][y_dl] + Matriz[x_dr][y_dr];
    return vizinhos;
}

int novoEstado(int linha, int coluna, int **grid)
{
    int i = linha;
    int j = coluna;
    if(grid[i][j] == 1)
    {
        if(numVizinhos(grid, i, j) < 2)
            return 0;
        else if( numVizinhos(grid, i, j)== 2 || numVizinhos(grid, i, j) == 3)
            return 1;
        else if(numVizinhos(grid, i, j) >= 4)
            return 0;
    }
    else
    {
        if(numVizinhos(grid, i, j) == 3)
            return 1;
        else
            return 0;
    }
}

void novaGeracao(int **grid, int **newgrid)
{
    int i, j, hcont;
    int cont=0;
    #pragma omp parallel for private(i, j, hcont) reduction(+:cont)
    for(i = 0; i < N; i++)
    {
        hcont = 0;
        for(j = 0; j < N; j++)
        {
            if(grid[i][j] == 1)
                hcont += 1;
            newgrid[i][j] = novoEstado(i, j, grid);
        }
        cont += hcont;
    }
    for(i = 0; i < N; i++)
    {
        for(j = 0; j < N; j ++)
            grid[i][j] = newgrid[i][j];
    }
}

void ImprimeVizinhos(int **Matriz)
{
    int i, j;
    for(i = 0; i < N; i++)
    {
        for(j = 0; j < N; j++)
            printf("\nNumero de vizinhos de (%d, %d): %d", i, j, numVizinhos(Matriz, i, j));
    }
}

int main(int argc, char * argv[])
{
    double start, end; //tempos
    int esteProc, numProc, sentinela = 0, proc_anterior, proc_sucessor, i, j, k, vivos, parte, iproc, tag_up = 10, tag_down = 11;
    int **grid, **newGrid;
    int *rec1, *send1, *rec2, *send2;

    for(i = 0; i < N; i++)
    {
        rec1 = (int *) malloc(sizeof(int)*N);
        send1 = (int *) malloc(sizeof(int)*N);
        rec2 = (int *) malloc(sizeof(int)*N);
        send2 = (int *) malloc(sizeof(int)*N);
    }

    MPI_Status status;
    start = clock();
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &esteProc);

    parte = N/numProc;
    vivos = 0;

    sentinela = esteProc*parte;
    if(sentinela + parte > N)
    {
        parte = N+sentinela-1;
    }

    proc_anterior = ((esteProc+numProc-1)%numProc);
    proc_sucessor = ((esteProc+1)%numProc);

    printf("Processo[%d][%d] - Tamanho %d\n", esteProc, numProc, parte);
    fflush(stdout);
    grid = (int **)malloc(sizeof(int *) * (parte+2));
    newGrid = (int **)malloc(sizeof(int *) * (parte+2));
    for(i = 0; i<parte+2; i++)
    {
        grid[i] = (int*)malloc(sizeof(int) * N+2);
        newGrid[i] = (int*)malloc(sizeof(int) * N+2);
    }


    srand(SRAND_VALUE);
    for(i = 0; i < esteProc*parte*N; i++)
    {
        rand();
    }

    for(i = 1; i<=parte; i++)
    {
        for(j = 1; j <= N; j++)
        {
            grid[i][j] = rand()%2;
        }
    }

    for(i = 0; i < EX; i++)
    {

        for(j = 0; j < N; j++)
        {
            send1[j] = grid[1][j+1];
            send2[j] = grid[parte][j+1];
        }

            //printf("MPI process %d sends value send2 to MPI process %d.\n", esteProc, proc_sucessor);
        MPI_Sendrecv (send2, N, MPI_INTEGER, proc_sucessor, tag_up,
                      rec1, N, MPI_INTEGER, proc_anterior, tag_up, MPI_COMM_WORLD, &status);
            //printf("MPI process %d received value rec1 from MPI process %d.\n", esteProc, proc_anterior);
        /*MPI_Send(send2, N, MPI_INTEGER, proc_sucessor, tag_up, MPI_COMM_WORLD);
        MPI_Recv(rec1, N, MPI_INTEGER, proc_anterior, tag_up, MPI_COMM_WORLD, &status);*/

            //printf("MPI process %d sends value send1 to MPI process %d.\n", esteProc, proc_anterior);
        MPI_Sendrecv (send1, N, MPI_INTEGER, proc_anterior, tag_down,
                      rec2,  N, MPI_INTEGER, proc_sucessor, tag_down, MPI_COMM_WORLD, &status);
            //printf("MPI process %d received value rec2 from MPI process %d.\n", esteProc, proc_sucessor);

        /*MPI_Send(send1, N, MPI_INTEGER, proc_anterior, tag_down, MPI_COMM_WORLD);
        MPI_Recv(rec2, N, MPI_INTEGER, proc_sucessor, tag_down, MPI_COMM_WORLD, &status);*/

        MPI_Barrier(MPI_COMM_WORLD);

        for(j = 0; j < N; j++)
        {
            grid[0][j+1] = rec1[j];
            grid[parte+1][j+1] = rec2[j];
        }

        for(k = 0; k <= parte+1; k++)
        {
            grid[k][N+1] = grid[k][1];
            grid[k][0] = grid[k][N];
        }

        for(k = 1; k <= parte; k++)
        {
            for(j = 1; j < N+1; j++)
            {
                newGrid[k][j] = novoEstado(k, j, grid);
            }
        }
        //////////////////////////////////
        MPI_Barrier(MPI_COMM_WORLD);
        vivos = 0;

        for(k = 1; k <= parte; k++)
        {
            for(j = 1; j <= N; j++)
            {
                vivos = vivos + newGrid[k][j];
                grid[k][j] = newGrid[k][j];
            }
        }
    }

    if(esteProc != 0)
    {/* Processos escravos */
        MPI_Send(&vivos, 1, MPI_INTEGER,0,tag_up, MPI_COMM_WORLD);
    }
    else
    {/* Processo Mestre */
        int results[numProc], total = 0;
        results[0] = vivos;

        for(iproc = 1; iproc < numProc; iproc++)
        {
            MPI_Recv(&results[iproc], 1, MPI_INTEGER,iproc,tag_up, MPI_COMM_WORLD, &status);
        }

        for(iproc = 0; iproc < numProc; iproc++)
        {
            total = total + results[iproc];
        }
        end = clock();
        printf("Celulas Vivas (Ultima Geracao): %d\n",total);

        printf("Tempo de execucao: %.3f s\n", ((double)end-start)/(double)CLOCKS_PER_SEC);

    }

    MPI_Finalize();
    free(rec1);
    free(send1);
    free(rec2);
    free(send2);

    return 0;
}
