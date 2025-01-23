#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "CGSolver.h"
#include "CRSmatfuncs.h"

void generate_esure(int grid_width, int grid_height, int *esure)
{
    int num_elements = grid_width * grid_height;
    for (int i = 0; i < num_elements; i++)
    {
        int row = i / grid_width;
        int col = i % grid_width;

        // Determine neighbors
        int bottom = (row > 0) ? (i - grid_width) : -1;            // Bottom boundary (-1)
        int right = (col < grid_width - 1) ? (i + 1) : -2;         // Right boundary (-2)
        int top = (row < grid_height - 1) ? (i + grid_width) : -3; // Top boundary (-3)
        int left = (col > 0) ? (i - 1) : -4;                       // Left boundary (-4)

        // Store neighbors in the esure array
        esure[i * 4 + 0] = bottom;
        esure[i * 4 + 1] = right;
        esure[i * 4 + 2] = top;
        esure[i * 4 + 3] = left;
    }
}
void generate_cell_stat(int nelem, int *esure, int nrpoin, int *cell_stat)
{
    for (int ele = 0; ele < nelem; ele++)
    {
        for (int i = 0; i < nrpoin; i++)
        {
            int neighbor = esure[nrpoin * ele + i];
            if (neighbor < 0)
            {
                cell_stat[ele] = neighbor;
                break;
            }
        }
    }
}
int main()
{
    // clock
    clock_t start_time, end_time;
    double cpu_time_used;
    // mesh
    int grid_width = 1000;
    int grid_height = 1000;
    int nelem = grid_height * grid_width;
    // int npoin = 25;
    int nrpoin = 4;
    int *esure = (int *)malloc((size_t)nelem * 4 * sizeof(int));
    int *cell_stat = calloc((size_t)nelem, sizeof(int));
    if (esure == NULL || cell_stat == NULL)
    {
        printf("Memory allocation failed.\n");
        return 1;
    }
    generate_esure(grid_width, grid_height, esure);
    generate_cell_stat(nelem, esure, nrpoin, cell_stat);
#ifdef DEBUG
    // Print the esure array
    printf("int esure[] = {\n");
    for (int i = 0; i < nelem; i++)
    {
        printf("    %d, %d, %d, %d", esure[i * 4 + 0], esure[i * 4 + 1], esure[i * 4 + 2], esure[i * 4 + 3]);
        if (i < nelem - 1)
        {
            printf(",\n");
        }
    }
    printf("\n};\n");
    printf("int cell_stat:\n");
    for (int i = 0; i < grid_height; i++)
    {
        for (int j = 0; j < grid_width; j++)
            printf("%d ", cell_stat[grid_height * i + j]);
        printf("\n");
    }
#endif
    // define sparse matrix of coefficient
    int *row_ptr = (int *)calloc((size_t)(nelem + 1), sizeof(int));
    int max_nnz = nelem * 6; // At most 4 non-zero entries per row
    double *val = (double *)malloc((size_t)max_nnz * sizeof(double));
    int *col_ind = (int *)malloc((size_t)max_nnz * sizeof(int));
    int nnz = 0;
    double *RHS = (double *)calloc((size_t)nelem, sizeof(double));
    if (row_ptr == NULL || val == NULL || col_ind == NULL || RHS == NULL)
    {
        printf("Memory allocation failed.\n");
        return 1;
    }
    // calculate the symmetrical Positive Definite Coefficient martix and Right-Hand-Sided vector
    for (int ele = 0; ele < nelem; ele++)
    {
        nnz++;
        int IDele = nnz - 1;
        // Known cells 
        if (cell_stat[ele]<0){
            // equation for this element is u=Tb
            val[IDele] = 1;
            col_ind[IDele] = ele;
            RHS[ele] = (cell_stat[ele] == -1) ? 5 : 1; 
        }else{
        // Unknown cells 
            for (int nei = 0; nei < nrpoin; nei++)
            {
                int neighbor = esure[nrpoin * ele + nei];
                if (neighbor >= 0)
                {   
                    // Internal flux  
                    // contribution for diagonal element of coefficient matrix 
                    val[IDele] += 1;
                    col_ind[IDele] = ele;
                    // contribution for off-diagonal element of coefficient matrix 
                    if (cell_stat[neighbor] <0 ){
                        // between known and unknown cells
                        RHS [ele] += (cell_stat[neighbor] == -1) ? 5 : 1;
                    }else{
                        // between 2known cells 
                        val[nnz] = -1;
                        col_ind[nnz] = neighbor;
                        nnz++; // cause creat new element in the coefficient matrix
                    }

                }
                else if (neighbor == -1)
                { // Boundary condition
                // coeff[nelem * ele + ele] += 2;               // Subtract halved contribution (distance halved)
                    val[IDele] += 2;
                    RHS[ele] += 5 * 2; // Add scaled boundary condition to RHS
                }
                else
                {
                    { // Boundary condition
                    // coeff[nelem * ele + ele] += 2;               // Subtract halved contribution (distance halved)
                        val[IDele] += 2;
                        RHS[ele] += 1 * 2; // Add scaled boundary condition to RHS
                    }
                }
            }
        }
        row_ptr[ele + 1] = nnz;
    }
    row_ptr[nelem] = nnz;

#ifdef DEBUG
    // Print the CRS representation
    printf("Values (val): ");
    for (int i = 0; i < nnz; i++)
    {
        printf("%.2lf ", val[i]);
    }
    printf("\n");
    printf("Column Indices (col_ind): ");
    for (int i = 0; i < nnz; i++)
    {
        printf("%d ", col_ind[i]);
    }
    printf("\n");
    printf("Row Pointers (row_ptr): ");
    for (int i = 0; i <= nelem; i++)
    {
        printf("%d ", row_ptr[i]);
    }
    printf("\n");
    printf("RHS: ");
    for (int i = 0; i < nelem; i++)
    {
        printf("%lf ", RHS[i]);
    }
    printf("\n");
#endif
    // Check the Symetrical Positive Difinte conditions
    if (nelem > 10000)
    {
        printf("! skipping the SPD condition because calculation load of this case (nelem=%d)>10000 is high\n", nelem);
    }
    else if (isPositiveDefinite(nelem, val, col_ind, row_ptr))
    {
        printf("The matrix is positive definite.\n");
    }
    else
    {
        fprintf(stderr, "The Coefficient matrix is not positive definite.\n");
        exit(EXIT_FAILURE);
    }

    // double b [] = {1, 1, 1, 1, 1, 2, 2, 1, 1, 6, 6, 1, 5, 5, 5, 5};
    // double vall[] = {1, 1, 1, 1, 1, 4, -1, -1, 4, -1, -1, 1, 1, 4, -1, -1, 4, -1, -1, 1, 1, 1, 1, 1};
    // int col_indx[] = {0, 1, 2, 3, 4, 5, 6, 9, 6, 5, 10, 7, 8, 9, 5, 10, 10, 6, 9, 11, 12, 13, 14, 15};
    // int row_p[] = {0, 1, 2, 3, 4, 5, 8, 11, 12, 13, 16, 19, 20, 21, 22, 23, 24};
    // int row = 16;
    // int nonzeros = 24;
    // CRSMatrix A;
    // A.n = row;
    // A.nnz = nonzeros;
    // A.row_ptr = row_p;
    // A.col_index = col_indx;
    // A.values = vall;

    CRSMatrix A;
    A.n = nelem;
    A.nnz = nnz;
    A.row_ptr = row_ptr;
    A.col_index = col_ind;
    A.values = val;

    // unknown vector
    double *u = (double *)calloc((size_t)A.n, sizeof(double));

    // Solve using CG
    SolverConfig config = {100000, 1e-14, true};
    solver_set_config(config);
    // precond_conjugate_gradient(&A, RHS, u);
    start_time = clock();
    conjugate_gradient(&A, RHS, u);
    end_time = clock();
    cpu_time_used = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("! CG Solver execution time : %.2f seconds\n", cpu_time_used);

    // // Output solution
    // printf("Solution u:\n");
    // for (int i = 0; i < grid_height; i++)
    // {
    //     for (int j = 0; j < grid_width; j++)
    //         printf("%f ", u[grid_height * i + j]);
    //     printf("\n");
    // }

    free(val);
    free(col_ind);
    free(row_ptr);
    free(RHS);
    free(u);
    free(esure);
    free(cell_stat);
    return 0;
}
