#include "mesh.h"
#include "common.h"
#include "vector.h"
#include <stdlib.h>

int main (void){
    int npoin, nelem, *elems;
    double *ptxyz;
    char path [] = {"temp/a06161.1.flds.zfem"};
    CHECK_ERROR(read_zfem(path,&npoin,&nelem,&ptxyz,&elems));
    if (ptxyz == NULL || elems == NULL)
    {
        fprintf(stderr,"Memory allocation (elems/ptxyz) failed.\n");
        return 1;
    }
    // created required data structure for mesh 
    int Nredge = 3;
    int *esurp,*esurp_ptr,*esure,*open;
    save_esurp(npoin,nelem,elems,&esurp,&esurp_ptr,Nredge);
    if (ptxyz == NULL || elems == NULL)
    {
        fprintf(stderr,"Memory allocation (esurp/esurp_ptr) failed.\n");
        return 1;
    }
    save_esure(nelem,elems,esurp_ptr,esurp,&esure,&open,Nredge);
    if (open == NULL || esure == NULL)
    {
        fprintf(stderr,"Memory allocation (esure/open) failed.\n");
        return 1;
    }
    int opencount = 0;
    for (int i=0;i<nelem;i++){
        if (open[i]==0) opencount++;
    }
    (opencount==0) ? printf("! this is open mesh.\n") : printf("* this is close mesh.\n");
    // creat mesh struct 
    mesh *M1 = (mesh *)malloc(sizeof(mesh));
    if (M1)
    {
        *M1 = (mesh){0}; // Set all integer and pointer fields to 0 or NULL
    }
    if (M1 == NULL)
    {
        fprintf(stderr, "Memory allocation failed for M1 pointer\n");
        exit(EXIT_FAILURE);
    }
    M1->elems =elems;M1->npoin =npoin;M1->nelem=nelem;M1->esure=esure;M1->esurp=esurp;M1->esurp_ptr=esurp_ptr;M1->nredge=Nredge;M1->open=open;M1->ptxyz=ptxyz;
    #ifdef DEBUG
           FunctionWithArgs prtelefield[] =
            {
                {"open", 1, nelem, open, SCA_int_VTK},
            };
        size_t countele = sizeof(prtelefield) / sizeof(prtelefield[0]);
        FunctionWithArgs prtpntfield[] = {0};
        size_t countpnt = 0;
        CHECK_ERROR(SaveVTK("./", "checkmesh", 0, M1, tri3funcVTK, prtelefield, countele, prtpntfield, countpnt));
    #endif
    // // define sparse matrix of coefficient
    // int *row_ptr = (int *)calloc((size_t)(nelem + 1), sizeof(int));
    // int max_nnz = nelem * 6; // At most 6 non-zero entries per row
    // double *val = (double *)malloc((size_t)max_nnz * sizeof(double));
    // int *col_ind = (int *)malloc((size_t)max_nnz * sizeof(int));
    // int nnz = 0;
    // double *RHS = (double *)calloc((size_t)nelem, sizeof(double));
    // if (row_ptr == NULL || val == NULL || col_ind == NULL || RHS == NULL)
    // {
    //     printf("Memory allocation failed.\n");
    //     return 1;
    // }




    // free dynamics arraies
    free (elems);
    free (ptxyz);
    free (open);
    free (esure);
    free (esurp);
    free (esurp_ptr);
    // free (row_ptr);
    // free (val);
    // free (col_ind);
    // free (RHS);
    return 0; // success signal
}