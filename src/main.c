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
    // calculate the centeroid of each element
    double *cen;
    CHECK_ERROR(save_centri3(nelem,elems,ptxyz,&cen));
    if (cen == NULL){
        fprintf(stderr,"Memory allocation (cen) failed.\n");
        return 1;        
    }
    // find a boundary condition
    int *cell_stat = (int *)calloc((size_t)nelem,sizeof(int));
    char bcpath [] = {"/dagon1/achitsaz/poisson/temp/a06161.1.flds.zfem.labels"};
    FILE *fptr = fopen(bcpath,"r");
    if (fptr == NULL){
        fprintf(stderr,"BC file does not open properly\n");
        return 1;
    }
    int buffer = 50, nscan;
    char line[buffer];
    char *str;
    for (int ele =0;ele<nelem;ele++){
        str=edit_endline_character(line,buffer,fptr);
        nscan = sscanf(str, "%d", &cell_stat[ele]);
        if (nscan!=1){
            fprintf(stderr,"! there is error in the line %d of BC file.\n",ele+1);
            return 1;
        }
    }
    fclose(fptr);
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
                {"BC", 1, nelem, cell_stat, SCA_int_VTK},
            };
        size_t countele = sizeof(prtelefield) / sizeof(prtelefield[0]);
        FunctionWithArgs prtpntfield[] = {0};
        size_t countpnt = 0;
        CHECK_ERROR(SaveVTK("./", "checkmesh", 0, M1, tri3funcVTK, prtelefield, countele, prtpntfield, countpnt));
    #endif
    // define sparse matrix of coefficient
    int *row_ptr = (int *)calloc((size_t)(nelem + 1), sizeof(int));
    int max_nnz = nelem * 6; // At most 6 non-zero entries per row
    double *val = (double *)malloc((size_t)max_nnz * sizeof(double));
    int *col_ind = (int *)malloc((size_t)max_nnz * sizeof(int));
    int nnz = 0;
    double *RHS = (double *)calloc((size_t)nelem, sizeof(double));
    if (row_ptr == NULL || val == NULL || col_ind == NULL || RHS == NULL)
    {
        printf("Memory allocation failed.\n");
        return 1;
    }

    // // calculate the symmetrical Positive Definite Coefficient martix and Right-Hand-Sided vector
    // for (int ele = 0; ele < nelem; ele++)
    // {
    //     nnz++;
    //     int IDele = nnz - 1;
    //     // Known cells 
    //     if (cell_stat[ele]<0){
    //         // equation for this element is u=Tb
    //         val[IDele] = 1;
    //         col_ind[IDele] = ele;
    //         RHS[ele] = (cell_stat[ele] == -1) ? 5 : 1; 
    //     }else{
    //     // Unknown cells 
    //         for (int nei = 0; nei < nrpoin; nei++)
    //         {
    //             int neighbor = esure[nrpoin * ele + nei];
    //             if (neighbor >= 0)
    //             {   
    //                 // Internal flux  
    //                 // contribution for diagonal element of coefficient matrix 
    //                 val[IDele] += 1;
    //                 col_ind[IDele] = ele;
    //                 // contribution for off-diagonal element of coefficient matrix 
    //                 if (cell_stat[neighbor] <0 ){
    //                     // between known and unknown cells
    //                     RHS [ele] += (cell_stat[neighbor] == -1) ? 5 : 1;
    //                 }else{
    //                     // between 2known cells 
    //                     val[nnz] = -1;
    //                     col_ind[nnz] = neighbor;
    //                     nnz++; // cause creat new element in the coefficient matrix
    //                 }

    //             }
    //             else if (neighbor == -1)
    //             { // Boundary condition
    //             // coeff[nelem * ele + ele] += 2;               // Subtract halved contribution (distance halved)
    //                 val[IDele] += 2;
    //                 RHS[ele] += 5 * 2; // Add scaled boundary condition to RHS
    //             }
    //             else
    //             {
    //                 { // Boundary condition
    //                 // coeff[nelem * ele + ele] += 2;               // Subtract halved contribution (distance halved)
    //                     val[IDele] += 2;
    //                     RHS[ele] += 1 * 2; // Add scaled boundary condition to RHS
    //                 }
    //             }
    //         }
    //     }
    //     row_ptr[ele + 1] = nnz;
    // }
    // row_ptr[nelem] = nnz;





    // free dynamics arraies
    free (elems);
    free (ptxyz);
    free (open);
    free (esure);
    free (esurp);
    free (esurp_ptr);
    free (row_ptr);
    free (val);
    free (col_ind);
    free (RHS);
    free (cen);
    free (M1);
    free (cell_stat);
    return 0; // success signal
}