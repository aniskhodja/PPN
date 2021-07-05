#include <stdio.h>
#include <mpi.h>
#include<omp.h>
#include <stdlib.h>
double *declare_vectors(unsigned long long length_vector)
{
    double *vect1=(double *) malloc(sizeof(double)*length_vector);
    for (int i = 0; i < length_vector; ++i) {
        //vect1[i]=rand();
        vect1[i]=i;
    }
    return vect1;
}
void print_vectors(double *vect1,double *vect2,unsigned long long length_vector)
{
    printf("vector 1 value:[");
    for (int i = 0; i < length_vector; ++i) {
        printf("%f,",vect1[i]);
    }
    printf("]\n vector 2 value:[");
    for (int i = 0; i < length_vector; ++i) {
        printf("%f,",vect2[i]);
    }
    printf("]\n");
}
int produit_scalaire(double *vect1,double *vect2,unsigned long long length_vector,int rank,int world_size)
{
    //length of process vector
    unsigned long long taille_vect_p;
    //vector of each process
    double *vect1_p=(double*) malloc(sizeof(double )*length_vector);
    double *vect2_p=(double*) malloc(sizeof(double )*length_vector);
    //send data to other process
    MPI_Scatter(vect1,length_vector/world_size,MPI_DOUBLE,vect1_p,length_vector/world_size,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Scatter(vect2,length_vector/world_size,MPI_DOUBLE,vect2_p,length_vector/world_size,MPI_DOUBLE,0,MPI_COMM_WORLD);

    taille_vect_p=length_vector/world_size;
    if(length_vector%world_size!=0)
    {

        if(rank==0)
        {
            MPI_Send(&vect1[length_vector-length_vector%world_size],length_vector%world_size,MPI_DOUBLE,world_size-1,0,MPI_COMM_WORLD);
            MPI_Send(&vect2[length_vector-length_vector%world_size],length_vector%world_size,MPI_DOUBLE,world_size-1,0,MPI_COMM_WORLD);
        }
        else if( rank==world_size-1)
        {
            MPI_Recv(&vect1_p[length_vector/world_size],length_vector%world_size,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            MPI_Recv(&vect2_p[length_vector/world_size],length_vector%world_size,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            taille_vect_p=taille_vect_p+length_vector%world_size;
        }


    }
    //local scaler product
    double sum=0,sum_global;
    for (int i = 0; i < taille_vect_p; ++i) {
        sum=sum+vect1_p[i]*vect2_p[i];
    }
    //sum local
    printf("sum of process %d=%f\n",rank,sum);
    //global scaler product
    MPI_Reduce(&sum,&sum_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    return sum_global;
}
double *declare_matrix(unsigned long long line_length,unsigned long long column_length)
{
    double *vect1=(double *) malloc(sizeof(double)*line_length*column_length);
    for (int i = 0; i < line_length*column_length; ++i) {
        //vect1[i]=rand();
        vect1[i]=i;
    }
    return vect1;
}

double *declare_matrix_column(unsigned long long line_length,unsigned long long column_length)
{
    double *vect1=(double *) malloc(sizeof(double)*line_length*column_length);
    int k=0;
    for (int i = 0; i < column_length; ++i) {
        //vect1[i]=rand();
        for (int j = 0; j < line_length; ++j)
        {
            vect1[i+j*column_length]=k;
            k=k+1;
        }
    }
    return vect1;
}
void print_matrix(double *matrix,double *matrix_1,unsigned long long line_matrix,unsigned long long line_matrix_1,
                  unsigned long long column_matrix,unsigned long long column_matrix_1)
{
    printf("matrix 1:\n[");
    for (int i = 0; i < line_matrix; ++i)
    {
        for (int j = 0; j < column_matrix; ++j)
        {
            printf("%f,",matrix[i*column_matrix+j]);
        }
        printf("]\n");
    }
    printf("matrix 2:\n[");
    for (int i = 0; i < column_matrix_1; ++i)
    {
        for (int j = 0; j < line_matrix_1; ++j)
        {
            printf("%f,",matrix_1[j*column_matrix+i]);
        }
        printf("]\n");
    }
}

void print_matrix_vector(double *vect1,double *vect2,unsigned long long line_matrix,unsigned long long column_matrix,unsigned long long length_vector_2)
{
    printf("matrix value:\n[");
    for (int i = 0; i < line_matrix; ++i) {
        for (int j = 0; j < column_matrix; ++j)
        {
            printf("%f,",vect1[i*column_matrix+j]);
        }
        printf("]\n[");

    }
    printf("]\n vector value:[");
    for (int i = 0; i < length_vector_2; ++i) {
        printf("%f,",vect2[i]);
    }
    printf("]\n");
}
double *produit_matrice_scalaire(unsigned long long line_matrix,unsigned long long column_matrix,unsigned long long length_vector,int world_size,int rank)
{
    double *matrix,*vector=(double*) malloc(sizeof(double)*length_vector),*sub_matrix=(double*) malloc(sizeof(double)*line_matrix*column_matrix);
    if(length_vector!=column_matrix && rank==0)
    {
        printf("produit matrice vecteur impossible %llu!=%llu \n",column_matrix,length_vector);
    }
    else
    {
        MPI_Barrier(MPI_COMM_WORLD);

        //decalaration des variable de taille pour les process
        unsigned long long nbr_ligne_p=line_matrix/world_size;
        unsigned long long line_matrix_process=nbr_ligne_p*column_matrix;

        //decalaration du vecteur resultat et initialisation de la matrice et du vecteur
        double *result=(double*) malloc(sizeof(double)*line_matrix);
        if(rank==0)
        {
            matrix= declare_matrix(line_matrix,column_matrix);
            vector= declare_vectors(length_vector);
            print_matrix_vector(matrix,vector,line_matrix,column_matrix,length_vector);
        }
        //envoi des données
        MPI_Scatter(matrix,nbr_ligne_p*column_matrix,MPI_DOUBLE,sub_matrix,nbr_ligne_p*column_matrix,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(vector,length_vector,MPI_DOUBLE,0,MPI_COMM_WORLD);
        //si le nombre de process ne divise pas le nombre de ligne alors en envoi le restant au dernier process
        if(line_matrix%world_size!=0)
        {

            if(rank==0)
            {
                MPI_Send(&matrix[line_matrix*column_matrix-(line_matrix%world_size)*column_matrix],
                         (line_matrix%world_size)*column_matrix,MPI_DOUBLE,world_size-1,0,MPI_COMM_WORLD);

            }
            else if( rank==world_size-1)
            {
                MPI_Recv(&sub_matrix[nbr_ligne_p*column_matrix],(line_matrix%world_size)*column_matrix,
                         MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            if(rank==world_size-1)
            {
                line_matrix_process=line_matrix_process+(line_matrix%world_size)*column_matrix;
            }

        }
        //calcul du produit local
        double *result_p=(double*) malloc(sizeof(double)*nbr_ligne_p),sum;
        for (int i = 0; i < nbr_ligne_p; ++i) {
            sum=0;
            for (int j = 0; j < column_matrix; ++j) {
                sum=sum+sub_matrix[i*column_matrix+j]*vector[j];
            }

            result_p[i]=sum;
        }
        //envoi des resultat au maitre
        MPI_Gather(result_p,nbr_ligne_p,MPI_DOUBLE,result,nbr_ligne_p,
                   MPI_DOUBLE,0,MPI_COMM_WORLD);
        if(rank==0)
        {
            MPI_Recv(&result[line_matrix-line_matrix%world_size],line_matrix%world_size,
                     MPI_DOUBLE,world_size-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        else if(rank==world_size-1)
        {
            for (int i = nbr_ligne_p; i < nbr_ligne_p+line_matrix%world_size; ++i) {
                sum=0;
                for (int j = 0; j < column_matrix; ++j) {
                    sum=sum+sub_matrix[i*column_matrix+j]*vector[j];
                }

                result_p[i]=sum;
            }

            MPI_Send(&result_p[nbr_ligne_p],
                     (line_matrix%world_size),MPI_DOUBLE,0,0,MPI_COMM_WORLD);

        }
        //affichage resultat
        if(rank==0)
        {
            return result;
        }

    }


    return NULL;

}
int main(int argc, char *argv[]) {
    int rank,world_size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    //main vector
    if(argc==2)
    {
        double *vect1,*vect2;
        //length of main vector
        unsigned long long length_vector=strtol(argv[1], NULL, 10);

        //initialisation of main vector

        if(rank==0)
        {
            vect1= declare_vectors(length_vector);
            vect2= declare_vectors(length_vector);
            print_vectors(vect1,vect2,length_vector);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        //produit scalaire
        int sum_global=produit_scalaire(vect1,vect2,length_vector,rank,world_size);
        //print the result
        MPI_Barrier(MPI_COMM_WORLD);

        if(rank==0)
        {
            printf("sum global %d \n",sum_global);

        }

    }
    else if(argc==4)
    {
        //decalaration des variable de taille pour le main
        unsigned long long line_matrix=strtol(argv[1], NULL, 10);
        unsigned long long column_matrix=strtol(argv[2], NULL, 10);
        unsigned long long length_vector=strtol(argv[3], NULL, 10);
        double *result= produit_matrice_scalaire(line_matrix,column_matrix,length_vector,world_size,rank);
        if(rank==0)
        {
            for (int i = 0; i < line_matrix; ++i) {
                printf("%f ",result[i]);
            }

        }
    }
    else if(argc==5)
    {
        unsigned long long line_matrix=strtol(argv[1], NULL, 10);
        unsigned long long column_matrix=strtol(argv[2], NULL, 10);
        unsigned long long line_matrix_1=strtol(argv[3], NULL, 10);
        unsigned long long column_matrix_1=strtol(argv[4], NULL, 10);
        double *matrix= declare_matrix(line_matrix,column_matrix);
        double *matrix1= declare_matrix_column(line_matrix,column_matrix);

        if(rank==0)
        {
            print_matrix(matrix,matrix1,line_matrix,column_matrix,line_matrix_1,column_matrix_1);
        }


    }
    MPI_Finalize();
    return 0;
}
