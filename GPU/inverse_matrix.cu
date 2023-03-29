#ifndef GPU_INVERSE_CU_
#define GPU_INVERSE_CU_

#include <iostream>
#include <cublas.h>
#include <stdio.h>
#include <fstream>
#include <string.h>
#include <omp.h>

using namespace std;

#include "struct.cu"
//#include "kernel.cu"
#include "MatrixVectorProduct.cu"

const int warp=32;

void InverseMatrix(
		CUSMA::CSRmatrix A,
		double* x);

void SInverseMatrix(
		CUSMA::CSRmatrix A,
		double* x);

extern "C" void InverseMatrix(int* row_offsets, int* column_indices, double* values, 
                              double* inv_values, double* inv_values1, int* num_rows, int* num_elem)
{
    CUSMA::CSRmatrix A;
    A.dim = *num_rows;
    A.NumEl = *num_elem;
    A.V = values;
    A.NC = column_indices;
    A.NL = row_offsets;

    int maxDimGrid;
	
    // определение количества и размера блоков
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties (&deviceProp,0);
    int maxDimBlock = deviceProp.maxThreadsDim[0];		// максимальный размер блока
    maxDimGrid = deviceProp.maxGridSize[1];				// максимальный размер сетки

    cout << "Matrix inversion of the method with the completion" <<endl
	 << "Dimension of a matrix " << A.dim << " Nonzero Element = " << A.NumEl << endl
	 << "MaxDimGrid " << maxDimGrid << " MaxDimBlock " << maxDimBlock << endl;


double time1,time2,time3, time4;
cout<<"InverseMatrix:"<<endl;
time1 = omp_get_wtime();
  InverseMatrix(A, inv_values);
time2 = omp_get_wtime();
cout<<"Time = "<<time2-time1<<endl;
cout<<"SInverseMatrix"<<endl;
time3 = omp_get_wtime();
  SInverseMatrix(A, inv_values1);
time4 = omp_get_wtime();
cout<<"Time = "<<time4-time3<<endl;
        
}
void InverseMatrix(
/*			int* row_offsets,
			int* column_indices,
			double* values,
			double* x,
			double* b,
			int* num_rows,
			int* num_elem,
			int* num_iter,
			double* tol,
			int* prectype)
*/
		CUSMA::CSRmatrix A,
		double* x)
{
// cudaThreadExit();
 cudaSetDevice(0);
	double *AV = A.V;//values;
	double *res = x;
	int* ANC = A.NC;//column_indices;
	int* ANL = A.NL;//row_offsets;
	int dim = A.dim;//*num_rows;
	int NumEl = A.NumEl;//*num_elem;
	// указатели на массивы в видеопамяти
	double *d_AV;//, *d_Inverse; 

	int *d_ANL; 
	int *d_ANC;
	int nwarp;
	int maxDimGrid, DimBlock, DimGrid;
	
	// определение количества и размера блоков
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties (&deviceProp,0);
	int maxDimBlock = deviceProp.maxThreadsDim[0];		// максимальный размер блока
	maxDimGrid = deviceProp.maxGridSize[1];				// максимальный размер сетки
/*	cout << "Matrix inversion of the method with the completion" <<endl
		 << "Dimension of a matrix " << dim << " Nonzero Element = " << NumEl << endl
		 << "MaxDimGrid " << maxDimGrid << " MaxDimBlock " << maxDimBlock << endl;
*/
	nwarp = 1 + dim/((maxDimGrid-1)*warp);
	DimBlock = nwarp * warp;
	DimGrid = dim / DimBlock +1;
	// установка количества блоков
	dim3 grid(DimGrid, 1, 1);
	// установка количества потоков в блоке
	dim3 threads(DimBlock, 1, 1);
	// выделение видеопамяти
	cudaMalloc((void **)&d_AV, sizeof(double)*NumEl);
//	cudaMalloc((void **)&d_Inverse, sizeof(double)*dim*dim);
	cudaMalloc((void **)&d_ANL, sizeof(int)*(dim+1));
	cudaMalloc((void **)&d_ANC, sizeof(int)*NumEl);
	
	// копирование из оперативной памяти в видеопамять
	cudaMemcpy(d_AV, AV, sizeof(double)*NumEl, cudaMemcpyHostToDevice);
	cudaMemcpy(d_ANL, ANL, sizeof(int)*(dim+1), cudaMemcpyHostToDevice);
	cudaMemcpy(d_ANC, ANC, sizeof(int)*NumEl, cudaMemcpyHostToDevice);

//	cout << "Memory allocation time - " << t0/1000 << endl;	// время выделения памяти GPU

	double *diagA;			// диагональный предобуславлеватель
	double *ak, *v;
	cudaMalloc((void **)&ak, sizeof(double)*dim);
	cudaMalloc((void **)&v, sizeof(double)*dim);
//	cout << "Time calculation diagonal preconditioner - " << t1 << endl;
	// время вычисления диагонального предобуславлевателя

	int k=0;
// **************************************** Ershow met
        unsigned int i,j;
	double tim1,tim2,result1,buf1,tim3,tim4,result2,buf2;
        ofstream yyy;
        
        double **d_Inverse;
        d_Inverse = new double* [dim];
            
        for(i=0 ; i<dim ; i++)
        {
            cudaMalloc((void **)&d_Inverse[i], sizeof(double)*dim);
            iden_vec <<<grid, threads >>> (d_Inverse[i], dim, i);			// обнуление вектора результата
        }
        double *loc = new double [dim];
/*        yyy.open("MatE");
        for(i=0;i<dim;i++)
        { 	
            cudaMemcpy(loc, d_Inverse[i], sizeof(double)*dim, cudaMemcpyDeviceToHost);

            for(j=0;j<dim;j++)
                yyy<<loc[j]<<"   ";
            yyy<<endl;
        }
        yyy.close();
*/
//	cublasHandle_t handle;
//	cublasCreate(&handle);   
     
/*	cudaMemcpy(res, d_Inverse, sizeof(double)*dim*dim, cudaMemcpyDeviceToHost);

        yyy.open("MatE");
        for(i=0;i<dim*dim;i++)
            yyy<<res[i]<<endl;
        yyy.close();
*/

        double alpha,alpha2,alpha1;
	
        double *loc1 = new double [dim];
        double *loc2 = new double [dim];

//        yyy.open("Vectors");
        for(i=0;i<dim;i++)
        {
            vec_def <<<grid, threads >>> (v , dim);
            SetV <<<grid, threads >>> (v, d_AV, d_ANC, ANL[i+1]-ANL[i], ANL[i], i);
            cublasDcopy(dim, d_Inverse[i], 1, ak, 1);
            for(j=0 ; j<dim ; j++)
            {
//                cudaMemcpy(loc, v, sizeof(double)*dim, cudaMemcpyDeviceToHost);
//                cudaMemcpy(loc1, ak, sizeof(double)*dim, cudaMemcpyDeviceToHost);
//                cudaMemcpy(loc2, d_Inverse[j], sizeof(double)*dim, cudaMemcpyDeviceToHost);
//                for(unsigned k=0;k<dim;k++)
//                   yyy<<loc[k]<<"   "<<loc1[k]<<"   "<<loc2[k]<<endl;
//                yyy<<"***"<<endl;
        tim1 = omp_get_wtime();
	alpha1=cublasDdot(dim, v, 1, d_Inverse[j], 1);
	tim2 = omp_get_wtime();
	buf1=tim2-tim1;
	result1+=buf1;

	tim3 = omp_get_wtime();
	alpha2=1+cublasDdot(dim, v, 1, ak, 1);
	tim4 = omp_get_wtime();
	buf2=tim4-tim3;
	result2+=buf2;

        alpha =alpha1/alpha2;
                cublasDaxpy(dim, -alpha, ak, 1, d_Inverse[j], 1);
            }
//            yyy<<"--------------------------"<<endl;
        
        }
//        yyy.close();
	yyy.open("time_vect_ershov.txt",ios::app);
	yyy<<"razmernost="<<dim<<endl;
	yyy<<"vectprod 1="<<result1<<endl;
	yyy<<"vectprod 2="<<result2<<endl;
	yyy<<"END"<<endl;
	yyy.close();
        maxDimGrid = deviceProp.maxGridSize[1];
        yyy.open("E");
          for(i=0;i<dim;i++)
        {
            MatrVectMul(maxDimGrid, v, d_Inverse[i], d_AV, d_ANC, d_ANL, dim, NumEl);
            cudaMemcpy(loc, v, sizeof(double)*dim, cudaMemcpyDeviceToHost);

            for(j=0;j<dim;j++)
                yyy<<loc[j]<<"   ";
            yyy<<endl;   
        }
        yyy.close();
	// копирование результата из видеопамяти в оперативную память
	for(i=0;i<dim;i++)
            cudaMemcpy(res+dim*i, d_Inverse[i], sizeof(double)*dim, cudaMemcpyDeviceToHost);
// *******************************************
//	cout << "Quentity of iterations " << k <<" Error: "<<nrmr/nrmb<< endl
//		 << "Time execution of the conjugate gradients method - " << t2/1000 << endl
// время решения системы методом сопряженных градиентов
//		 << "Total time - " << t/1000 << endl;
// общее время работы функции CG
//  освобождение памяти

        delete [] loc1;
        delete [] loc2;
        delete [] loc;
        for (i=0;i<dim;i++)
        cudaFree(d_Inverse[i]);

	delete [] d_Inverse;
	cudaFree(d_AV);
        cudaFree(d_ANL);
	cudaFree(d_ANC);
	cudaFree(ak);
	cudaFree(v);
//	cublasDestroy(handle);
	
	cudaDeviceReset();
}
// end InverseMatrix
// *********************************************

void MatrVectMul (const int maxDimGrid, 
	      double* Ax , 
	  const double *p, 
	  const int *d_ANL, 
	  const int dim,
	  const int NumEl)
{
    
    const size_t dimBlock = 128;                                              //количество нитей в блоке
    
    int nnz_per_row = NumEl / dim;                                            //среднее количество ненулевых элементов в строке
    unsigned int thr_per_vec;                                                 //количество нитей для вычисления одной координаты вектора
    
    if (nnz_per_row <=  2) thr_per_vec=2;
    else
    if (nnz_per_row <=  4)  thr_per_vec=4;
    else
    if (nnz_per_row <=  8)  thr_per_vec=8;
    else
    if (nnz_per_row <=  16)  thr_per_vec=16;
    else
    thr_per_vec=32;
    
    const size_t VECTORS_PER_BLOCK = dimBlock / thr_per_vec;         // количество векторов в блоке
    
    const size_t DimGrid = std::min<int>(maxDimGrid, (dim + VECTORS_PER_BLOCK-1)/VECTORS_PER_BLOCK); //количество блоков
    
    cudaBindTexture(0, tex_b, p, sizeof(double)*dim);                                                   //"Привязка текстуры"
    
    dev_MatrVectMul <<<DimGrid, dimBlock>>> (Ax ,  d_ANL, dim, VECTORS_PER_BLOCK, thr_per_vec);  //выхов функции на девайсе
        //dev_MatrVectMul < VECTORS_PER_BLOCK, THREADS_PER_VECTOR > <<<DimGrid, dimBlock>>> (Ax ,  d_ANL, dim);
    cudaUnbindTexture(tex_b);                                                                           //"Отвязка текстуры"
}


void CG(int* d_ANL,
    double* d_b, 
    int dim, 
    int NumEl, 
    int iter, 
    double acc,
    double *diagA,
    double *r,
    double *Ax,
    double *p,
    double *z,
    double *d_res,
    dim3 grid,
    dim3 threads,
    int maxDimGrid,
    int *k,
    double *buf1,double *buf2,double *buf3, double *buf4)
{
//cout<<"check1"<<endl;

    /*cublasHandle_t handle;
    cublasCreate(&handle);*/
    double    buff1,buff2,buff3,buff4;
    double ro1, ro0, alpha, alpha1, beta, nrmr, nrmb;
    int kk=0;buff1=0;buff2=0;buff3=0;buff4=0;
    double time1, time2,time3,time4,time5,time6,time7,time8,result1,result2,result3,result4;

//****************************************
//метод сопряжённых градиентов

/*    double *local = new double [dim];
    for(int j=0;j<dim;j++)
        local[j] = 0;
    local[number]=1;
    

    double *d_b;
    cudaMalloc((void **)&d_b, sizeof(double)*dim);

    cudaMemcpy(d_b, local, sizeof(double)*dim, cudaMemcpyHostToDevice);

    cudaMemcpy(local, d_b, sizeof(double)*dim, cudaMemcpyDeviceToHost);
    for(int j=0;j<dim;j++)
        cout<<local[j]<<"   ";
    cout<<endl;*/

    cublasDcopy(dim,d_b,1,r,1);                              // r = d_b

    double nrmb1, nrmr1;
    nrmb1 = cublasDdot(dim,d_b,1,d_b,1);                      // nrmb1 = (d_b, d_b)
    nrmb=sqrt(nrmb1);
    nrmr=nrmb;
//cout<<"check2"<<endl;

    DiagVectMul<<< grid, threads >>>(z ,diagA, r, dim);             // z = M*r, M - диагональный предобуславлеватель
//cout<<"check3"<<endl;
    cublasDcopy (dim,z,1,p,1);                               // p = z
    ro1 = cublasDdot (dim,z,1,r,1);                           //ro1 = (z,r)
//    cout<<"nrmb = "<<nrmb<<endl;
//cout<<"check4"<<endl;

    if (nrmb)
    while ( ((nrmr/nrmb) > acc) && (kk <= iter ) )
    {
	time1 = omp_get_wtime();
	MatrVectMul(maxDimGrid, Ax , p, d_ANL, dim, NumEl);	 //матрично-векторное произведение  Ax = A * p
	time2 = omp_get_wtime();
	result1=time2-time1;
	buff1+=result1;

	alpha1=0;

	time3 = omp_get_wtime();
	alpha1 = cublasDdot(dim,p,1,Ax,1);                 //alpha1 = (p, Ax)
	time4 = omp_get_wtime();
	result2=time4-time3;
	buff2+=result2;

	alpha=ro1 / alpha1;		
	cublasDaxpy(dim,alpha,p,1,d_res,1) ;             //d_res = d_res + alpha * p
	alpha1=-alpha;
	cublasDaxpy(dim,alpha1,Ax,1,r,1) ;               //r = r - alpha * Ax
	
	time5 = omp_get_wtime();
	nrmr1 = cublasDdot(dim,r,1,r,1) ;                  // nrmr1 = (r,r)
	time6 = omp_get_wtime();
	result3=time6-time5;
	buff3+=result3;

	nrmr=sqrt(nrmr1);
	DiagVectMul<<< grid, threads>>>(z ,diagA, r, dim);      //z = M*r, M - диагональный предобуславлеватель
	ro0 = ro1 ;

	time7 = omp_get_wtime();
	ro1 = cublasDdot (dim , z , 1 , r , 1 );          //ro1 = (z,r)
	time8 = omp_get_wtime();
	result4=time8-time7;
	buff4+=result4;

	beta=ro1 / ro0 ;
	cublasDscal (dim , beta , p , 1 ) ;             //p = beta*p
	alpha1=1.0;
	cublasDaxpy (dim ,alpha1 , z , 1 , p , 1 ) ;    //p = p + z
	kk++;
    }
//конец метода сопряжённых градиентов
//*******************************************
//cout<<"check5"<<endl;
*k=kk;
*buf1=buff1;
*buf3=buff3;
*buf2=buff2;
*buf4=buff4;

// освобождение памяти
//	cublasDestroy(handle);
}


void SInverseMatrix(
		CUSMA::CSRmatrix A,
		double* x)
{
// cudaThreadExit();
 cudaSetDevice(0);
	double *AV = A.V;//values;
	double *res = x;
	int* ANC = A.NC;//column_indices;
	int* ANL = A.NL;//row_offsets;
	int dim = A.dim;//*num_rows;
	int NumEl = A.NumEl;//*num_elem;
	
	long Norma_A=0,Sum=0;
	for(int i=0;i<A.dim;i++)
	{
	    for(int j=A.NL[i]; j<A.NL[i+1]; j++)
	    {
		Sum+=abs(A.V[j]) * abs(A.V[j]);
	    }
	}
	Norma_A = sqrt(Sum);

	// указатели на массивы в видеопамяти
	double *d_AV;//, *d_Inverse; 
	int *d_ANL; 
	int *d_ANC;
	int nwarp;
	int maxDimGrid, DimBlock, DimGrid;
	
	// определение количества и размера блоков
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties (&deviceProp,0);
	int maxDimBlock = deviceProp.maxThreadsDim[0];		// максимальный размер блока
	maxDimGrid = deviceProp.maxGridSize[1];				// максимальный размер сетки

/*	cout << "Matrix inversion of the method with the completion" <<endl
		 << "Dimension of a matrix " << dim << " Nonzero Element = " << NumEl << endl
		 << "MaxDimGrid " << maxDimGrid << " MaxDimBlock " << maxDimBlock << endl;
*/
	nwarp = 1 + dim/((maxDimGrid-1)*warp);
	DimBlock = nwarp * warp;
	DimGrid = dim / DimBlock +1;
	
	// установка количества блоков
	dim3 grid(DimGrid, 1, 1);
	// установка количества потоков в блоке
	dim3 threads(DimBlock, 1, 1);
	// выделение видеопамяти
	cudaMalloc((void **)&d_AV, sizeof(double)*NumEl);
//	cudaMalloc((void **)&d_Inverse, sizeof(double)*dim*dim);
	cudaMalloc((void **)&d_ANL, sizeof(int)*(dim+1));
	cudaMalloc((void **)&d_ANC, sizeof(int)*NumEl);
	
	// копирование из оперативной памяти в видеопамять
	cudaMemcpy(d_AV, AV, sizeof(double)*NumEl, cudaMemcpyHostToDevice);
	cudaMemcpy(d_ANL, ANL, sizeof(int)*(dim+1), cudaMemcpyHostToDevice);
	cudaMemcpy(d_ANC, ANC, sizeof(int)*NumEl, cudaMemcpyHostToDevice);

        cudaBindTexture(0, tex_AV, d_AV, sizeof(double)*NumEl);
        cudaBindTexture(0, tex_ANC, d_ANC, sizeof(int)*NumEl);

//	cout << "Memory allocation time - " << t0/1000 << endl;	// время выделения памяти GPU

//	cout << "Time calculation diagonal preconditioner - " << t1 << endl;	// время вычисления диагонального предобуславлевателя

	//int k=0;
// ****************************************CG metod
        unsigned int i,j;
	int k,buf_k;
        ofstream yyy;
        
        double **d_Inverse;
        d_Inverse = new double* [dim];
            
        for(i=0 ; i<dim ; i++)
        {
            cudaMalloc((void **)&d_Inverse[i], sizeof(double)*dim);
            vec_def <<<grid, threads >>> (d_Inverse[i], dim);			// обнуление вектора результата
        }

/*        yyy.open("MatE");
        for(i=0;i<dim;i++)
        { 	
            cudaMemcpy(loc, d_Inverse[i], sizeof(double)*dim, cudaMemcpyDeviceToHost);

            for(j=0;j<dim;j++)
                yyy<<loc[j]<<"   ";
            yyy<<endl;
        }
        yyy.close();
*/
//	cublasHandle_t handle;
//	cublasCreate(&handle);   
     
/*	cudaMemcpy(res, d_Inverse, sizeof(double)*dim*dim, cudaMemcpyDeviceToHost);

        yyy.open("MatE");
        for(i=0;i<dim*dim;i++)
            yyy<<res[i]<<endl;
        yyy.close();
*/

	double *diagA;			// диагональный предобуславлеватель
        double *r, *Ax, *p, *z, *d_b;
        cudaMalloc((void **)&r, sizeof(double)*dim);
        cudaMalloc((void **)&Ax, sizeof(double)*dim);
        cudaMalloc((void **)&p, sizeof(double)*dim);
        cudaMalloc((void **)&z, sizeof(double)*dim);
        cudaMalloc((void **)&diagA, sizeof(double)*dim);
        cudaMalloc((void **)&d_b, sizeof(double)*dim);

	double buf1,buf2,buf3,buf4,tim1,tim2,tim3,tim4;
    k=0;buf1=0;buf2=0;buf3=0;buf4=0;
        diag<<< grid, threads >>>(diagA, d_ANL, dim);                   //вычисление диагонального предобуславлевателя


        double acc = 1e-08;
//        double *local = new double [dim];
        for(i=0;i<dim;i++)
        {
//	    cudaMemcpy(d_b, local, sizeof(double)*dim, cudaMemcpyHostToDevice);
            iden_vec <<<grid, threads >>> (d_b, dim, i);			//identifix ed vect
//            line_def <<<grid, threads, 0, stream >>> (d_Inverse[i], dim);
	    CG(d_ANL, d_b, dim, NumEl, dim, acc, diagA, r, Ax, p, z, d_Inverse[i], grid, threads, maxDimGrid,&k,&buf1,&buf2,&buf3,&buf4);
            tim1+=buf1;
	    tim2+=buf2;
	    tim3+=buf3;
	    tim4+=buf4;
	    buf_k+=k;
        }

ofstream print;
print.open("time_CG.txt",ios::app);
print<<"razmernost="<<dim<<endl;
print<<"time for matrix vect mult="<<tim1<<endl;
print<<"time for vectprod ="<<tim2<<endl;
print<<"time for vectprod ="<<tim3<<endl;
print<<"time for vectprod ="<<tim4<<endl;
print<<"END"<<endl;
print.close();
//delete [] local;

        maxDimGrid = deviceProp.maxGridSize[1];	

        double *loc = new double [dim];

        yyy.open("SE");

        for(i=0;i<dim;i++)
        {
            MatrVectMul(maxDimGrid, r, d_Inverse[i], d_AV, d_ANC, d_ANL, dim, NumEl);
            cudaMemcpy(loc, r, sizeof(double)*dim, cudaMemcpyDeviceToHost);

            for(j=0;j<dim;j++)
                yyy<<loc[j]<<"   ";
            yyy<<endl;   
        }
                  
        yyy.close();
ofstream write;
write.open("iter.txt",ios::app);
write<<"razmernost="<<dim<<endl;
write<<buf_k<<endl;
write<<"END"<<endl;
write.close();
                                                                                                                                                                                                                                                                
	// копирование результата из видеопамяти в оперативную память
	for(i=0;i<dim;i++)
            cudaMemcpy(res+dim*i, d_Inverse[i], sizeof(double)*dim, cudaMemcpyDeviceToHost);
// *******************************************
//	cout << "Quentity of iterations " << k <<" Error: "<<nrmr/nrmb<< endl
//		 << "Time execution of the conjugate gradients method - " << t2/1000 << endl	// время решения системы методом сопряженных градиентов
//		 << "Total time - " << t/1000 << endl;											// общее время работы функции CG
	//  освобождение памяти
        cudaUnbindTexture(tex_AV);
        cudaUnbindTexture(tex_ANC);
	cudaFree(d_AV);
	for(i=0;i<dim;i++)
            cudaFree(d_Inverse[i]);
        delete [] d_Inverse;
	cudaFree(d_ANL);
	cudaFree(d_ANC);
	cudaFree(r);
	cudaFree(Ax);
	cudaFree(p);
	cudaFree(z);
	cudaFree(diagA);
	cudaFree(d_b);

//	cublasDestroy(handle);
	
	cudaDeviceReset();
/*
	long Norma_invA=0;
	Sum = 0;
	for(int i=0;i<dim;i++)
	{
	    for(int j=0; j<dim; j++)
	    {
		if (i==j) res[j+dim*i]-=1;
		Sum+=abs(res[j+dim*i]) * abs(res[j+dim*i]);
	    }
	}
	Norma_invA = sqrt(Sum);
*/
	//cout<<"Obuslovl. A = "<<Norma_A * Norma_invA<<" norma A= "<<Norma_A<<;
//	cout<<" norma INV ="<<Norma_invA<<endl;
}
// end SInverseMatrix
// *********************************************


#endif //GPU_INVERSE_CU_
