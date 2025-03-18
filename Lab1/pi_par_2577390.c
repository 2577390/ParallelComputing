#include <stdio.h>
#include <omp.h>
#define NUM_THREADS 2

/*
long num_steps = 1000000;
double step = 1.0 / (double) num_steps;
*/
double serialPi(){

    long num_steps = 1000000;
    double step = 1.0 / (double) num_steps;
    double x, pi, sum = 0.0;

    for(int  i = 0; i < num_steps; i++){
        x = (i + 0.5 ) * step;
        sum += 4.0 / (1 + (x*x));
    }

    pi = step * sum;
    printf("Method 0 = %f \n",pi);
    return pi;
}


double method1Pi(){


    long num_steps = 1000000;
    double step = 1.0 / (double) num_steps;

    int nthreads;
    double pi;
    double sum[NUM_THREADS];
    
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel
    {
        int id, tthreads; 
        double x;

        tthreads = omp_get_num_threads();
        id = omp_get_thread_num();
        if(id == 0){
            nthreads = tthreads;
        }
        sum[id] = 0;
        for(int i = id; i < num_steps; i += tthreads){
            x = (i + 0.5 ) * step;
            sum[id] += 4.0 / (1 + (x*x));
        }
    }

    for (int i = 0; i < nthreads; i++)
    {
        /* code */
        pi += sum[i] * step;
    }
    
    
    printf("Method 1 = %f \n",pi);
    return pi;
}


double method2Pi(){method1Pi();


    long num_steps = 1000000;
    double step = 1.0 / (double) num_steps;
    
    double pi = 0.0;
    double sum = 0.0;

    int id, tthreads; 
    double x;

    #pragma omp parallel
    {  
        tthreads = omp_get_num_threads();
        id = omp_get_thread_num();
        for(int i = id; i < num_steps; i += tthreads){
            x = (i + 0.5 ) * step;
            sum += 4.0 / (1.0 + (x*x));
        }
    }

    pi = sum * step;    
    
    printf("Method 2 = %f \n",pi);
    return pi;
}

int main(int argc, char const *argv[])
{
    /* code */

    serialPi();
    method1Pi();
    method2Pi();
    return 0;
}

