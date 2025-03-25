#include <stdio.h>
#include <stdlib.h>
#include <omp.h>


void usage(char prog_name[]);
double serialPi(double step);
double falseSharingPi(double step);
double raceConditionPi(double step);
double noRaceConditionPi(double step);
static long num_steps = 1000000;
int NUM_THREADS;

int main(int argc, char const *argv[])
{
    double start_time, run_time=0, pi, step, tSerial,tParallel,speedup;
	int iter;


	if (argc != 2) {
		usage(argv[0]);
		exit (-1);
	}
	
    NUM_THREADS = omp_get_max_threads();
    omp_set_num_threads(NUM_THREADS);

    run_time=0;
    iter=atoi(argv[1]);
	step = 1.0/(double)num_steps;
	for(int i=0; i<iter; i++){
        // Record the time for staring the computation of pi 
		start_time = omp_get_wtime();
		pi=serialPi(step);
        // Record the time for ending of the computation of pi and
		//find the elapsed time for the computation of pi
		//Note we are running the same process for multiple time and taking the average here 
		run_time += omp_get_wtime() - start_time;
	}
    tSerial = run_time/iter;
	printf("\n Sequential: pi with %ld steps is %f in %f seconds\n",num_steps,pi,tSerial);
    
    run_time=0;
    iter=atoi(argv[1]);
	step = 1.0/(double)num_steps;
    for(int i=0; i<iter; i++){
        // Record the time for staring the computation of pi 
		start_time = omp_get_wtime();
		pi=falseSharingPi(step);
        // Record the time for ending of the computation of pi and
		//find the elapsed time for the computation of pi
		//Note we are running the same process for multiple time and taking the average here 
		run_time += omp_get_wtime() - start_time;
	}
    tParallel = run_time/iter;
    speedup = tSerial / tParallel;
	printf("\n Parallel with false sharing: pi with %ld steps is %f in %f seconds using %d threads\n",num_steps,pi,tParallel,NUM_THREADS);
    printf("Speedup: %f\n",speedup);
    
    run_time=0;
    iter=atoi(argv[1]);
	step = 1.0/(double)num_steps;
    for(int i=0; i<iter; i++){
        // Record the time for staring the computation of pi 
		start_time = omp_get_wtime();
		pi=raceConditionPi(step);
        // Record the time for ending of the computation of pi and
		//find the elapsed time for the computation of pi
		//Note we are running the same process for multiple time and taking the average here 
		run_time += omp_get_wtime() - start_time;
	}
    tParallel =run_time/iter;
    speedup = tSerial / tParallel;
	printf("\n Parallel with race Condtion: pi with %ld steps is %f in %f seconds using %d threads\n",num_steps,pi,tParallel,NUM_THREADS);
    printf("Speedup: %f\n",speedup);

    run_time=0;
    iter=atoi(argv[1]);
	step = 1.0/(double)num_steps;
    for(int i=0; i<iter; i++){
        // Record the time for staring the computation of pi 
		start_time = omp_get_wtime();
		pi=noRaceConditionPi(step);
        // Record the time for ending of the computation of pi and
		//find the elapsed time for the computation of pi
		//Note we are running the same process for multiple time and taking the average here 
		run_time += omp_get_wtime() - start_time;
	}
    tParallel =run_time/iter;
    speedup = tSerial / tParallel;
	printf("\n Parallel with no race Condtion: pi with %ld steps is %f in %f seconds using %d threads\n",num_steps,pi,tParallel,NUM_THREADS);
    printf("Speedup: %f\n",speedup);
    return EXIT_SUCCESS; 
}

double serialPi(double step){
    double x, pi, sum = 0.0;

    for(int  i = 0; i < num_steps; i++){
        x = (i + 0.5 ) * step;
        sum += 4.0 / (1 + (x*x));
    }

    pi = step * sum;
    return pi;
}

double falseSharingPi(double step){

    int nthreads;
    double pi;
    double sum[NUM_THREADS];
    
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel
    {
        int i,id, tthreads; 
        double x;

        tthreads = omp_get_num_threads();
        id = omp_get_thread_num();
        if(id == 0){
            nthreads = tthreads;
        }
        
        for(i = id,sum[id] = 0.0; i < num_steps; i += tthreads){
            x = (i + 0.5 ) * step;
            sum[id] += 4.0 / (1 + (x*x));
        }
    }
    pi = 0.0;
    for (int i = 0; i < nthreads; i++)
    {
        /* code */
        pi += sum[i] * step;
    }
    
    return pi;
}

double raceConditionPi(double step){

    
    double pi = 0.0;
    double sum = 0.0;    

    #pragma omp parallel
    {  
        int id, tthreads; 
        double x;
        tthreads = omp_get_num_threads();
        id = omp_get_thread_num();
        for(int i = id; i < num_steps; i += tthreads){
            x = (i + 0.5 ) * step;
            sum += 4.0 / (1.0 + (x*x));
        }

        pi += step * sum;   
    }

     
    
    return pi;
}

double noRaceConditionPi(double step){
    
    double pi = 0.0;

    #pragma omp parallel
    {  
        int id, tthreads; 
        double x;
        tthreads = omp_get_num_threads();
        id = omp_get_thread_num();
        double sum = 0.0;
        for(int i = id; i < num_steps; i += tthreads){
            x = (i + 0.5 ) * step;
            sum += 4.0 / (1.0 + (x*x));
        }

        #pragma omp atomic
        pi += step * sum;    

    } 
    
    return pi;
}

/*--------------------------------------------------------------------
 * Function:    usage
 * Purpose:     Print command line for function
 * In arg:      prog_name
 */
void usage(char prog_name[]) {
    fprintf(stderr, "usage:  %s <number of times to run>\n", prog_name);
 } /* usage */
 