/*****************************************************
 *
 * Gaussian elimination
 *
 * sequential version
 *
 *****************************************************/

#include <stdio.h>
#include <pthread.h>
#define MAX_SIZE 4096

typedef double matrix[MAX_SIZE][MAX_SIZE];

int	N;		/* matrix size		*/
int	maxnum;		/* max number of element*/
char	*Init;		/* matrix init type	*/
int	PRINT;		/* print switch		*/
matrix	A;		/* matrix A		*/
double	b[MAX_SIZE];	/* vector b             */
double	y[MAX_SIZE];	/* vector y             */
#define NUM_THREADS 16
pthread_t threads[NUM_THREADS];
int threadID[NUM_THREADS];

long double divider; // Tried with only one divide but due to floatingpoint arithmetics the answers diviated a bit.
double temp[NUM_THREADS][MAX_SIZE];
pthread_barrier_t bar;



/* forward declarations */
void* work(void* arg);
void Init_Matrix(void);
void Print_Matrix(void);
void Init_Default(void);
int Read_Options(int, char **);

int 
main(int argc, char **argv)
{
    int i, timestart, timeend, iter;

	
    Init_Default();		/* Init default values	*/
    Read_Options(argc,argv);	/* Read arguments	*/
    Init_Matrix();		/* Init the matrix	*/
    
	pthread_barrier_init(&bar, NULL, NUM_THREADS);

	for (i = 0; i < NUM_THREADS; i++)
	{
		threadID[i] = i;
		pthread_create(&threads[i], NULL, work, (void*)&threadID[i]);
	}
	for (i = 0; i < NUM_THREADS; i++)
	{
		pthread_join(threads[i], NULL);
	}
	
    if (PRINT == 1)
	Print_Matrix();
}

void*
work(void* arg)
{
    int i, j, k;
	double div;
	int myID = *(int*)arg;

	for (i = 0; i < N; i++)
	{	
		for (k = myID; k < N; k += NUM_THREADS) // The rows the thread should work on
		{
			if(k == i) // Skip all rows less than i(completed rows)
			{
				for (j = i + 1; j < N; j++)
					temp[myID][j] =A[i][j]/ A[i][i]; // Division step 
				y[i] = b[i] / A[i][i];
			}
			else if(k > i)
			{
				for (j = i + 1; j < N; j++)
					A[k][j] = A[k][j] - A[k][i]* (A[i][j]/ A[i][i]);// Division and Elimination step
				b[k] = b[k] - A[k][i]*(b[i] / A[i][i]);
				A[k][i] = 0.0;
			}				
		}
			
		pthread_barrier_wait(&bar); // wait for everyone to finished this iteration
		if(myID == i % NUM_THREADS) // If the current pivot row belongs to this thread 
		{
			// Copy from temp
			for(j = i + 1; j < N; j++)
				A[i][j] = temp[myID][j];
			A[i][i] = 1.0;		
		}
	}


	
}

void
Init_Matrix()
{
    int i, j;
 
    printf("\nsize      = %dx%d ", N, N);
    printf("\nmaxnum    = %d \n", maxnum);
    printf("Init	  = %s \n", Init);
    printf("Initializing matrix...");
 
    if (strcmp(Init,"rand") == 0) {
	for (i = 0; i < N; i++){
	    for (j = 0; j < N; j++) {
		if (i == j) /* diagonal dominance */
		    A[i][j] = (double)(rand() % maxnum) + 5.0;
		else
		    A[i][j] = (double)(rand() % maxnum) + 1.0;
	    }
	}
    }
    if (strcmp(Init,"fast") == 0) {
	for (i = 0; i < N; i++) {
	    for (j = 0; j < N; j++) {
		if (i == j) /* diagonal dominance */
		    A[i][j] = 5.0;
		else
		    A[i][j] = 2.0;
	    }
	}
    }

    /* Initialize vectors b and y */
    for (i = 0; i < N; i++) {
	b[i] = 2.0;
	y[i] = 1.0;
    }

    printf("done \n\n");
    if (PRINT == 1)
	Print_Matrix();
}

void
Print_Matrix()
{
    int i, j;
 
    printf("Matrix A:\n");
    for (i = 0; i < N; i++) {
	printf("[");
	for (j = 0; j < N; j++)
	    printf(" %5.2f,", A[i][j]);
	printf("]\n");
    }
    printf("Vector b:\n[");
    for (j = 0; j < N; j++)
	printf(" %5.2f,", b[j]);
    printf("]\n");
    printf("Vector y:\n[");
    for (j = 0; j < N; j++)
	printf(" %5.2f,", y[j]);
    printf("]\n");
    printf("\n\n");
}

void 
Init_Default()
{
    N = 2048;
    Init = "rand";
    maxnum = 15.0;
    PRINT = 0;
}
 
int
Read_Options(int argc, char **argv)
{
    char    *prog;
 
    prog = *argv;
    while (++argv, --argc > 0)
	if (**argv == '-')
	    switch ( *++*argv ) {
	    case 'n':
		--argc;
		N = atoi(*++argv);
		break;
	    case 'h':
		printf("\nHELP: try sor -u \n\n");
		exit(0);
		break;
	    case 'u':
		printf("\nUsage: sor [-n problemsize]\n");
		printf("           [-D] show default values \n");
		printf("           [-h] help \n");
		printf("           [-I init_type] fast/rand \n");
		printf("           [-m maxnum] max random no \n");
		printf("           [-P print_switch] 0/1 \n");
		exit(0);
		break;
	    case 'D':
		printf("\nDefault:  n         = %d ", N);
		printf("\n          Init      = rand" );
		printf("\n          maxnum    = 5 ");
		printf("\n          P         = 0 \n\n");
		exit(0);
		break;
	    case 'I':
		--argc;
		Init = *++argv;
		break;
	    case 'm':
		--argc;
		maxnum = atoi(*++argv);
		break;
	    case 'P':
		--argc;
		PRINT = atoi(*++argv);
		break;
	    default:
		printf("%s: ignored option: -%s\n", prog, *argv);
		printf("HELP: try %s -u \n\n", prog);
		break;
	    } 
}
