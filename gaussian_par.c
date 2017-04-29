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
#define NUM_THREADS 4
pthread_t threads[NUM_THREADS];
int threadID[NUM_THREADS];

pthread_mutex_t counterMutex;
pthread_cond_t counterCond;
int counter;
int read;

double divider;
matrix tempMatrix;


void ReadLock(int iter)
{
	pthread_mutex_lock(&counterMutex);
	while(iter != counter)
	{
		pthread_cond_wait(&counterCond, &counterMutex);
		printf("Awakened\n");
		
	}
	pthread_mutex_unlock(&counterMutex);
}
void ReadUnlock()
{
	pthread_mutex_lock(&counterMutex);
	read++;
	printf("Read = %d\n", read);
	if(read == NUM_THREADS)
	{
		printf("Signaling\n");
		
		pthread_cond_broadcast(&counterCond);
		printf("Signaling done\n");
	}
		
	pthread_mutex_unlock(&counterMutex);
}
void WriteLock()
{
	pthread_mutex_lock(&counterMutex);
	while(read != NUM_THREADS)
	{
		printf("Waiting, read = %d\n", read);
		pthread_cond_wait(&counterCond, &counterMutex);
		printf("Signaled, read = %d\n", read);
		
	}
	pthread_mutex_unlock(&counterMutex);
}

void WriteUnlock()
{
	pthread_mutex_lock(&counterMutex);
	counter++;
	read = 0;
	pthread_cond_broadcast(&counterCond);
	pthread_mutex_unlock(&counterMutex);
}

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
    
	pthread_mutex_init(&counterMutex, NULL);
	pthread_cond_init(&counterCond, NULL);
	counter = 0;
	read = NUM_THREADS;
	
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
	int myID = *(int*)arg;
	
	printf("Thread %d starting...\n", myID);
	
	sleep(2);
	
	for (i = 0; i < 1; i++)
	{	
		if(myID == i % NUM_THREADS) // If the current row to be divided belongs to this thread 
		{
			printf("Thread %d wants to lock, counter = %d, read = %d\n", myID, counter, read);
			WriteLock();
			
			if(i > 0)
			{
				// Copy from temp matrix to A
				for(k = i; k < N; k++)
					A[i][k] = tempMatrix[i][k];
			}
			
			
			divider = 1.0 / A[i][i];	// Calc divider
			A[i][i] = 1.0; 
			printf("Thread %d is writing\n", myID);
			WriteUnlock();		
		}


		ReadLock(i + 1);				// Lock for reading the divider, this is blocked until the writer unlocks increasing the counter to i + 1
		
		printf("Thread %d is reading\n", myID);
		
	
	

		for (k = myID; k < N; k += NUM_THREADS) // The rows the thread should work on
		{
			if(k >= i) // If the row is complete skip it.
			{
				for (j = counter; j < N; j++)
					tempMatrix[k][j] = A[i][j] * divider; // Division step 
				y[k] = b[k]*divider;
			}	
			if(k >= counter) // If the row is complete skip it.
			{
				for (j = counter; j < N; j++)
					A[k][j] = A[i][j] * divider;
				//tempMatrix[k][j];
				//A[k][j] - tempMatrix[k][j] * A[k][i]; // Elimination step 	
				A[k][i] = 0.0;
				b[k] = b[k] - A[k][i]*y[k];
			}					
		}

		ReadUnlock();
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
