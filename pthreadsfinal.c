/* Gaussian elimination*/

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#define MAX_SIZE 4096
#define NUMOFCORES  1


typedef double matrix[MAX_SIZE][MAX_SIZE];


int	SIZE;		            // SIZE is the size of the matrix
int	maxnum;		        /* maxnum is the maximum number of elements */
char	*Init;	        // matrix init type	
int	PRINT_SWITCH;		        // PRINT_SWITCH switch		
matrix	A;		        // matrix A		
double	b[MAX_SIZE];	// vector b  
double	y[MAX_SIZE];	// vector y 
pthread_barrier_t barrier;
pthread_mutex_t mutex_variable = PTHREAD_MUTEX_INITIALIZER;
pthread_once_t once_control = PTHREAD_ONCE_INIT;


/* these are forward declarations */
void cal(void*);
void Init_Matrix(void);
void Print_Matrix(void);
void Default_Init(void);
int Read_Options(int, char **);
void Syn_Init(void);

int main(int argc, char **argv)
{
	long i;
	pthread_t cal_thrs[NUMOFCORES];

	Default_Init();
	Read_Options(argc, argv);
	Init_Matrix();
	Syn_Init();

	pthread_barrier_init(&barrier, NULL, NUMOFCORES);

	// creation of threads

	for (i = 0; i < NUMOFCORES; i++)
		pthread_create(&cal_thrs[i], NULL, (void *)&cal, (void *)i);

	// termination of all threads

	for (i = 0; i < NUMOFCORES; i++)
		pthread_join(cal_thrs[i], NULL);

	pthread_barrier_destroy(&barrier);

	if (PRINT_SWITCH == 1)
		Print_Matrix();

}
void  cal(void *thread_id)
{
	long tid = (long)thread_id;
	int i, j, k;
	int c = 0;
	

	pthread_mutex_init(&mutex_variable, NULL);
	pthread_mutex_lock(&mutex_variable);
	if (c)
		pthread_once(&once_control, Syn_Init);
	pthread_mutex_unlock(&mutex_variable);


	/*  Algorithm 8.4 from Grama */
	/* Outer loop */

	for (k = 0; k < SIZE; k++)
	{
		if (tid == (k % NUMOFCORES))                  /* assign each row to a thread in a cyclic fashion */
		{
			for (j = k + 1; j < SIZE; j++)
				A[k][j] = A[k][j] / A[k][k];             /* Division step */
			y[k] = b[k] / A[k][k];
			A[k][k] = 1.0;
		}
		pthread_barrier_wait(&barrier);

		for (i = k + 1; i < SIZE; i++)
		{
			if (tid == (i % NUMOFCORES))               /* assign each row to a thread in a cyclic fashion */
			{
				for (j = k + 1; j < SIZE; j++)
					A[i][j] = A[i][j] - A[i][k] * A[k][j]; /* Elimination step */
				b[i] = b[i] - A[i][k] * y[k];
				A[i][k] = 0.0;
			}
		}
		pthread_barrier_wait(&barrier);
	}

}
void
Syn_Init()
{

	pthread_mutex_init(&mutex_variable, NULL);

}

void
Init_Matrix()
{
	int i, j;

	printf("\nsize      = %dx%d ", SIZE, SIZE);
	printf("\nmaxnum    = %d \n", maxnum);
	printf("Init	  = %s \n", Init);
	printf("Initializing matrix...");

	if (strcmp(Init, "rand") == 0)
	{
		for (i = 0; i < SIZE; i++)
		{
			for (j = 0; j < SIZE; j++)
			{
				if (i == j) /* diagonal dominance */
					A[i][j] = (double)(rand() % maxnum) + 5.0;
				else
					A[i][j] = (double)(rand() % maxnum) + 1.0;
			}
		}
	}
	if (strcmp(Init, "fast") == 0)
	{
		for (i = 0; i < SIZE; i++)
		{
			for (j = 0; j < SIZE; j++)
			{
				if (i == j) /* diagonal dominance */
					A[i][j] = 5.0;
				else
					A[i][j] = 2.0;
			}
		}
	}

	// Initiation of vectors b and y
	for (i = 0; i < SIZE; i++)
	{
		b[i] = 2.0;
		y[i] = 1.0;
	}
	printf("done \n\n");
	if (PRINT_SWITCH == 1)
		Print_Matrix();
}

void
Print_Matrix()
{
	int i, j;

	printf("Matrix A:\n");
	for (i = 0; i < SIZE; i++)
	{
		printf("[");
		for (j = 0; j < SIZE; j++)
			printf(" %5.2f,", A[i][j]);
		printf("]\n");
	}

	printf("Vector b:\n[");
	for (j = 0; j < SIZE; j++)
		printf(" %5.2f,", b[j]);
	printf("]\n");

	printf("Vector y:\n[");
	for (j = 0; j < SIZE; j++)
		printf(" %5.2f,", y[j]);
	printf("]\n");
	printf("\n\n");
}

void
Default_Init()
{
	SIZE = 2048;
	Init = "rand";
	maxnum = 15.0;
	PRINT_SWITCH = 0;
}

int
Read_Options(int argc, char **argv)
{
	char    *prog;
	prog = *argv;

	while (++argv, --argc > 0)
		if (**argv == '-')
			switch (*++*argv)
		{
			case 'n':
				--argc;
				SIZE = atoi(*++argv);
				break;
			case 'h':
				printf("\nHELP: try gauss -u \n\n");
				exit(0);
				break;
			case 'u':
				printf("\nUsage: gauss [-n problemsize]\n");
				printf("           [-D] show default values \n");
				printf("           [-h] help \n");
				printf("           [-I init_type] fast/rand \n");
				printf("           [-m maxnum] max random no \n");
				printf("           [-P print_switch] 0/1 \n");
				exit(0);
				break;
			case 'D':
				printf("\nDefault:  n         = %d ", SIZE);
				printf("\n          Init      = rand");
				printf("\n          maxnum    = 15 ");
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
				PRINT_SWITCH = atoi(*++argv);
				break;
			default:
				printf("%s: ignored option: -%s\n", prog, *argv);
				printf("HELP: try %s -u \n\n", prog);
				break;
		}

}

