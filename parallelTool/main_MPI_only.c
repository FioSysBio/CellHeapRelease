#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>
#include <wait.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include "mpi.h"
#include <time.h>
#include "file_manager.h"
#define MSG_FROM_GLOBAL_MASTER 999
#define NTHREADS 4
#define DONE 0
#define ALMOST_DONE 1

#define PROCESSING 0
#define GO 1
#define PAUSE 2
#define END 3

extern int errno ;


int main(int argc, char *argv[]) {

    /*
    * Input files
    */
    FILE *commandsFile;
    FILE *samplesFile;
    int errnum;

    if ((commandsFile = fopen(argv[1], "r"))==NULL){
        errnum = errno;
        fprintf(stderr,"error while opening commands file, error: %s\n", strerror(errnum));
        exit(1);
    }

    if ((samplesFile = fopen(argv[2], "r"))==NULL){
        errnum = errno;
        fprintf(stderr,"error while opening samples file, error: %s\n", strerror(errnum));
        exit(1);
    }

    int i;
    int j;
    int k;
    char buffer[MAX_BUFFER_CHAR];

    //commands input
	int countCommands=0;
    while (fscanf(commandsFile, "%[^\n] ", buffer) != EOF) {
        countCommands++;
    }
    fprintf(stderr,"number of commands = %d\n", countCommands);
    rewind(commandsFile);
    
    char commandsMatrix[countCommands][MAX_BUFFER_CHAR];

    for(j=0; j<countCommands; j++){
        for(k; k< MAX_BUFFER_CHAR; k++){
            buffer[k]='\0';
        }
		fscanf(commandsFile, "%[^\n] ", buffer);
		strcpy(commandsMatrix[j], buffer);
	}


    //samples input
	int countSamples=0;
    while (fscanf(samplesFile, "%[^\n] ", buffer) != EOF) {
        countSamples++;
    }
    fprintf(stderr,"number os samples ids .sra = %d\n", countSamples);
    rewind(samplesFile);

    char samplesMatrix[countSamples][MAX_BUFFER_CHAR];

	for(j=0; j<countSamples; j++){
        for(k; k< MAX_BUFFER_CHAR; k++){
            buffer[k]='\0';
        }
        fscanf(samplesFile, "%[^\n] ", buffer);
		strcpy(samplesMatrix[j], buffer);
	}

    char executeMatrix[countSamples][countCommands][MAX_BUFFER_CHAR];

    for(j=0; j<countSamples; j++){
        for(k=0; k< countCommands; k++){
                makeQueueOutOfCommandsAndSample (commandsMatrix, countCommands, samplesMatrix[j], executeMatrix[j]);
        }
    }

    /***********************************************************************************/

    int numprocs, rank, namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int init;
    int provided, required=MPI_THREAD_MULTIPLE;
    init=MPI_Init(&argc, &argv);
    if(init!=MPI_SUCCESS){
        fprintf(stderr,"Error while initializing MPI_Init_thread: %d, rank: %d\n", init, rank);
        exit(1);
    }
    init=MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    if(init!=MPI_SUCCESS){
        fprintf(stderr,"Error in Comm_size: %d, rank: %d\n", init, rank);
    }
    init=MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(init!=MPI_SUCCESS){
        fprintf(stderr,"Error in Comm_rank: %d, rank: %d\n", init, rank);
    }
    init=MPI_Get_processor_name(processor_name, &namelen);
    if(init!=MPI_SUCCESS){
        fprintf(stderr,"Error in MPI_Get_processor_name: %d, rank: %d\n", init, rank);
    }
    fprintf(stderr,"at MPI, rank: %d, MPI_COMM_WORLD: %d, processors name: %s\n", rank, MPI_COMM_WORLD, processor_name);


    /***********************************************************************************/

    //take initial time
    time_t initial_time, current_time;
    initial_time = time(NULL);
    fprintf(stderr,"current time inicial: %d\n", initial_time);

    /***********************************************************************************/


    /***********************************************************************************/

        //enviromnet dependent code. Used for eventuals retry during fasterq-dump
		int sysReturn;
		double duration;
        int cmd_rm[] = "rm ";
        int path[] = "/scratch/inova-covd19/helena.silva/MeusTestes/experimento01/";
        int content_rm_fastq = "*.fastq";
        int content_rm_gz_1 = "*sra_1.fastq";
        int content_rm_gz_2 = "*sra_2.fastq";
        int buffer_rewind[MAX_BUFFER_CHAR] = "";

		int iterator=0;
        int position = rank + numprocs*iterator;
		time_t current_time2;

		while (position < countSamples)
		{

			//execute until the penultimate one
			for (i=0; i<countCommands-1; i++){
				//printf("dentro da thread %d, processo %d, execute[%d][%d]: %s\n", iam, rank, position, iterator, executeMatrix[position][i]);
				sysReturn = system(executeMatrix[position][i]);
				if(sysReturn!=0){
				    errnum=errno;
				    fprintf(stderr,"error %s (%d) status system = %d, rank %d\n=> %s\n", strerror(errnum), errnum, sysReturn, rank, executeMatrix[position][i]);
			}


			//execute the last one
			sysReturn = system(executeMatrix[position][i]);
			if(sysReturn!=0){
			        fprintf(stderr,"error in rank %d %s\n", rank, executeMatrix[position][i]);
			}
			current_time = time(NULL);
			duration = difftime(current_time, initial_time);
			fprintf(stderr,"finished at rank %d, duration:  %lf, %s\n", rank, duration, executeMatrix[position][0]);                
			iterator++;
			position = rank + numprocs*iterator;
			
		}


    /***********************************************************************************/

  MPI_Finalize();
  return 0;
}



