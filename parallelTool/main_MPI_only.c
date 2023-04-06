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
    //FILE *confFile;
    FILE *commandsFile;
    FILE *samplesFile;
    int errnum;

    if ((commandsFile = fopen(argv[1], "r"))==NULL){
        errnum = errno;
        fprintf(stderr,"erro ao abrir arquivo comandos, erro: %s\n", strerror(errnum));
        exit(1);
    }

    if ((samplesFile = fopen(argv[2], "r"))==NULL){
        errnum = errno;
        fprintf(stderr,"erro ao abrir arquivo samples, erro: %s\n", strerror(errnum));
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
    fprintf(stderr,"numero de comandos = %d\n", countCommands);
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
    fprintf(stderr,"numero de amostras .sra = %d\n", countSamples);
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
        fprintf(stderr,"Erro ao inicializar o MPI_Init_thread: %d, rank: %d\n", init, rank);
        exit(1);
    }
    init=MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    if(init!=MPI_SUCCESS){
        fprintf(stderr,"Erro no Comm_size: %d, rank: %d\n", init, rank);
    }
    init=MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(init!=MPI_SUCCESS){
        fprintf(stderr,"Erro no Comm_rank: %d, rank: %d\n", init, rank);
    }
    init=MPI_Get_processor_name(processor_name, &namelen);
    if(init!=MPI_SUCCESS){
        fprintf(stderr,"Erro no MPI_Get_processor_name: %d, rank: %d\n", init, rank);
    }
    fprintf(stderr,"no MPI, rank: %d, MPI_COMM_WORLD: %d, processors name: %s\n", rank, MPI_COMM_WORLD, processor_name);


    /* Check that the MPI implementation supports MPI_THREAD_MULTIPLE 
    if (provided < MPI_THREAD_MULTIPLE) {
        fprintf(stderr,"MPI does not support MPI_THREAD_MULTIPLE\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
        return 0;
    }*/

    /***********************************************************************************/



    //pega o tempo atual
    time_t initial_time, current_time;
    initial_time = time(NULL);
    fprintf(stderr,"current time inicial: %d\n", initial_time);

    /***********************************************************************************/


    /***********************************************************************************/

        char recebendo[50];
        char devolta[] = "voltei\n";
        char id[40];
        int statusRun=0;
        int iteration = 0;

    /***********************************************************************************/

		int sysReturn;
		double duration;
        int cmd_rm[] = "rm ";
        int path[] = "/scratch/inova-covd19/helena.silva/MeusTestes/experimento01/";
        int content_rm_fastq = "*.fastq";
        int content_rm_gz_1 = "*sra_1.fastq";
        int content_rm_gz_2 = "*sra_2.fastq";
        int buffer_rewind[MAX_BUFFER_CHAR] = "";

		int iterator=0;
		//int offset = (rank==0) ? 2 : 1;
		int offset = 0;
		//int position = (iam-offset)*numprocs + rank + (NTHREADS*numprocs-2*offset)*iterator;
		//int position = (iam-offset)*numprocs + rank + (NTHREADS-offset)*numprocs*iterator;
        int position = rank + numprocs*iterator;
		time_t current_time2;

		while (position < countSamples)
		{

			/*
			*	Nem sempre ele consegue terminar o processo chamado pelo system quando são alocados mais de dois nodes no SDumont
			* 	Porém, quando são apenas dois nodes, esse problema nunca aconteceu nas inúmeras vezes que foi testado
			*/
			//execute until the penultimate one
			for (i=0; i<countCommands-1; i++){
				//printf("dentro da thread %d, processo %d, execute[%d][%d]: %s\n", iam, rank, position, iterator, executeMatrix[position][i]);
				sysReturn = system(executeMatrix[position][i]);
	/*			if(sysReturn!=0){
				    errnum=errno;
				    //current_time2 = time(NULL);
				    //fprintf(stderr, "tempo inicial: %d, tempo final: %d, thread %d, rank %d\n", current_time, current_time2, iam, rank);
				    fprintf(stderr,"error %s (%d) status system = %d, rank %d\n=> %s\n", strerror(errnum), errnum, sysReturn, rank, executeMatrix[position][i]);
                    ///codigo engessado -> refinar*/
                    //se o erro for no fasterq-dump:
                    if (position == 1){
                        strcat(buffer_rewind, cmd_rm);
                        strcat(buffer_rewind, path);
                        strcat(buffer_rewind, samplesMatrix[position]);
                        strcat(buffer_rewind, content_rm_fastq);
                        //sysReturn = system(buffer_rewind);
                        print("rewind: %s\n", buffer_rewind);
                    }
                    if (position == 2){
                        strcat(buffer_rewind, cmd_rm);
                        strcat(buffer_rewind, path);
                        strcat(buffer_rewind, samplesMatrix[position]);
                        strcat(buffer_rewind, content_rm_gz_1);
                        //sysReturn = system(buffer_rewind);
                        print("rewind: %s\n", buffer_rewind);
                    }
                    if (position == 3){
                        strcat(buffer_rewind, cmd_rm);
                        strcat(buffer_rewind, path);
                        strcat(buffer_rewind, samplesMatrix[position]);
                        strcat(buffer_rewind, content_rm_gz_2);
                        //sysReturn = system(buffer_rewind);
                        print("rewind: %s\n", buffer_rewind);
                    }
				//}
			}


			//execute the last one
			sysReturn = system(executeMatrix[position][i]);
			if(sysReturn!=0){
			        fprintf(stderr,"error in rank %d %s\n", rank, executeMatrix[position][i]);
			}
			current_time = time(NULL);
			duration = difftime(current_time, initial_time);
			fprintf(stderr,"termino no rank %d, duration:  %lf, %s\n", rank, duration, executeMatrix[position][0]);                
			iterator++;
			position = rank + numprocs*iterator;
			
		}



    /***********************************************************************************/

  MPI_Finalize();
  return 0;
}



