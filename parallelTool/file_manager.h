# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# define MAX_BUFFER_CHAR 512
# define MAX_LINE 512

//Estruturas de dados
struct queueNode{
	char nome[50];
	struct queueNode *next;
} queueNode;


//declaracoes de funcoes
 char* readlinefile (FILE *, int, char *);
// struct queueNode insertElem (struct queueNode *, char *);
// struct queueNode* trataSamples (FILE *, struct queueNode *);
void printQueue(struct queueNode *queue);
struct queueNode* retornaElemN(struct queueNode *, int );
void insertVariableValue(char *, char *, char *);
void makeQueueOutOfCommandsAndSample(char **, int , char *, char** );