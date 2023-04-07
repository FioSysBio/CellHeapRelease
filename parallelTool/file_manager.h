# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# define MAX_BUFFER_CHAR 512
# define MAX_LINE 512

//data structs
struct queueNode{
	char nome[50];
	struct queueNode *next;
} queueNode;


//functions declarations
 char* readlinefile (FILE *, int, char *);
void printQueue(struct queueNode *queue);
struct queueNode* retornaElemN(struct queueNode *, int );
void insertVariableValue(char *, char *, char *);
void makeQueueOutOfCommandsAndSample(char **, int , char *, char** );
