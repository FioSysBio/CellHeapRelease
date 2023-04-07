# include "file_manager.h"
//definicoes das funcoes

char* readlinefile (FILE *f, int length, char *c){
	//char c;
	char *string, *aux;
	int i=0;
	//buffer de tamanho excedente
	printf("antes do malloc de readlinefile\n");
	string = (char*)calloc(MAX_BUFFER_CHAR, sizeof(char));

	if(string==NULL)
		printf("readlinefile falhou no calloc\n");
	printf("depois do malloc de readlinefile\n");

	*c = fgetc(f);
	printf("depois do primeiro fget de readlinefile\n");
    while(*c!='\n' && *c!='\r' && *c!= EOF){
            string[i]=*c;
            i++;
			*c = fgetc(f);
			printf("em cada fgetc do while: %c\n", *c);
	}
	string[i]='\0';

	printf("q q tem depois do while: %s\n", string);
	//aux de tamanho apropriado para que nao haja lixo
	//aux = (char*)malloc(sizeof(char)*(length));
	//strncpy(aux, string, length);
	return string;
}




void printQueue(struct queueNode *queue){
	struct queueNode *aux;
	aux=queue;
	while(aux!=NULL){
		printf("%s\n", aux->nome);
		aux=aux->next;
	}
}

struct queueNode* retornaElemN(struct queueNode *queue, int n){
	int i=0;
	struct queueNode *aux = queue;
	while(i<n && aux!=NULL){
		aux=aux->next;
		i++;
	}
	return aux;
 }

//insert a string value when it finds the symbol $ in string sentence
void insertVariableValue(char *sentence, char *value, char *result){
	int sentenceSize = strlen(sentence);
	int valueSize = strlen(value);
	int i = 0;
	char j = 0;
	int k;
	char c='-';
	char *aux;
	char aux2[MAX_BUFFER_CHAR];
	char aux3[MAX_BUFFER_CHAR];
	char *aux4;
	for(k= 0; k<MAX_BUFFER_CHAR; k++){
		aux2[k]='\0';
	}

	strcpy(aux3, sentence);
	aux4 = strtok (aux3,"$");
	while (aux4 != NULL)
		{
			aux4 = strtok (NULL, "$");
			i++;
		}

	aux = strtok (sentence,"$");
	while (i>1)
		{
			strcat(aux2, aux);
			strcat(aux2, value);
			//printf ("%s\n",aux2);
			aux = strtok (NULL, "$");
			i--;
		}

	strcat(aux2, aux);
	strcpy(result, aux2);
}



void makeQueueOutOfCommandsAndSample(char **commandsMatrix, int rows, char* sampleElem, char** result){
	//char aux[rows][MAX_BUFFER_CHAR];
	char **aux;
	int i = 0;
	int j, k;
	char cpyCommands[rows][MAX_BUFFER_CHAR];
	char bufferResult[rows][MAX_BUFFER_CHAR];
	memcpy(cpyCommands, commandsMatrix, (rows*MAX_BUFFER_CHAR));
	char bufCom[MAX_BUFFER_CHAR];
	char bufRes[MAX_BUFFER_CHAR];

	//let's clean everybody

	for(j=0; j< rows; j++){
			for(k=0; k< MAX_BUFFER_CHAR; k++){
				bufferResult[j][k]='\0';
		}
	}

	while(i<rows){
		for(j=0; j< MAX_BUFFER_CHAR; j++){
			bufCom[j]='\0';
			bufRes[j]='\0';
		}
		//printf("makeQueue antes do memcpy cpyCommands: %s, bufCom: %s\n", cpyCommands[i], bufCom);
		memcpy(bufCom, cpyCommands[i], MAX_BUFFER_CHAR);
		//printf("makeQueue depois do memcpy cpyCommands: %s, bufCom: %s\n", cpyCommands[i], bufCom);
		insertVariableValue(bufCom, sampleElem, bufRes);
		//printf("makeQueue depois do insertVariable bufCom: %s, bufRes: %s\n", bufCom, bufRes);
		memcpy(bufferResult[i], bufRes, MAX_BUFFER_CHAR);
		//printf("makeQueue depois do outro memcpy bufferResult[i]: %s, bufRes: %s\n", bufferResult[i], bufRes);
		i++;
	}
	memcpy(result, bufferResult, rows*MAX_BUFFER_CHAR);
	//printf("result in makequeue: %s\n", )
	
}

