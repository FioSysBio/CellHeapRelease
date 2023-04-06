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

//void insertElem (struct queueNode *, char *);
// struct queueNode * insertElem (struct queueNode * queue, char * buffer){

// 	//create the new Element
// 	struct queueNode new;
// 	struct queueNode *aux;

// 	//new = (struct queueNode*) malloc (sizeof(struct queueNode));
// 	strcpy(new.nome, buffer);
// 	new.next=NULL;

// 	//not the first element ever
// 	if(queue!=NULL){
// 		aux=queue;
// 		while(aux->next!=NULL){
// 			aux=aux->next;
// 		}

// 		//update before the last one
// 		aux->next=new;

// 	//the first element ever
// 		return queue;
// 	}else{
// 		return &new;
// 	}
// }

// struct queueNode * trataSamples (FILE *f, struct queueNode *queue){

// 	char c;
// 	char *buffer;
// 	int i=0;
// 	printf("chegou em  trataSample\n");
// 	//buffer = (char*)malloc(sizeof(char)*(MAX_BUFFER_CHAR));
// 	buffer = (char*)malloc(MAX_LINE);
// 	if(buffer==NULL)
// 		printf("impossivel alocar\n");
// 	printf("malocou buffer de trataSample\n");

// 	while(c!=EOF){ 
// 		printf("antes do da chamada de readline de trataSamples\n");
// 		//buffer=readlinefile(f, MAX_BUFFER_CHAR, &c);
// 		fscanf(f, "%[^\b]", buffer);
// 		printf("depois do da chamada de readline de trataSamples\n");
// 		queue = insertElem (queue,buffer);
// 		printf("depois do da chamada de insertElem de trataSamples\n");
// 		for (int i = 0; i < MAX_BUFFER_CHAR; i++)
// 		{
// 			buffer='\0';
// 		}	

// 		printf("depois de limpar o buffer no trataSample\n");
// 	}
// 	return queue;
// }

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

//insere a string value onde encontrar o simbolo $ na string sentence
// char * insertVariableValue(char *sentence, char *value){
// 	int sentenceSize = strlen(sentence);
// 	int valueSize = strlen(value);
// 	int i = 0;
// 	char j = 0;
// 	char c='-';
// 	char *buffer = (char*) malloc (sizeof(char)*MAX_BUFFER_CHAR);
// 	while(c!='\0'){
// 		c=sentence[i];
// 		if (c=='$'){
// 			strcat(buffer, value);
// 			j=j+valueSize;
// 			i++;
// 		}else{
// 			buffer[j]=c;
// 			i++;
// 			j++;
// 		}
// 		if (i>=MAX_BUFFER_CHAR){
// 			printf("extrapolou tamanho da string\n");
// 			return NULL;
// 		}		
// 	}
// 	return buffer;
// }

//insere a string value onde encontrar o simbolo $ na string sentence
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
	//strcpy(aux, sentence);
	// printf("entrada do insertVariable: %s\n", sentence);
	// while(c!='\0'){
	// 	c=sentence[i];
	// 	if (c=='$'){
	// 		strcat(aux, value);
	// 		//strncpy(aux, value, valueSize-1);
	// 		j=j+valueSize;
	// 		i++;
	// 		printf("aux quando $: -%s-\n", aux);
	// 	}else{
	// 		aux[j]=c;
	// 		i++;
	// 		j++;
	// 		printf("aux normal: -%s-\n", aux);
	// 	}
	// 	if (i>=MAX_BUFFER_CHAR){
	// 		printf("extrapolou tamanho da string\n");
	// 		break;
	// 	}		
	// }
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
	//printf("aux: %s, result: %s\n\n", aux, result);
}

// struct queueNode * makeQueueOutOfCommandsAndSample(struct queueNode *commandsQueue , struct queueNode * samplesQueue){
//     struct queueNode *commands = commandsQueue;
// 	struct queueNode *sample = samplesQueue;
// 	struct queueNode *aux = NULL;
// 	char *buffer = (char *) malloc (sizeof(char)*MAX_BUFFER_CHAR);

// 	while(commands!=NULL){
// 		buffer=insertVariableValue(commands->nome, sample->nome);
// 		aux= insertElem(aux, buffer);
// 		commands=commands->next;
// 	}
// 	return aux;
// }

// char ** makeQueueOutOfCommandsAndSample(char **commandsMatrix, int rows, char* sampleElem){
// 	//char aux[rows][MAX_BUFFER_CHAR];
// 	char **aux;
// 	char *buffer = (char *) malloc (sizeof(char)*MAX_BUFFER_CHAR);
// 	int i = 0;
// 	aux=(char**) malloc (sizeof(char*)*rows);
// 	while(i<rows){
// 		buffer=insertVariableValue(commandsMatrix[i], sampleElem);
// 		aux[i] = (char*) malloc (sizeof(char)*MAX_BUFFER_CHAR);
// 		strcpy(aux[i], buffer);
// 		i++;
// 	}
// 	return aux;
// }

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
	//printf("resultado no makequeue: %s\n", )
	
}

// void makeQueueOutOfCommandsAndSample(char **commandsMatrix, int rows, char* sampleElem, char** result){
// 	//char aux[rows][MAX_BUFFER_CHAR];
// 	int i = 0;
// 	int j;
// 	// char cpyCommands[rows][MAX_BUFFER_CHAR];
// 	// char bufferResult[rows][MAX_BUFFER_CHAR];
// 	// memcpy(cpyCommands, commandsMatrix, (rows*MAX_BUFFER_CHAR));
// 	// char bufCom[MAX_BUFFER_CHAR];
// 	char bufRes[MAX_BUFFER_CHAR];

// 	while(i<rows){
// 		for(j=0; j< MAX_BUFFER_CHAR; j++){
// 			bufRes[j]='\0';
// 		}
// 		printf("commandsMatrix[%d]: %s\n", i, (commandsMatrix+i*MAX_BUFFER_CHAR));
// 		insertVariableValue((commandsMatrix+i*MAX_BUFFER_CHAR), sampleElem, bufRes);
// 		printf("resultado no makequeue: %s\n", bufRes );
// 		memcpy((result+i*MAX_BUFFER_CHAR), bufRes, MAX_BUFFER_CHAR);
// 		i++;
// 	}
// 	//printf("resultado no makequeue: %s\n", )
	
// }