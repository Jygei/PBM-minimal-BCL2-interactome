/*************************************
 * Struct Type Definitions            *
 *************************************/

struct xmachine_Receptor
{
    int state;
    double x;
    double y;
    double z;
	double collision_counter_ligand;
	int mitosize;
	int mitoid;

    struct xmachine_Receptor * next;
};

typedef struct xmachine_Receptor xmachine_Receptor;

struct xmachine_Ligand
{
    int state;
    double x;
    double y;
    double z;
	double collision_counter_ligand;
    
    struct xmachine_Ligand * next;
};

typedef struct xmachine_Ligand xmachine_Ligand;

struct xmachine_Mito
{
    int mito_collision_counter;
	int mitosize;
    
    struct xmachine_Mito * next;
};

typedef struct xmachine_Mito xmachine_Mito;

/*****************************************************************
 * FUNCTIONS: linked list functions                               *
 * PURPOSE: to allocate and free memory in linked lists           *
 *****************************************************************/
void freeReceptors(xmachine_Receptor * head)
{
    /* Tempory element needed for loop */
    xmachine_Receptor * tmp;
    
    /* Loop while head is not NULL */
    while(head)
    {
        /* Make next in list tmp */
        tmp = head->next;
        /* Free memory of head */
        free(head);
        /* Make head the next in the list */
        head = tmp;
    }
}

void freeLigands(xmachine_Ligand * head)
{
    /* Tempory element needed for loop */
    xmachine_Ligand * tmp;
    
    /* Loop while head is not NULL */
    while(head)
    {
        /* Make next in list tmp */
        tmp = head->next;
        /* Free memory of head */
        free(head);
        /* Make head the next in the list */
        head = tmp;
    }
}

void freeMitos(xmachine_Mito * head)
{
    /* Tempory element needed for loop */
    xmachine_Mito * tmp;
    
    /* Loop while head is not NULL */
    while(head)
    {
        /* Make next in list tmp */
        tmp = head->next;
        /* Free memory of head */
        free(head);
        /* Make head the next in the list */
        head = tmp;
    }
}


xmachine_Receptor * addReceptor(xmachine_Receptor ** pointer_to_Receptors, xmachine_Receptor * current)
{
    /* The new tail of the linked list */
    xmachine_Receptor * tail;
    
    /* Allocate memory for new neighbour data */
    
    tail = (xmachine_Receptor *)malloc(sizeof(xmachine_Receptor));
    
    /* Check if current is not NULL */
    if(current)
    {
        /* Current exists therefore make its next point to tail */
        current->next = tail;
    }
    else
    {
        /* Current is NULL therefore make the cell neighbour_head point to tail */
        *pointer_to_Receptors = tail;
    }
    /* Point next to NULL */
    tail->next = NULL;
    /* Return new neighbour data */
    return tail;
}


xmachine_Ligand * addLigand(xmachine_Ligand ** pointer_to_Ligands, xmachine_Ligand * current)
{
    /* The new tail of the linked list */
    xmachine_Ligand * tail;
    
    /* Allocate memory for new neighbour data */
    
    tail = (xmachine_Ligand *)malloc(sizeof(xmachine_Ligand));
    
    /* Check if current is not NULL */
    if(current)
    {
        /* Current exists therefore make its next point to tail */
        current->next = tail;
    }
    else
    {
        /* Current is NULL therefore make the cell neighbour_head point to tail */
        *pointer_to_Ligands = tail;
    }
    /* Point next to NULL */
    tail->next = NULL;
    /* Return new neighbour data */
    return tail;
}

xmachine_Mito * addMito(xmachine_Mito ** pointer_to_Mitos, xmachine_Mito * current)
{
    /* The new tail of the linked list */
    xmachine_Mito * tail;
    
    /* Allocate memory for new neighbour data */
    
    tail = (xmachine_Mito *)malloc(sizeof(xmachine_Mito));
    
    /* Check if current is not NULL */
    if(current)
    {
        /* Current exists therefore make its next point to tail */
        current->next = tail;
    }
    else
    {
        /* Current is NULL therefore make the cell neighbour_head point to tail */
        *pointer_to_Mitos = tail;
    }
    /* Point next to NULL */
    tail->next = NULL;
    /* Return new neighbour data */
    return tail;
}


/*****************************************************************
 * FUNCTION: getIteration                                         *
 * PURPOSE: read iteration xml file                               *
 *****************************************************************/
int getiteration(char * filepath, int itno, xmachine_Receptor ** pointer_to_Receptors , xmachine_Ligand ** pointer_to_Ligands, xmachine_Mito ** pointer_to_Mitos)
{
    /* Pointer to file */
    FILE *file;
    /* Char and char buffer for reading file to */
    char c = ' ';
    char buffer[30];
    /* Variable to keep reading file */
    int reading = 1;
    /* Variables for checking tags */
    int i, initeration, inagent;
    int intag, instate, inname, inSd, inx, iny, inz, inmito_collision_counter, incurrent_mitosize, incollision_counter_ligand, inmitosize, inmitoid;
    /* Variables for model data */
    int state, mito_collision_counter, current_mitosize, collision_counter_ligand, mitosize, mitoid;
    double Sd, x, y, z;
    char name[100];
    
    xmachine_Receptor * current_Receptor, * tail_Receptor;
    tail_Receptor = *pointer_to_Receptors;
    current_Receptor = NULL;
    
    xmachine_Ligand * current_Ligand, * tail_Ligand;
    tail_Ligand = *pointer_to_Ligands;
    current_Ligand = NULL;
	
	xmachine_Mito * current_Mito, * tail_Mito;
    tail_Mito = *pointer_to_Mitos;
    current_Mito = NULL;
    
    /* Open config file to read-only */
    char data[200];
    sprintf(data, "%s%i%s", filepath, itno, ".xml");
    printf("%s", data);
    file = fopen("getdatalog.txt", "a");
    fputs(data, file);
    fputs("\n", file);
    fclose (file);
    if((file = fopen(data, "r"))==NULL)
    {
        printf(" nent\n");
        return 0;
        /*exit(0);*/
    }
    else
    {
        
        printf("\n");
        /* Initialise variables */
        i = 0;
        intag = 0;
        initeration = 0;
        inagent = 0;
        instate = 0;
        inname = 0;
        state = 0;
        inx = 0;
        iny = 0;
        inz = 0;
        inx = 0;
        iny = 0;
        inz = 0;
		inmito_collision_counter = 0;
		mito_collision_counter = 0;
		incurrent_mitosize = 0;
		current_mitosize = 0;
		incollision_counter_ligand = 0;
		collision_counter_ligand = 0;
		inmitosize = 0;
		mitosize = 0;
		inmitoid = 0;
		mitoid = 0;
        
        /* Read characters until end of xml found */
        while(reading == 1)
        {
            /* Get the next char from the file */
            c = fgetc(file);
            
            /* If the end of a tag */
            if(c == '>')
            {
                /* Place 0 at end of buffer to make chars a string */
                buffer[i] = 0;
                if(strcmp(buffer, "states") == 0) { initeration = 1; }
                if(strcmp(buffer, "/states") == 0) { initeration = 0; reading = 0; }
                if(strcmp(buffer, "xagent") == 0) { inagent = 1; }
                if(strcmp(buffer, "/xagent") == 0)
                {
                    inagent = 0;
                    
                    //printf("Found agent ... ");
                    
                    if(strcmp(name, "Receptor") == 0)
                    {
                        //printf("Adding agent\n");
                        
                        /* check if tail is NULL */
                        if(tail_Receptor == NULL)
                        {
                            //printf("tail is null allocate more memory\n");
                            /* Allocate memory */
                            tail_Receptor = addReceptor(pointer_to_Receptors, current_Receptor);
                        }
                        //else printf("tail exisits\n");
                        
                        /* Make tail the current element to add to */
                        current_Receptor = tail_Receptor;
                        
                        current_Receptor->x = x;
                        current_Receptor->y = y;
                        current_Receptor->z = z;
                        current_Receptor->state = state;
						current_Receptor->mitosize = mitosize;
						current_Receptor->mitoid = mitoid;
                        
                        /* Make tail the next element in the linked list */
                        tail_Receptor = current_Receptor->next;
                    }
                    
                    if(strcmp(name, "Ligand") == 0)
                    {
                        //printf("Adding agent\n");
                        
                        /* check if tail is NULL */
                        if(tail_Ligand == NULL)
                        {
                            //printf("tail is null allocate more memory\n");
                            /* Allocate memory */
                            tail_Ligand = addLigand(pointer_to_Ligands, current_Ligand);
                        }
                        //else printf("tail exisits\n");
                        
                        /* Make tail the current element to add to */
                        current_Ligand = tail_Ligand;
                        
                        current_Ligand->state = state;
                        current_Ligand->x = x;
                        current_Ligand->y = y;
                        current_Ligand->z = z;
						current_Ligand->collision_counter_ligand = collision_counter_ligand;
                        
                        /* Make tail the next element in the linked list */
                        tail_Ligand = current_Ligand->next;
                    }
					
					if(strcmp(name, "Mito") == 0)
                    {
                        //printf("Adding agent\n");
                        
                        /* check if tail is NULL */
                        if(tail_Mito == NULL)
                        {
                            //printf("tail is null allocate more memory\n");
                            /* Allocate memory */
                            tail_Mito = addMito(pointer_to_Mitos, current_Mito);
                        }
                        //else printf("tail exisits\n");
                        
                        /* Make tail the current element to add to */
                        current_Mito = tail_Mito;
                        
                        current_Mito->mito_collision_counter = mito_collision_counter;
						current_Mito->mitosize = current_mitosize;
                        
                        /* Make tail the next element in the linked list */
                        tail_Mito = current_Mito->next;
                    }
                    
                    
                    
                    //else printf("Not adding agent\n");
                }
                if(strcmp(buffer, "name") == 0) { inname = 1; }
                if(strcmp(buffer, "/name") == 0) { inname = 0; }
                if(strcmp(buffer, "state") == 0) { instate = 1; }
                if(strcmp(buffer, "/state") == 0) { instate = 0; }
                if(strcmp(buffer, "x") == 0) { inx = 1; }
                if(strcmp(buffer, "/x") == 0) { inx = 0; }
                if(strcmp(buffer, "y") == 0) { iny = 1; }
                if(strcmp(buffer, "/y") == 0) { iny = 0; }
                if(strcmp(buffer, "z") == 0) { inz = 1; }
                if(strcmp(buffer, "/z") == 0) { inz = 0; }
				if(strcmp(buffer, "mito_collision_counter") == 0) { inmito_collision_counter = 1; }
                if(strcmp(buffer, "/mito_collision_counter") == 0) { inmito_collision_counter = 0; }
				if(strcmp(buffer, "current_mitosize") == 0) { incurrent_mitosize = 1; }
                if(strcmp(buffer, "/current_mitosize") == 0) { incurrent_mitosize = 0; }
				if(strcmp(buffer, "collision_counter_ligand") == 0) { incollision_counter_ligand = 1; }
                if(strcmp(buffer, "/collision_counter_ligand") == 0) { incollision_counter_ligand = 0; }
				if(strcmp(buffer, "mitosize") == 0) { inmitosize = 1; }
                if(strcmp(buffer, "/mitosize") == 0) { inmitosize = 0; }
				if(strcmp(buffer, "mitoid") == 0) { inmitoid = 1; }
                if(strcmp(buffer, "/mitoid") == 0) { inmitoid = 0; }
				
                /* End of tag and reset buffer */
                intag = 0;
                i = 0;
            }
            /* If start of tag */
            else if(c == '<')
            {
                /* Place /0 at end of buffer to end numbers */
                buffer[i] = 0;
                /* Flag in tag */
                intag = 1;
                /* If just read data into buffer retrieve it */
                if(inagent && inname)  { strcpy(name, buffer); }
                if(inagent && instate)  { state = atoi(buffer); }
                if(inagent && inx)  { x = atof(buffer); }
                if(inagent && iny)  { y = atof(buffer); }
                if(inagent && inz)  { z = atof(buffer); }
				if(inagent && inmito_collision_counter)  { mito_collision_counter = atoi(buffer); }
				if(inagent && incurrent_mitosize)  { current_mitosize = atoi(buffer); }
				if(inagent && incollision_counter_ligand)  { collision_counter_ligand = atoi(buffer); }
				if(inagent && inmitosize)  { mitosize = atoi(buffer); }
				if(inagent && inmitoid)  { mitoid = atoi(buffer); }
                
                
                /* Reset buffer */
                i = 0;
            }
            /* If in tag put read char into buffer */
            else if(intag)
            {
                buffer[i] = c;
                i++;
            }
            /* If in data read char into buffer */
            else if(inname || instate || inx || iny || inz || inmito_collision_counter || incurrent_mitosize || incollision_counter_ligand || inmitosize || inmitoid)
            {
                buffer[i] = c;
                i++;
            }
        }
        
        /* Free memory for unused linked list elements */
        if(tail_Receptor)
        {
            freeReceptors(tail_Receptor);
            /* Make pointer to tail equal NULL */
            if(current_Receptor) { current_Receptor->next = NULL; }
        }
        
        /* Free memory for unused linked list elements */
        if(tail_Ligand)
        {
            freeLigands(tail_Ligand);
            /* Make pointer to tail equal NULL */
            if(current_Ligand) { current_Ligand->next = NULL; }
        }
        
		/* Free memory for unused linked list elements */
        if(tail_Mito)
        {
            freeMitos(tail_Mito);
            /* Make pointer to tail equal NULL */
            if(current_Mito) { current_Mito->next = NULL; }
        }
		
        /* Close the file */
        fclose(file);
        
        
        return 1;
    }
}
