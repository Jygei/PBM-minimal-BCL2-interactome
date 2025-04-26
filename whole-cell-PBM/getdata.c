#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ReadData.h"
#include "WriteData.h"



int main( int argc, char **argv )
{
    FILE *file;
    /* Check the Number of Arguments is Correct */
    if (argc != 5)
    {
        printf("Wrong number of arguments (./getdata [Number of Iterations] [Iteration Step Size] [xml directory] [excelfile.xls]\n");
        exit (0);
    }
    char * filepath;
	char * excel_filepath;
    
    xmachine_Receptor * Receptors, * current_Receptor;
    xmachine_Receptor ** p_Receptors;
    xmachine_Ligand * Ligands, * current_Ligand;
    xmachine_Ligand ** p_Ligands;
	xmachine_Mito * Mitos, * current_Mito;
    xmachine_Mito ** p_Mitos;
    
    /* Initialise Variables with Command Line Arguments */
    int current_iteration = 0;
    int total_iterations = atoi(argv[1]);
    int iteration_step = atoi(argv[2]);
    filepath = argv[3];
	excel_filepath = argv[4];
    
    int stilldata = 1;
    
    printf("Output dir: %s\n", filepath);
    
    /* Initialise pointers */
    Receptors = NULL;
    p_Receptors = &Receptors;
    Ligands = NULL;
    p_Ligands = &Ligands;
	Mitos = NULL;
    p_Mitos = &Mitos;
    
    while(current_iteration <= total_iterations)
    {
        /* Read Iteration Data */
        stilldata = getiteration(filepath, current_iteration, p_Receptors, p_Ligands, p_Mitos);
        
        /* Write Iteration Data */
        if(stilldata) savedatatofile(current_iteration, p_Receptors, p_Ligands, p_Mitos, excel_filepath);
        current_iteration += iteration_step;
    }
    
    /* Should never get here */
    return( 0 );
}
