void savedatatofile(int itno, xmachine_Receptor ** pointer_to_Receptors, xmachine_Ligand ** pointer_to_Ligands, xmachine_Mito ** pointer_to_Mitos, char * excel_filepath) //Changed by Gavin
{
FILE *file;
xmachine_Receptor * current_Receptor;
current_Receptor = *pointer_to_Receptors;

xmachine_Ligand * current_Ligand;
current_Ligand = *pointer_to_Ligands;

xmachine_Mito * current_Mito;
current_Mito = *pointer_to_Mitos;

/* Declare Initial Variables */
int total = 0;
int totalReceptor = 0;
int totalLigand = 0;

int State_1 = 0;
int State_2 = 0;
int State_20 = 0;

//int mito_collision_counter_total = 0;
int mito_number = 198;										/* has to be changed to actual mitochondria number */
double collision_counter_ligand_total = 0;
double collisions_total = 0;
double small_mito_bound = 0;
double normal_mito_bound = 0;
int current_mitonum = 0;
int mitoid[198];											/* has to be changed to actual mitochondria number */
for (int i = 0; i < mito_number; i++)	{	mitoid[i]=0;	}
int current_mitosize[198];									/* has to be changed to actual mitochondria number */
for (int i = 0; i < mito_number; i++)	{	current_mitosize[i]=2;	}

if (itno == 0)
{
    file = fopen(excel_filepath, "w");
    fprintf(file, "Time(s)\t");
    fprintf(file, "MCL-1 Cyto\t");
    fprintf(file, "MCL-1 Mito\t");
//	fprintf(file, "mito_collision_counter\t");
//	fprintf(file, "collision_counter_ligand\t");
	fprintf(file, "collisions_total\t");
	fprintf(file, "small_mito_bound\t");
	fprintf(file, "normal_mito_bound\t");
	for(int i = 0; i < mito_number; i++)
	{ fprintf(file, "mito ID%d\t ", i); }
	
    fprintf(file, "\n");
    fclose(file);
}

current_Receptor = *pointer_to_Receptors;
while(current_Receptor)
{
    total++;
    totalReceptor++;
	collisions_total += current_Receptor->collision_counter_ligand;
	if (current_Receptor->mitosize == 0)
	{
		small_mito_bound++;
	}
	if (current_Receptor->mitosize == 1)
	{
		normal_mito_bound++;
	}

	for (int i = 0; i < mito_number; i++)
	{
		if (current_Receptor->mitoid == i)
		{
			mitoid[i]++;
		}
	}

    if (current_Receptor->state == 20)
    {
        State_20++;
    }

    current_Receptor = current_Receptor->next;
}

current_Ligand = *pointer_to_Ligands;
while(current_Ligand)
{
    total++;
    totalLigand++;
//	collision_counter_ligand_total += current_Ligand->collision_counter_ligand;
	collisions_total += current_Ligand->collision_counter_ligand;

    if (current_Ligand->state == 1)
    {
        State_1++;
    }
    if (current_Ligand->state == 2)
    {
        State_2++;
    }
    current_Ligand = current_Ligand->next;
}

current_Mito = *pointer_to_Mitos;
while(current_Mito)
{
	//for (int i = 0; i < mito_number; i++)
	//{
		if (current_Mito->mitosize == 1)
		{
			current_mitosize[current_mitonum] = 1;
		}
		if (current_Mito->mitosize == 0)
		{
			current_mitosize[current_mitonum] = 0;
		}
	//}
	
	//mito_collision_counter_total = current_Mito->mito_collision_counter;
	current_mitonum++;
    current_Mito = current_Mito->next;
 }
 
if (itno == 0)
{
	file = fopen(excel_filepath, "a");
	for (int i = 0; i < 6; i++)
	{ fprintf(file, "\t"); }
	for (int i = 0; i < mito_number; i++)
	{ fprintf(file, "%d\t", current_mitosize[i]); }
	fprintf(file, "\n");
	fclose(file);
}

file = fopen(excel_filepath, "a");
    fprintf(file, "%f\t", itno*0.050000);
    fprintf(file, "%d\t", State_2);
    fprintf(file, "%d\t", State_20);
//	fprintf(file, "%d\t", mito_collision_counter_total);
//	fprintf(file, "%f\t", collision_counter_ligand_total);     
	fprintf(file, "%f\t", collisions_total);
	fprintf(file, "%f\t", small_mito_bound);
	fprintf(file, "%f\t", normal_mito_bound);
	for (int i = 0; i < mito_number; i++)
	{ fprintf(file, "%d\t", mitoid[i]); }
    fprintf(file, "\n");

fclose(file);
}
