/*

Jenny Geiger
Institute of Cell Biology and Immunology (IZI) at the University of Stuttgart

FaST generated Codes were changed and modified to simulate a single mitochondrial environment with protein interactions
of a reduced interactome of the BCL-2 family (MCL-1, BAK, tBID)


Code version information: this code creates the 0.xml file
- code of creating other members of the BCL-2 protein family are remains of the FaST generated code and was not removed for completeness
- changing the seed of the random number generator accordingly to the wished mitochondria positions: srand(time(NULL)) -> random positions; srand(1) -> not at random positions
- change in Mito creation the minimum distance to the cell boundaries depends on the Dc_state

*/


#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "AgentVariables.h"
#include "globals.h"


int main()
{
FILE *file;

/* Seed the random number generator */
srand (3);								/************************************** ATTENTION: change between srand(time(NULL)) -> mitos at random positions or srand(1) -> mitos at the same position **************************************/

/* Initialise the 0.xml file */
printf ("Initialising 0.xml File\n");

file = fopen("xml/0.xml", "w");
fputs("<states>\n", file);
fputs("<itno>0</itno>\n", file);
fputs("<environment>\n", file);
fputs("</environment>\n", file);

/* Declare agent counters for id */
int Receptor_Agents = 0;
int Ligand_Agents = 0;
int Mito_Agents = 0;



/********************************************************************************************************/
/********************************************* Mitochondria *********************************************/
/********************************************************************************************************/
/* Create Mitos */
printf ("Creating Mito Ligand Agents\n");
int counter = 1;
int mito_counter = 0;

// mito sizes in G1 phase of fragmented mitos
double length, height, width, radius, radius2;
length = 580.003;		//580.003;		//in nm
radius = length/2;
height = 580.003;		//580.003;		//in nm
width = 25428.1;		//25428.1;		//in nm
radius2 = width/2;

// new total mito_number because fragmentation changes mito_number
int mito_normal_number = mito_mean_number*(1.0-fragmentation_percent);

// matrices of x-, y- and z-coordinates of the mitos; flexible row-number which is the same as the total number of mitos; variables are set free in the end of Mito creation
// x
double **mitox = (double **)malloc(mito_number * sizeof(double *));
    for (int i = 0; i < mito_number; i++)
	{
        mitox[i] = (double *)malloc(8 * sizeof(double));  // 8 columns, mito_number rows
    }
// y 
double **mitoy = (double **)malloc(mito_number * sizeof(double *));
    for (int i = 0; i < mito_number; i++)
	{
        mitoy[i] = (double *)malloc(8 * sizeof(double));  // 8 columns, mito_number rows
    }	
// z
double **mitoz = (double **)malloc(mito_number * sizeof(double *));
    for (int i = 0; i < mito_number; i++)
	{
        mitoz[i] = (double *)malloc(8 * sizeof(double));  // 8 columns, mito_number rows
    }

// defining lists for all mitochondria variables with size of mito_number; variables are set free in the end of Mito creation
double *mitoradius = (double *)malloc(mito_number * sizeof(double));
double *mitoradius2 = (double *)malloc(mito_number * sizeof(double));
double *mitolength = (double *)malloc(mito_number * sizeof(double));
double *mitoheight = (double *)malloc(mito_number * sizeof(double));
double *mitowidth = (double *)malloc(mito_number * sizeof(double));
double *mitocentrex = (double *)malloc(mito_number * sizeof(double));
double *mitocentrey = (double *)malloc(mito_number * sizeof(double));
double *mitocentrez = (double *)malloc(mito_number * sizeof(double));
double *mito_collision_counter = (double *)malloc(mito_number * sizeof(double));
double *mito_energetics = (double *)malloc(mito_number * sizeof(double));
int *current_mitosize = (int *)malloc(mito_number * sizeof(int));

// variable for the biggest diffusion distance that is possible, it is necessary between mitos and the boundaries
// in this case its created for MCL-1 distance to have the exact same Mito arrangement compared to 24_019 and 24_021
double diffusion_distance = Dc_state_2*3;

// ID-list of mitos
int small_mito_check;		// 0 -> small_mito; 1 -> normal_mito
int *normal_mitos = NULL;
int *small_mitos = NULL;

// fills list with IDs of normal mitos (=1) -> normal mitochondria are created first (!)
while (mito_counter < mito_normal_number)
	{
		normal_mitos = (int *)realloc(normal_mitos, (mito_counter + 1) * sizeof(int));
		normal_mitos[mito_counter] = mito_counter;
		current_mitosize[mito_counter] = 1;
		mito_counter++;
	}
// fills list with IDs of fragmented mitos (=0)
while ((mito_counter >= mito_normal_number) && (mito_counter < mito_number))
	{
		small_mitos = (int *)realloc(small_mitos, (mito_counter + 1) * sizeof(int));
		small_mitos[mito_counter] = mito_counter;
		current_mitosize[mito_counter] = 0;
		mito_counter++;
	}

// Mito size, start coordinates, random distribution w/o intersection
mito_counter = 0;
	while (mito_counter < mito_number)
	{
		// mito sizes in G1 phase
		if (current_mitosize[mito_counter] == 1)	// normal mito
		{
			mitolength[mito_counter] = length;
			mitoradius[mito_counter] = radius;
			mitoheight[mito_counter] = height;
			mitowidth[mito_counter] = width;
			mitoradius2[mito_counter] = radius2;
		}
		
		if (current_mitosize[mito_counter] == 0)	// small mito
		{
			mitolength[mito_counter] = length;
			mitoradius[mito_counter] = radius;
			mitoheight[mito_counter] = height;
			mitowidth[mito_counter] = width/fragmentation;
			mitoradius2[mito_counter] = radius2/fragmentation;
		}
		
		// check for intersections; intersection = 1, no intersection = 0
		int intersecting = 1;
		while (intersecting == 1)
		{
			counter = 0;
			intersecting = 0;
			
			// Set Mito start coordinates (random startpoints) with minimal distance of MCL-1 diffusion distance to boundaries
			mitox[mito_counter][0] = minintracellular_xboundary + diffusion_distance + (rand()/(double)(RAND_MAX)*(maxintracellular_xboundary - minintracellular_xboundary - 2*diffusion_distance - (mitolength[mito_counter])));
			mitoy[mito_counter][0] = minintracellular_yboundary + diffusion_distance + (rand()/(double)(RAND_MAX)*(maxintracellular_yboundary - minintracellular_xboundary - 2*diffusion_distance - (mitoheight[mito_counter])));
			mitoz[mito_counter][0] = minintracellular_zboundary + diffusion_distance + (rand()/(double)(RAND_MAX)*(maxintracellular_zboundary - minintracellular_xboundary - 2*diffusion_distance - (mitowidth[mito_counter])));
			
			// Set Mito centre points
			mitocentrex[mito_counter] = mitox[mito_counter][0] + mitoradius[mito_counter];
			mitocentrey[mito_counter] = mitoy[mito_counter][0] + mitoradius[mito_counter];
			mitocentrez[mito_counter] = mitoz[mito_counter][0] + mitoradius2[mito_counter];
			
			// Set Mito cornerpoints
			mitox[mito_counter][1] = mitox[mito_counter][0] + mitolength[mito_counter];
			mitoy[mito_counter][1] = mitoy[mito_counter][0];
			mitoz[mito_counter][1] = mitoz[mito_counter][0];
			mitox[mito_counter][2] = mitox[mito_counter][0];
			mitoy[mito_counter][2] = mitoy[mito_counter][0] + mitoheight[mito_counter];
			mitoz[mito_counter][2] = mitoz[mito_counter][0];
			mitox[mito_counter][3] = mitox[mito_counter][0] + mitolength[mito_counter];
			mitoy[mito_counter][3] = mitoy[mito_counter][0] + mitoheight[mito_counter];
			mitoz[mito_counter][3] = mitoz[mito_counter][0];
			mitox[mito_counter][4] = mitox[mito_counter][0];
			mitoy[mito_counter][4] = mitoy[mito_counter][0];
			mitox[mito_counter][5] = mitox[mito_counter][0] + mitolength[mito_counter];
			mitoy[mito_counter][5] = mitoy[mito_counter][0];
			mitox[mito_counter][6] = mitox[mito_counter][0];
			mitoy[mito_counter][6] = mitoy[mito_counter][0] + mitoheight[mito_counter];
			mitox[mito_counter][7] = mitox[mito_counter][0] + mitolength[mito_counter];
			mitoy[mito_counter][7] = mitoy[mito_counter][0] + mitoheight[mito_counter];
			
			// difference between small and normal mito size in mitowidth; sets z-coordinates and centre of Mito
			mitoz[mito_counter][4] = mitoz[mito_counter][0] + mitowidth[mito_counter];
			mitoz[mito_counter][5] = mitoz[mito_counter][0] + mitowidth[mito_counter];
			mitoz[mito_counter][6] = mitoz[mito_counter][0] + mitowidth[mito_counter];
			mitoz[mito_counter][7] = mitoz[mito_counter][0] + mitowidth[mito_counter];
			
			// Prevent mito-creation in the blackbox/nucleus; no minimal distance needed if ligand diffuses through blackbox, otherwise add +/- diffusion_distance -> if mito is in blackbox (intersection = 1) it generates new random start coordinates
			/*if (intersecting == 0)
			{
				if ((((mitox[mito_counter][0] >= minnucleus_xboundary) && (mitox[mito_counter][0] <= maxnucleus_xboundary))	
					&& ((mitoy[mito_counter][0] >= minnucleus_yboundary) && (mitoy[mito_counter][0] <= maxnucleus_yboundary))
					&& ((mitoz[mito_counter][0] >= minnucleus_zboundary) && (mitoz[mito_counter][0] <= maxnucleus_zboundary)))
					
					|| (((mitox[mito_counter][0] + mitolength[mito_counter] >= minnucleus_xboundary) && (mitox[mito_counter][0] + mitolength[mito_counter] <= maxnucleus_xboundary))
					&& ((mitoy[mito_counter][0] + mitoheight[mito_counter] >= minnucleus_yboundary) && (mitoy[mito_counter][0] + mitoheight[mito_counter] <= maxnucleus_yboundary))
					&& ((mitoz[mito_counter][0] + mitowidth[mito_counter] >= minnucleus_zboundary) && (mitoz[mito_counter][0] + mitowidth[mito_counter] <= maxnucleus_zboundary)))	)																															   
				{
					intersecting = 1;
				}
			}*/
			
			// Prevent intersection of current mito with other mitos with a minimal diffusion distance -> if there is a intersection (intersection = 1) it generates new random start coordinates
			while (intersecting == 0 && (counter < mito_counter))
			{
				// calculation of (always positive) distances between two mitos
				double distancex = mitocentrex[mito_counter] - mitocentrex[counter];
				double distancey = mitocentrey[mito_counter] - mitocentrey[counter];
				double distancez = mitocentrez[mito_counter] - mitocentrez[counter];
				double distancex_pos = distancex;
				double distancey_pos = distancey;
				double distancez_pos = distancez;
				if (distancex < 0)	{	distancex_pos = fabs(distancex); 	}
				if (distancey < 0)	{	distancey_pos = fabs(distancey);	}
				if (distancez < 0)	{	distancez_pos = fabs(distancez);	}
				double min_diff_distancex = mitoradius[mito_counter] + mitoradius[counter] + diffusion_distance;
				double min_diff_distancey = mitoradius[mito_counter] + mitoradius[counter] + diffusion_distance;
				double min_diff_distancez = mitoradius2[mito_counter] + mitoradius2[counter] + diffusion_distance;
				
				// intersection if distance between two mitos is to greater than the diffusion distance; to small -> intersection = 1, restarts mito creation
				if	(	(	distancex_pos < min_diff_distancex	)
					&& 	(	distancey_pos < min_diff_distancey	)
					&& 	(	distancez_pos < min_diff_distancez	)	)
				{
					intersecting = 1;
				}
				
				if (intersecting == 0)	{	counter++;	}
			}
		}
		mito_counter++;
	}

// creating & saving mitochondria data
mito_counter = 0;
while (mito_counter < mito_number)
{
	mito_collision_counter[mito_counter] = 0.0;
	mito_energetics[mito_counter] = 0.0;

	// Write mito variables to 0.xml file
	fprintf(file, "<xagent>\n");
	fprintf(file, "<name>Mito</name>\n");
	fprintf(file, "<id>%d</id>\n", mito_counter);
	fprintf(file, "<x1>%f</x1>\n", mitox[mito_counter][0]);
	fprintf(file, "<y1>%f</y1>\n", mitoy[mito_counter][0]);
	fprintf(file, "<z1>%f</z1>\n", mitoz[mito_counter][0]);
	fprintf(file, "<x2>%f</x2>\n", mitox[mito_counter][1]);
	fprintf(file, "<y2>%f</y2>\n", mitoy[mito_counter][1]);
	fprintf(file, "<z2>%f</z2>\n", mitoz[mito_counter][1]);
	fprintf(file, "<x3>%f</x3>\n", mitox[mito_counter][2]);
	fprintf(file, "<y3>%f</y3>\n", mitoy[mito_counter][2]);
	fprintf(file, "<z3>%f</z3>\n", mitoz[mito_counter][2]);
	fprintf(file, "<x4>%f</x4>\n", mitox[mito_counter][3]);
	fprintf(file, "<y4>%f</y4>\n", mitoy[mito_counter][3]);
	fprintf(file, "<z4>%f</z4>\n", mitoz[mito_counter][3]);
	fprintf(file, "<x5>%f</x5>\n", mitox[mito_counter][4]);
	fprintf(file, "<y5>%f</y5>\n", mitoy[mito_counter][4]);
	fprintf(file, "<z5>%f</z5>\n", mitoz[mito_counter][4]);
	fprintf(file, "<x6>%f</x6>\n", mitox[mito_counter][5]);
	fprintf(file, "<y6>%f</y6>\n", mitoy[mito_counter][5]);
	fprintf(file, "<z6>%f</z6>\n", mitoz[mito_counter][5]);
	fprintf(file, "<x7>%f</x7>\n", mitox[mito_counter][6]);
	fprintf(file, "<y7>%f</y7>\n", mitoy[mito_counter][6]);
	fprintf(file, "<z7>%f</z7>\n", mitoz[mito_counter][6]);
	fprintf(file, "<x8>%f</x8>\n", mitox[mito_counter][7]);
	fprintf(file, "<y8>%f</y8>\n", mitoy[mito_counter][7]);
	fprintf(file, "<z8>%f</z8>\n", mitoz[mito_counter][7]);
	fprintf(file, "<mitocentrex>%f</mitocentrex>\n", mitocentrex[mito_counter]);
	fprintf(file, "<mitocentrey>%f</mitocentrey>\n", mitocentrey[mito_counter]);
	fprintf(file, "<mitocentrez>%f</mitocentrez>\n", mitocentrez[mito_counter]);
	fprintf(file, "<mitolength>%f</mitolength>\n", mitolength[mito_counter]);
	fprintf(file, "<mitoheight>%f</mitoheight>\n", mitoheight[mito_counter]);
	fprintf(file, "<mitowidth>%f</mitowidth>\n", mitowidth[mito_counter]);
	fprintf(file, "<radius>%f</radius>\n", mitoradius[mito_counter]);
	fprintf(file, "<radius2>%f</radius2>\n", mitoradius2[mito_counter]);
	fprintf(file, "<mito_collision_counter>%f</mito_collision_counter>\n", mito_collision_counter[mito_counter]);
	fprintf(file, "<mito_energetics>%f</mito_energetics>\n", mito_energetics[mito_counter]);
	fprintf(file, "<current_mitosize>%d</current_mitosize>\n", current_mitosize[mito_counter]);
	
	fprintf(file, "</xagent>\n");
	
	mito_counter++;
	Mito_Agents++;
}


// used variables are set free
	for (int i = 0; i < mito_number; i++) {
        free(mitox[i]);
    }
	for (int i = 0; i < mito_number; i++) {
			free(mitoy[i]);
		}
		free(mitoy);
	for (int i = 0; i < mito_number; i++) {
        free(mitoz[i]);
    }
    free(mitoz);	

free(mitoradius);
free(mitoradius2);
free(mitolength);
free(mitoheight);
free(mitowidth);
free(mitocentrex);
free(mitocentrey);
free(mitocentrez);
free(mito_collision_counter);
free(mito_energetics);
free(current_mitosize);


// creates random seed on mitochondria
double mitoSA = 4*(mitolength[0] * mitowidth[0]) + 2*(mitolength[0]*mitoheight[0]);
double prob_smallface = (mitolength[0]*mitoheight[0]) / mitoSA;
double prob_bigface = (mitowidth[0]*mitoheight[0]) / mitoSA;
double probability_boundaries[7];
probability_boundaries[0] = 0; // Min edge1
probability_boundaries[1] = probability_boundaries[0] + prob_bigface; // Max edge1, Min edge2
probability_boundaries[2] = probability_boundaries[1] + prob_bigface; // Max edge2, Min edge3
probability_boundaries[3] = probability_boundaries[2] + prob_bigface; // Max edge3, Min edge4
probability_boundaries[4] = probability_boundaries[3] + prob_bigface; // Max edge4, Min edge5
probability_boundaries[5] = probability_boundaries[4] + prob_smallface; // Max edge5, Min edge6
probability_boundaries[6] = probability_boundaries[5] + prob_smallface; // Max edge6: should be 1

/* Seed the random number generator */
srand (time (NULL));

/********************************************************************************************************/
/************************************************ Ligands ***********************************************/
/********************************************************************************************************/
/* Create MCL-1 Ligand Agents */
printf ("Creating MCL-1 Ligand Agents\n");
int state2_counter = 0;
double x, y, z;
while (state2_counter < number_state_2)
{
	double collision_counter_ligand = 0;
	int check_collision = 0;
	int foundedge = 0;
	int mitoid = 0;
	int mito_counter = 0;
	int intersection = 0;
	
	//Prevents intersections
	while (intersection == 0)
	{
		mito_counter = 0;
		x = minintracellular_xboundary + (rand()/(double)(RAND_MAX)*(maxintracellular_xboundary-minintracellular_xboundary));
		y = minintracellular_yboundary + (rand()/(double)(RAND_MAX)*(maxintracellular_yboundary-minintracellular_yboundary));
		z = minintracellular_zboundary + (rand()/(double)(RAND_MAX)*(maxintracellular_zboundary-minintracellular_zboundary));
		intersection = 1;
		// check-up that Ligand is not spawning in the blackbox/nucleus
		/*if (intersection == 1)
		{
			if ((x >= minnucleus_xboundary) && (x <= maxnucleus_xboundary) && ((y >= minnucleus_yboundary) && (y <= maxnucleus_yboundary)) && ((z > minnucleus_zboundary) && (z <= maxnucleus_zboundary)))
				{
					intersection = 0;
				}
		}*/
		// check-up that Ligand is not spawning in mitochondria
		while (intersection == 1 && (mito_counter < mito_number))
		{
			if (((x >= mitox[mito_counter][0]) && (x <= mitox[mito_counter][1])) && ((y >= mitoy[mito_counter][0]) && (y <= mitoy[mito_counter][2])) && ((z >= mitoz[mito_counter][0]) && (z <= mitoz[mito_counter][4])))
			{
				intersection = 0;
			}
			mito_counter++;
			if (intersection == 0) { mito_counter = mito_number; }
		}
	}

    fprintf(file, "<xagent>\n");
    fprintf(file, "<name>Ligand</name>\n");
    fprintf(file, "<id>%d</id>\n", Ligand_Agents);
    fprintf(file, "<state>2</state>\n");
    fprintf(file, "<Pending>-1</Pending>\n");
    fprintf(file, "<active>1</active>\n");
    fprintf(file, "<x>%f</x>\n", x);
    fprintf(file, "<y>%f</y>\n", y);
    fprintf(file, "<z>%f</z>\n", z);
    fprintf(file, "<Dc>%f</Dc>\n", Dc_state_2);
	fprintf(file, "<collision_counter_ligand>%f</collision_counter_ligand>\n", collision_counter_ligand);
	fprintf(file, "<check_collision>%d</check_collision>\n", check_collision);
	fprintf(file, "<foundedge>%d</foundedge>\n", foundedge);
	fprintf(file, "<mitoid>%d</mitoid>\n", mitoid);
    fprintf(file, "</xagent>\n");
state2_counter++;
Ligand_Agents++;
}

/* Create Mcl1Mito Receptor Agents */
printf ("Creating Mcl1Mito Receptor Agents\n");
int state20_counter = 0;
while (state20_counter < number_state_20)
{
	double collision_counter_ligand = 0;
	int check_collision = 0;
	int foundedge = 0;
	
	// creates random seed on mitochondria
	int mitoid = mito_number + 1;
	while (mitoid > mito_number || mitoid == mito_number)
	{
	mitoid = (int) rand()/(double)(RAND_MAX)*mito_number; // Rounds down a random double value between 0 and mito_number
	}
	
	double random_edge = rand()/(double)(RAND_MAX);
	// Edge 1: Formed by V1,V2,V5 and V6
	if (random_edge >= probability_boundaries[0] && random_edge < probability_boundaries[1])
	{
		foundedge = 1;
		x = mitox[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][1] - mitox[mitoid][0]));
		y = mitoy[mitoid][0];
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][4] - mitoz[mitoid][0]));
	}
	// Edge 2: Formed by V2,V4,V6 and V8
	if (random_edge >= probability_boundaries[1] && random_edge < probability_boundaries[2])
	{
		foundedge = 2;
		x = mitox[mitoid][1];
		y = mitoy[mitoid][1] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][3] - mitoy[mitoid][1]));
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][5] - mitoz[mitoid][1]));
	}
	// Edge 3: Formed by V3,V4,V7 and V8
	if (random_edge >= probability_boundaries[2] && random_edge < probability_boundaries[3])
	{
		foundedge = 3;
		x = mitox[mitoid][2] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][3] - mitox[mitoid][2]));
		y = mitoy[mitoid][2];
		z = mitoz[mitoid][2] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][6] - mitoz[mitoid][2]));
	}
	// Edge 4: Formed by V1,V3,V5 and V7
	if (random_edge >= probability_boundaries[3] && random_edge < probability_boundaries[4])
	{
		foundedge = 4;
		x = mitox[mitoid][0];
		y = mitoy[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][2] - mitoy[mitoid][0]));
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][4] - mitoz[mitoid][0]));
	}
	// Edge 5: Formed by V1,V2,V3 and V4
	if (random_edge >= probability_boundaries[4] && random_edge < probability_boundaries[5])
	{
		foundedge = 5;
		x = mitox[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][1] - mitox[mitoid][0]));
		y = mitoy[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][2] - mitoy[mitoid][0]));
		z = mitoz[mitoid][0];
	}
	// Edge 6: Formed by V5,V6,V7 and V8
	if (random_edge >= probability_boundaries[5] && random_edge < probability_boundaries[6])
	{
		foundedge = 6;
		x = mitox[mitoid][4] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][5] - mitox[mitoid][4]));
		y = mitoy[mitoid][4] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][6] - mitoy[mitoid][4]));
		z = mitoz[mitoid][4];
	}

	fprintf(file, "<xagent>\n");
    fprintf(file, "<name>Receptor</name>\n");
    fprintf(file, "<id>%d</id>\n", Receptor_Agents);
    fprintf(file, "<state>20</state>\n");
    fprintf(file, "<Pending>-1</Pending>\n");
    fprintf(file, "<active>1</active>\n");
    fprintf(file, "<x>%f</x>\n", x);
    fprintf(file, "<y>%f</y>\n", y);
    fprintf(file, "<z>%f</z>\n", z);
    fprintf(file, "<Dc>%f</Dc>\n", Dc_state_20);
	fprintf(file, "<collision_counter_ligand>%f</collision_counter_ligand>\n", collision_counter_ligand);
	fprintf(file, "<check_collision>%d</check_collision>\n", check_collision);
	fprintf(file, "<foundedge>%d</foundedge>\n", foundedge);
	fprintf(file, "<mitoid>%d</mitoid>\n", mitoid);
    fprintf(file, "</xagent>\n");


state20_counter++;
Receptor_Agents++;
}


/* Create Bclxl Ligand Agents */
printf ("Creating Bclxl Ligand Agents\n");
int state3_counter = 0;
while (state3_counter < number_state_3)
{
	double collision_counter_ligand = 0;
	int check_collision = 0;
	int foundedge = 0;
	int mitoid = 0;
	int mito_counter = 0;
	int intersection = 0;
	
	while (intersection == 0)
	{
		mito_counter = 0;
		x = minintracellular_xboundary + (rand()/(double)(RAND_MAX)*(maxintracellular_xboundary-minintracellular_xboundary));
		y = minintracellular_yboundary + (rand()/(double)(RAND_MAX)*(maxintracellular_yboundary-minintracellular_yboundary));
		z = minintracellular_zboundary + (rand()/(double)(RAND_MAX)*(maxintracellular_zboundary-minintracellular_zboundary));
		intersection = 1;
		while (mito_counter < mito_number)
		{
			if (((x > mitox[mito_counter][0]) && (x < mitox[mito_counter][1])) && ((y > mitoy[mito_counter][0]) && (y < mitoy[mito_counter][2])) && ((z > mitoz[mito_counter][0]) && (z < mitoz[mito_counter][4])))
			{
				intersection = 0;
			}
			mito_counter++;
			if (intersection == 0) { mito_counter = mito_number; }
		}
	}

    fprintf(file, "<xagent>\n");
    fprintf(file, "<name>Ligand</name>\n");
    fprintf(file, "<id>%d</id>\n", Ligand_Agents);
    fprintf(file, "<state>3</state>\n");
    fprintf(file, "<Pending>-1</Pending>\n");
    fprintf(file, "<active>1</active>\n");
    fprintf(file, "<x>%f</x>\n", x);
    fprintf(file, "<y>%f</y>\n", y);
    fprintf(file, "<z>%f</z>\n", z);
    fprintf(file, "<Dc>%f</Dc>\n", Dc_state_2);
	fprintf(file, "<collision_counter_ligand>%f</collision_counter_ligand>\n", collision_counter_ligand);
	fprintf(file, "<check_collision>%d</check_collision>\n", check_collision);
	fprintf(file, "<foundedge>%d</foundedge>\n", foundedge);
	fprintf(file, "<mitoid>%d</mitoid>\n", mitoid);
    fprintf(file, "</xagent>\n");
state3_counter++;
Ligand_Agents++;
}

/* Create Bax Ligand Agents */
printf ("Creating Bax Ligand Agents\n");
int state4_counter = 0;
while (state4_counter < number_state_4)
{
	double collision_counter_ligand = 0;
	int check_collision = 0;
	int foundedge = 0;
	int mitoid = 0;
	int mito_counter = 0;
	int intersection = 0;
	
	while (intersection == 0)
	{
		mito_counter = 0;
		x = minintracellular_xboundary + (rand()/(double)(RAND_MAX)*(maxintracellular_xboundary-minintracellular_xboundary));
		y = minintracellular_yboundary + (rand()/(double)(RAND_MAX)*(maxintracellular_yboundary-minintracellular_yboundary));
		z = minintracellular_zboundary + (rand()/(double)(RAND_MAX)*(maxintracellular_zboundary-minintracellular_zboundary));
		intersection = 1;
		while (mito_counter < mito_number)
		{
			if (((x > mitox[mito_counter][0]) && (x < mitox[mito_counter][1])) && ((y > mitoy[mito_counter][0]) && (y < mitoy[mito_counter][2])) && ((z > mitoz[mito_counter][0]) && (z < mitoz[mito_counter][4])))
			{
				intersection = 0;
			}
			mito_counter++;
			if (intersection == 0) { mito_counter = mito_number; }
		}
	}

    fprintf(file, "<xagent>\n");
    fprintf(file, "<name>Ligand</name>\n");
    fprintf(file, "<id>%d</id>\n", Ligand_Agents);
    fprintf(file, "<state>4</state>\n");
    fprintf(file, "<Pending>-1</Pending>\n");
    fprintf(file, "<active>1</active>\n");
    fprintf(file, "<x>%f</x>\n", x);
    fprintf(file, "<y>%f</y>\n", y);
    fprintf(file, "<z>%f</z>\n", z);
    fprintf(file, "<Dc>%f</Dc>\n", Dc_state_2);
	fprintf(file, "<collision_counter_ligand>%f</collision_counter_ligand>\n", collision_counter_ligand);
	fprintf(file, "<check_collision>%d</check_collision>\n", check_collision);
	fprintf(file, "<foundedge>%d</foundedge>\n", foundedge);
	fprintf(file, "<mitoid>%d</mitoid>\n", mitoid);
    fprintf(file, "</xagent>\n");
state4_counter++;
Ligand_Agents++;
}

/* Create aBax Receptor Agents */
printf ("Creating aBax Receptor Agents\n");
int state41_counter = 0;
while (state41_counter < number_state_41)
{
	double collision_counter_ligand = 0;
	int check_collision = 0;
	int foundedge = 0;

    // creates random seed on mitochondria
	int mitoid = mito_number + 1;
	while (mitoid > mito_number || mitoid == mito_number)
	{
	mitoid = (int) rand()/(double)(RAND_MAX)*mito_number; // Rounds down a random double value between 0 and mito_number
	}
	
	double random_edge = rand()/(double)(RAND_MAX);
	// Edge 1: Formed by V1,V2,V5 and V6
	if (random_edge >= probability_boundaries[0] && random_edge < probability_boundaries[1])
	{
		foundedge = 1;
		x = mitox[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][1] - mitox[mitoid][0]));
		y = mitoy[mitoid][0];
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][4] - mitoz[mitoid][0]));
	}
	// Edge 2: Formed by V2,V4,V6 and V8
	if (random_edge >= probability_boundaries[1] && random_edge < probability_boundaries[2])
	{
		foundedge = 2;
		x = mitox[mitoid][1];
		y = mitoy[mitoid][1] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][3] - mitoy[mitoid][1]));
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][5] - mitoz[mitoid][1]));
	}
	// Edge 3: Formed by V3,V4,V7 and V8
	if (random_edge >= probability_boundaries[2] && random_edge < probability_boundaries[3])
	{
		foundedge = 3;
		x = mitox[mitoid][2] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][3] - mitox[mitoid][2]));
		y = mitoy[mitoid][2];
		z = mitoz[mitoid][2] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][6] - mitoz[mitoid][2]));
	}
	// Edge 4: Formed by V1,V3,V5 and V7
	if (random_edge >= probability_boundaries[3] && random_edge < probability_boundaries[4])
	{
		foundedge = 4;
		x = mitox[mitoid][0];
		y = mitoy[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][2] - mitoy[mitoid][0]));
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][4] - mitoz[mitoid][0]));
	}
	// Edge 5: Formed by V1,V2,V3 and V4
	if (random_edge >= probability_boundaries[4] && random_edge < probability_boundaries[5])
	{
		foundedge = 5;
		x = mitox[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][1] - mitox[mitoid][0]));
		y = mitoy[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][2] - mitoy[mitoid][0]));
		z = mitoz[mitoid][0];
	}
	// Edge 6: Formed by V5,V6,V7 and V8
	if (random_edge >= probability_boundaries[5] && random_edge < probability_boundaries[6])
	{
		foundedge = 6;
		x = mitox[mitoid][4] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][5] - mitox[mitoid][4]));
		y = mitoy[mitoid][4] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][6] - mitoy[mitoid][4]));
		z = mitoz[mitoid][4];
	}

    fprintf(file, "<xagent>\n");
    fprintf(file, "<name>Receptor</name>\n");
    fprintf(file, "<id>%d</id>\n", Receptor_Agents);
    fprintf(file, "<state>41</state>\n");
    fprintf(file, "<Pending>-1</Pending>\n");
    fprintf(file, "<active>1</active>\n");
    fprintf(file, "<x>%f</x>\n", x);
    fprintf(file, "<y>%f</y>\n", y);
    fprintf(file, "<z>%f</z>\n", z);
    fprintf(file, "<Dc>%f</Dc>\n", Dc_state_41);
	fprintf(file, "<collision_counter_ligand>%f</collision_counter_ligand>\n", collision_counter_ligand);
	fprintf(file, "<check_collision>%d</check_collision>\n", check_collision);
	fprintf(file, "<foundedge>%d</foundedge>\n", foundedge);
	fprintf(file, "<mitoid>%d</mitoid>\n", mitoid);
    fprintf(file, "</xagent>\n");
state41_counter++;
Receptor_Agents++;
}

/* Create aBax2 Receptor Agents */
printf ("Creating aBax2 Receptor Agents\n");
int state42_counter = 0;
while (state42_counter < number_state_42)
{
	double collision_counter_ligand = 0;
	int check_collision = 0;
	int foundedge = 0;
	
	// creates random seed on mitochondria
	int mitoid = mito_number + 1;
	while (mitoid > mito_number || mitoid == mito_number)
	{
	mitoid = (int) rand()/(double)(RAND_MAX)*mito_number; // Rounds down a random double value between 0 and mito_number
	}
	
	double random_edge = rand()/(double)(RAND_MAX);
	// Edge 1: Formed by V1,V2,V5 and V6
	if (random_edge >= probability_boundaries[0] && random_edge < probability_boundaries[1])
	{
		foundedge = 1;
		x = mitox[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][1] - mitox[mitoid][0]));
		y = mitoy[mitoid][0];
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][4] - mitoz[mitoid][0]));
	}
	// Edge 2: Formed by V2,V4,V6 and V8
	if (random_edge >= probability_boundaries[1] && random_edge < probability_boundaries[2])
	{
		foundedge = 2;
		x = mitox[mitoid][1];
		y = mitoy[mitoid][1] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][3] - mitoy[mitoid][1]));
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][5] - mitoz[mitoid][1]));
	}
	// Edge 3: Formed by V3,V4,V7 and V8
	if (random_edge >= probability_boundaries[2] && random_edge < probability_boundaries[3])
	{
		foundedge = 3;
		x = mitox[mitoid][2] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][3] - mitox[mitoid][2]));
		y = mitoy[mitoid][2];
		z = mitoz[mitoid][2] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][6] - mitoz[mitoid][2]));
	}
	// Edge 4: Formed by V1,V3,V5 and V7
	if (random_edge >= probability_boundaries[3] && random_edge < probability_boundaries[4])
	{
		foundedge = 4;
		x = mitox[mitoid][0];
		y = mitoy[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][2] - mitoy[mitoid][0]));
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][4] - mitoz[mitoid][0]));
	}
	// Edge 5: Formed by V1,V2,V3 and V4
	if (random_edge >= probability_boundaries[4] && random_edge < probability_boundaries[5])
	{
		foundedge = 5;
		x = mitox[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][1] - mitox[mitoid][0]));
		y = mitoy[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][2] - mitoy[mitoid][0]));
		z = mitoz[mitoid][0];
	}
	// Edge 6: Formed by V5,V6,V7 and V8
	if (random_edge >= probability_boundaries[5] && random_edge < probability_boundaries[6])
	{
		foundedge = 6;
		x = mitox[mitoid][4] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][5] - mitox[mitoid][4]));
		y = mitoy[mitoid][4] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][6] - mitoy[mitoid][4]));
		z = mitoz[mitoid][4];
	}

    fprintf(file, "<xagent>\n");
    fprintf(file, "<name>Receptor</name>\n");
    fprintf(file, "<id>%d</id>\n", Receptor_Agents);
    fprintf(file, "<state>42</state>\n");
    fprintf(file, "<Pending>-1</Pending>\n");
    fprintf(file, "<active>1</active>\n");
    fprintf(file, "<x>%f</x>\n", x);
    fprintf(file, "<y>%f</y>\n", y);
    fprintf(file, "<z>%f</z>\n", z);
    fprintf(file, "<Dc>%f</Dc>\n", Dc_state_42);
	fprintf(file, "<collision_counter_ligand>%f</collision_counter_ligand>\n", collision_counter_ligand);
	fprintf(file, "<check_collision>%d</check_collision>\n", check_collision);
	fprintf(file, "<foundedge>%d</foundedge>\n", foundedge);
	fprintf(file, "<mitoid>%d</mitoid>\n", mitoid);
    fprintf(file, "</xagent>\n");
state42_counter++;
Receptor_Agents++;
}

/* Create aBax4 Receptor Agents */
printf ("Creating aBax4 Receptor Agents\n");
int state44_counter = 0;
while (state44_counter < number_state_44)
{
	double collision_counter_ligand = 0;
	int check_collision = 0;
	int foundedge = 0;
	
	// creates random seed on mitochondria
	int mitoid = mito_number + 1;
	while (mitoid > mito_number || mitoid == mito_number)
	{
	mitoid = (int) rand()/(double)(RAND_MAX)*mito_number; // Rounds down a random double value between 0 and mito_number
	}
	
	double random_edge = rand()/(double)(RAND_MAX);
	// Edge 1: Formed by V1,V2,V5 and V6
	if (random_edge >= probability_boundaries[0] && random_edge < probability_boundaries[1])
	{
		foundedge = 1;
		x = mitox[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][1] - mitox[mitoid][0]));
		y = mitoy[mitoid][0];
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][4] - mitoz[mitoid][0]));
	}
	// Edge 2: Formed by V2,V4,V6 and V8
	if (random_edge >= probability_boundaries[1] && random_edge < probability_boundaries[2])
	{
		foundedge = 2;
		x = mitox[mitoid][1];
		y = mitoy[mitoid][1] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][3] - mitoy[mitoid][1]));
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][5] - mitoz[mitoid][1]));
	}
	// Edge 3: Formed by V3,V4,V7 and V8
	if (random_edge >= probability_boundaries[2] && random_edge < probability_boundaries[3])
	{
		foundedge = 3;
		x = mitox[mitoid][2] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][3] - mitox[mitoid][2]));
		y = mitoy[mitoid][2];
		z = mitoz[mitoid][2] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][6] - mitoz[mitoid][2]));
	}
	// Edge 4: Formed by V1,V3,V5 and V7
	if (random_edge >= probability_boundaries[3] && random_edge < probability_boundaries[4])
	{
		foundedge = 4;
		x = mitox[mitoid][0];
		y = mitoy[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][2] - mitoy[mitoid][0]));
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][4] - mitoz[mitoid][0]));
	}
	// Edge 5: Formed by V1,V2,V3 and V4
	if (random_edge >= probability_boundaries[4] && random_edge < probability_boundaries[5])
	{
		foundedge = 5;
		x = mitox[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][1] - mitox[mitoid][0]));
		y = mitoy[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][2] - mitoy[mitoid][0]));
		z = mitoz[mitoid][0];
	}
	// Edge 6: Formed by V5,V6,V7 and V8
	if (random_edge >= probability_boundaries[5] && random_edge < probability_boundaries[6])
	{
		foundedge = 6;
		x = mitox[mitoid][4] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][5] - mitox[mitoid][4]));
		y = mitoy[mitoid][4] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][6] - mitoy[mitoid][4]));
		z = mitoz[mitoid][4];
	}

    fprintf(file, "<xagent>\n");
    fprintf(file, "<name>Receptor</name>\n");
    fprintf(file, "<id>%d</id>\n", Receptor_Agents);
    fprintf(file, "<state>44</state>\n");
    fprintf(file, "<Pending>-1</Pending>\n");
    fprintf(file, "<active>1</active>\n");
    fprintf(file, "<x>%f</x>\n", x);
    fprintf(file, "<y>%f</y>\n", y);
    fprintf(file, "<z>%f</z>\n", z);
    fprintf(file, "<Dc>%f</Dc>\n", Dc_state_44);
	fprintf(file, "<collision_counter_ligand>%f</collision_counter_ligand>\n", collision_counter_ligand);
	fprintf(file, "<check_collision>%d</check_collision>\n", check_collision);
	fprintf(file, "<foundedge>%d</foundedge>\n", foundedge);
	fprintf(file, "<mitoid>%d</mitoid>\n", mitoid);
    fprintf(file, "</xagent>\n");
state44_counter++;
Receptor_Agents++;
}

/* Create aBax6 Receptor Agents */
printf ("Creating aBax6 Receptor Agents\n");
int state46_counter = 0;
while (state46_counter < number_state_46)
{
	double collision_counter_ligand = 0;
	int check_collision = 0;
	int foundedge = 0;
	
    // creates random seed on mitochondria
	int mitoid = mito_number + 1;
	while (mitoid > mito_number || mitoid == mito_number)
	{
	mitoid = (int) rand()/(double)(RAND_MAX)*mito_number; // Rounds down a random double value between 0 and mito_number
	}
	
	double random_edge = rand()/(double)(RAND_MAX);
	// Edge 1: Formed by V1,V2,V5 and V6
	if (random_edge >= probability_boundaries[0] && random_edge < probability_boundaries[1])
	{
		foundedge = 1;
		x = mitox[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][1] - mitox[mitoid][0]));
		y = mitoy[mitoid][0];
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][4] - mitoz[mitoid][0]));
	}
	// Edge 2: Formed by V2,V4,V6 and V8
	if (random_edge >= probability_boundaries[1] && random_edge < probability_boundaries[2])
	{
		foundedge = 2;
		x = mitox[mitoid][1];
		y = mitoy[mitoid][1] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][3] - mitoy[mitoid][1]));
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][5] - mitoz[mitoid][1]));
	}
	// Edge 3: Formed by V3,V4,V7 and V8
	if (random_edge >= probability_boundaries[2] && random_edge < probability_boundaries[3])
	{
		foundedge = 3;
		x = mitox[mitoid][2] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][3] - mitox[mitoid][2]));
		y = mitoy[mitoid][2];
		z = mitoz[mitoid][2] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][6] - mitoz[mitoid][2]));
	}
	// Edge 4: Formed by V1,V3,V5 and V7
	if (random_edge >= probability_boundaries[3] && random_edge < probability_boundaries[4])
	{
		foundedge = 4;
		x = mitox[mitoid][0];
		y = mitoy[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][2] - mitoy[mitoid][0]));
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][4] - mitoz[mitoid][0]));
	}
	// Edge 5: Formed by V1,V2,V3 and V4
	if (random_edge >= probability_boundaries[4] && random_edge < probability_boundaries[5])
	{
		foundedge = 5;
		x = mitox[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][1] - mitox[mitoid][0]));
		y = mitoy[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][2] - mitoy[mitoid][0]));
		z = mitoz[mitoid][0];
	}
	// Edge 6: Formed by V5,V6,V7 and V8
	if (random_edge >= probability_boundaries[5] && random_edge < probability_boundaries[6])
	{
		foundedge = 6;
		x = mitox[mitoid][4] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][5] - mitox[mitoid][4]));
		y = mitoy[mitoid][4] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][6] - mitoy[mitoid][4]));
		z = mitoz[mitoid][4];
	}

    fprintf(file, "<xagent>\n");
    fprintf(file, "<name>Receptor</name>\n");
    fprintf(file, "<id>%d</id>\n", Receptor_Agents);
    fprintf(file, "<state>46</state>\n");
    fprintf(file, "<Pending>-1</Pending>\n");
    fprintf(file, "<active>1</active>\n");
    fprintf(file, "<x>%f</x>\n", x);
    fprintf(file, "<y>%f</y>\n", y);
    fprintf(file, "<z>%f</z>\n", z);
    fprintf(file, "<Dc>%f</Dc>\n", Dc_state_46);
	fprintf(file, "<collision_counter_ligand>%f</collision_counter_ligand>\n", collision_counter_ligand);
	fprintf(file, "<check_collision>%d</check_collision>\n", check_collision);
	fprintf(file, "<foundedge>%d</foundedge>\n", foundedge);
	fprintf(file, "<mitoid>%d</mitoid>\n", mitoid);
    fprintf(file, "</xagent>\n");
state46_counter++;
Receptor_Agents++;
}


/* Create tBid Ligand Agents */
printf ("Creating tBid Ligand Agents\n");
int state5_counter = 0;
while (state5_counter < number_state_5)
{
	double collision_counter_ligand = 0;
	int check_collision = 0;
	int foundedge = 0;
	int mitoid = 0;
	int mito_counter = 0;
	int intersection = 0;
	
	while (intersection == 0)
	{
		mito_counter = 0;
		x = minintracellular_xboundary + (rand()/(double)(RAND_MAX)*(maxintracellular_xboundary-minintracellular_xboundary));
		y = minintracellular_yboundary + (rand()/(double)(RAND_MAX)*(maxintracellular_yboundary-minintracellular_yboundary));
		z = minintracellular_zboundary + (rand()/(double)(RAND_MAX)*(maxintracellular_zboundary-minintracellular_zboundary));
		intersection = 1;
		while (mito_counter < mito_number)
		{
			if (((x > mitox[mito_counter][0]) && (x < mitox[mito_counter][1])) && ((y > mitoy[mito_counter][0]) && (y < mitoy[mito_counter][2])) && ((z > mitoz[mito_counter][0]) && (z < mitoz[mito_counter][4])))
			{
				intersection = 0;
			}
			mito_counter++;
			if (intersection == 0) { mito_counter = mito_number; }
		}
	}

    fprintf(file, "<xagent>\n");
    fprintf(file, "<name>Ligand</name>\n");
    fprintf(file, "<id>%d</id>\n", Ligand_Agents);
    fprintf(file, "<state>5</state>\n");
    fprintf(file, "<Pending>-1</Pending>\n");
    fprintf(file, "<active>1</active>\n");
    fprintf(file, "<x>%f</x>\n", x);
    fprintf(file, "<y>%f</y>\n", y);
    fprintf(file, "<z>%f</z>\n", z);
    fprintf(file, "<Dc>%f</Dc>\n", Dc_state_5);
	fprintf(file, "<collision_counter_ligand>%f</collision_counter_ligand>\n", collision_counter_ligand);
	fprintf(file, "<check_collision>%d</check_collision>\n", check_collision);
	fprintf(file, "<foundedge>%d</foundedge>\n", foundedge);
	fprintf(file, "<mitoid>%d</mitoid>\n", mitoid);
    fprintf(file, "</xagent>\n");
state5_counter++;
Ligand_Agents++;
}

/* Create Mcl1Mito Receptor Agents */
printf ("Creating Mcl1Mito Receptor Agents\n");
int state56_counter = 0;
while (state56_counter < number_state_56)
{
	double collision_counter_ligand = 0;
	int check_collision = 0;
	int foundedge = 0;
	
	// creates random seed on mitochondria
	int mitoid = mito_number + 1;
	while (mitoid > mito_number || mitoid == mito_number)
	{
	mitoid = (int) rand()/(double)(RAND_MAX)*mito_number; // Rounds down a random double value between 0 and mito_number
	}
	
	double random_edge = rand()/(double)(RAND_MAX);
	// Edge 1: Formed by V1,V2,V5 and V6
	if (random_edge >= probability_boundaries[0] && random_edge < probability_boundaries[1])
	{
		foundedge = 1;
		x = mitox[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][1] - mitox[mitoid][0]));
		y = mitoy[mitoid][0];
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][4] - mitoz[mitoid][0]));
	}
	// Edge 2: Formed by V2,V4,V6 and V8
	if (random_edge >= probability_boundaries[1] && random_edge < probability_boundaries[2])
	{
		foundedge = 2;
		x = mitox[mitoid][1];
		y = mitoy[mitoid][1] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][3] - mitoy[mitoid][1]));
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][5] - mitoz[mitoid][1]));
	}
	// Edge 3: Formed by V3,V4,V7 and V8
	if (random_edge >= probability_boundaries[2] && random_edge < probability_boundaries[3])
	{
		foundedge = 3;
		x = mitox[mitoid][2] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][3] - mitox[mitoid][2]));
		y = mitoy[mitoid][2];
		z = mitoz[mitoid][2] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][6] - mitoz[mitoid][2]));
	}
	// Edge 4: Formed by V1,V3,V5 and V7
	if (random_edge >= probability_boundaries[3] && random_edge < probability_boundaries[4])
	{
		foundedge = 4;
		x = mitox[mitoid][0];
		y = mitoy[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][2] - mitoy[mitoid][0]));
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][4] - mitoz[mitoid][0]));
	}
	// Edge 5: Formed by V1,V2,V3 and V4
	if (random_edge >= probability_boundaries[4] && random_edge < probability_boundaries[5])
	{
		foundedge = 5;
		x = mitox[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][1] - mitox[mitoid][0]));
		y = mitoy[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][2] - mitoy[mitoid][0]));
		z = mitoz[mitoid][0];
	}
	// Edge 6: Formed by V5,V6,V7 and V8
	if (random_edge >= probability_boundaries[5] && random_edge < probability_boundaries[6])
	{
		foundedge = 6;
		x = mitox[mitoid][4] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][5] - mitox[mitoid][4]));
		y = mitoy[mitoid][4] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][6] - mitoy[mitoid][4]));
		z = mitoz[mitoid][4];
	}

	fprintf(file, "<xagent>\n");
    fprintf(file, "<name>Receptor</name>\n");
    fprintf(file, "<id>%d</id>\n", Receptor_Agents);
    fprintf(file, "<state>56</state>\n");
    fprintf(file, "<Pending>-1</Pending>\n");
    fprintf(file, "<active>1</active>\n");
    fprintf(file, "<x>%f</x>\n", x);
    fprintf(file, "<y>%f</y>\n", y);
    fprintf(file, "<z>%f</z>\n", z);
    fprintf(file, "<Dc>%f</Dc>\n", Dc_state_56);
	fprintf(file, "<collision_counter_ligand>%f</collision_counter_ligand>\n", collision_counter_ligand);
	fprintf(file, "<check_collision>%d</check_collision>\n", check_collision);
	fprintf(file, "<foundedge>%d</foundedge>\n", foundedge);
	fprintf(file, "<mitoid>%d</mitoid>\n", mitoid);
    fprintf(file, "</xagent>\n");
state56_counter++;
Receptor_Agents++;
}


/* Create tBidBax Receptor Agents */
printf ("Creating tBidBax Receptor Agents\n");
int state51_counter = 0;
while (state51_counter < number_state_51)
{
	double collision_counter_ligand = 0;
	int check_collision = 0;
	int foundedge = 0;
	
    // creates random seed on mitochondria
	int mitoid = mito_number + 1;
	while (mitoid > mito_number || mitoid == mito_number)
	{
	mitoid = (int) rand()/(double)(RAND_MAX)*mito_number; // Rounds down a random double value between 0 and mito_number
	}
	
	double random_edge = rand()/(double)(RAND_MAX);
	// Edge 1: Formed by V1,V2,V5 and V6
	if (random_edge >= probability_boundaries[0] && random_edge < probability_boundaries[1])
	{
		foundedge = 1;
		x = mitox[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][1] - mitox[mitoid][0]));
		y = mitoy[mitoid][0];
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][4] - mitoz[mitoid][0]));
	}
	// Edge 2: Formed by V2,V4,V6 and V8
	if (random_edge >= probability_boundaries[1] && random_edge < probability_boundaries[2])
	{
		foundedge = 2;
		x = mitox[mitoid][1];
		y = mitoy[mitoid][1] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][3] - mitoy[mitoid][1]));
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][5] - mitoz[mitoid][1]));
	}
	// Edge 3: Formed by V3,V4,V7 and V8
	if (random_edge >= probability_boundaries[2] && random_edge < probability_boundaries[3])
	{
		foundedge = 3;
		x = mitox[mitoid][2] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][3] - mitox[mitoid][2]));
		y = mitoy[mitoid][2];
		z = mitoz[mitoid][2] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][6] - mitoz[mitoid][2]));
	}
	// Edge 4: Formed by V1,V3,V5 and V7
	if (random_edge >= probability_boundaries[3] && random_edge < probability_boundaries[4])
	{
		foundedge = 4;
		x = mitox[mitoid][0];
		y = mitoy[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][2] - mitoy[mitoid][0]));
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][4] - mitoz[mitoid][0]));
	}
	// Edge 5: Formed by V1,V2,V3 and V4
	if (random_edge >= probability_boundaries[4] && random_edge < probability_boundaries[5])
	{
		foundedge = 5;
		x = mitox[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][1] - mitox[mitoid][0]));
		y = mitoy[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][2] - mitoy[mitoid][0]));
		z = mitoz[mitoid][0];
	}
	// Edge 6: Formed by V5,V6,V7 and V8
	if (random_edge >= probability_boundaries[5] && random_edge < probability_boundaries[6])
	{
		foundedge = 6;
		x = mitox[mitoid][4] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][5] - mitox[mitoid][4]));
		y = mitoy[mitoid][4] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][6] - mitoy[mitoid][4]));
		z = mitoz[mitoid][4];
	}

    fprintf(file, "<xagent>\n");
    fprintf(file, "<name>Receptor</name>\n");
    fprintf(file, "<id>%d</id>\n", Receptor_Agents);
    fprintf(file, "<state>51</state>\n");
    fprintf(file, "<Pending>-1</Pending>\n");
    fprintf(file, "<active>1</active>\n");
    fprintf(file, "<x>%f</x>\n", x);
    fprintf(file, "<y>%f</y>\n", y);
    fprintf(file, "<z>%f</z>\n", z);
    fprintf(file, "<Dc>%f</Dc>\n", Dc_state_51);
	fprintf(file, "<collision_counter_ligand>%f</collision_counter_ligand>\n", collision_counter_ligand);
	fprintf(file, "<check_collision>%d</check_collision>\n", check_collision);
	fprintf(file, "<foundedge>%d</foundedge>\n", foundedge);
	fprintf(file, "<mitoid>%d</mitoid>\n", mitoid);
    fprintf(file, "</xagent>\n");
state51_counter++;
Receptor_Agents++;
}

/* Create tBidaBax Receptor Agents */
printf ("Creating tBidaBax Receptor Agents\n");
int state52_counter = 0;
while (state52_counter < number_state_52)
{
	double collision_counter_ligand = 0;
	int check_collision = 0;
	int foundedge = 0;
	
    // creates random seed on mitochondria
	int mitoid = mito_number + 1;
	while (mitoid > mito_number || mitoid == mito_number)
	{
	mitoid = (int) rand()/(double)(RAND_MAX)*mito_number; // Rounds down a random double value between 0 and mito_number
	}
	
	double random_edge = rand()/(double)(RAND_MAX);
	// Edge 1: Formed by V1,V2,V5 and V6
	if (random_edge >= probability_boundaries[0] && random_edge < probability_boundaries[1])
	{
		foundedge = 1;
		x = mitox[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][1] - mitox[mitoid][0]));
		y = mitoy[mitoid][0];
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][4] - mitoz[mitoid][0]));
	}
	// Edge 2: Formed by V2,V4,V6 and V8
	if (random_edge >= probability_boundaries[1] && random_edge < probability_boundaries[2])
	{
		foundedge = 2;
		x = mitox[mitoid][1];
		y = mitoy[mitoid][1] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][3] - mitoy[mitoid][1]));
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][5] - mitoz[mitoid][1]));
	}
	// Edge 3: Formed by V3,V4,V7 and V8
	if (random_edge >= probability_boundaries[2] && random_edge < probability_boundaries[3])
	{
		foundedge = 3;
		x = mitox[mitoid][2] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][3] - mitox[mitoid][2]));
		y = mitoy[mitoid][2];
		z = mitoz[mitoid][2] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][6] - mitoz[mitoid][2]));
	}
	// Edge 4: Formed by V1,V3,V5 and V7
	if (random_edge >= probability_boundaries[3] && random_edge < probability_boundaries[4])
	{
		foundedge = 4;
		x = mitox[mitoid][0];
		y = mitoy[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][2] - mitoy[mitoid][0]));
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][4] - mitoz[mitoid][0]));
	}
	// Edge 5: Formed by V1,V2,V3 and V4
	if (random_edge >= probability_boundaries[4] && random_edge < probability_boundaries[5])
	{
		foundedge = 5;
		x = mitox[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][1] - mitox[mitoid][0]));
		y = mitoy[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][2] - mitoy[mitoid][0]));
		z = mitoz[mitoid][0];
	}
	// Edge 6: Formed by V5,V6,V7 and V8
	if (random_edge >= probability_boundaries[5] && random_edge < probability_boundaries[6])
	{
		foundedge = 6;
		x = mitox[mitoid][4] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][5] - mitox[mitoid][4]));
		y = mitoy[mitoid][4] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][6] - mitoy[mitoid][4]));
		z = mitoz[mitoid][4];
	}

    fprintf(file, "<xagent>\n");
    fprintf(file, "<name>Receptor</name>\n");
    fprintf(file, "<id>%d</id>\n", Receptor_Agents);
    fprintf(file, "<state>52</state>\n");
    fprintf(file, "<Pending>-1</Pending>\n");
    fprintf(file, "<active>1</active>\n");
    fprintf(file, "<x>%f</x>\n", x);
    fprintf(file, "<y>%f</y>\n", y);
    fprintf(file, "<z>%f</z>\n", z);
    fprintf(file, "<Dc>%f</Dc>\n", Dc_state_52);
	fprintf(file, "<collision_counter_ligand>%f</collision_counter_ligand>\n", collision_counter_ligand);
	fprintf(file, "<check_collision>%d</check_collision>\n", check_collision);
	fprintf(file, "<foundedge>%d</foundedge>\n", foundedge);
	fprintf(file, "<mitoid>%d</mitoid>\n", mitoid);
    fprintf(file, "</xagent>\n");
state52_counter++;
Receptor_Agents++;
}

/* Create tBidBclxl Receptor Agents */
printf ("Creating tBidBclxl Receptor Agents\n");
int state53_counter = 0;
while (state53_counter < number_state_53)
{
	double collision_counter_ligand = 0;
	int check_collision = 0;
	int foundedge = 0;
	
    // creates random seed on mitochondria
	int mitoid = mito_number + 1;
	while (mitoid > mito_number || mitoid == mito_number)
	{
	mitoid = (int) rand()/(double)(RAND_MAX)*mito_number; // Rounds down a random double value between 0 and mito_number
	}
	
	double random_edge = rand()/(double)(RAND_MAX);
	// Edge 1: Formed by V1,V2,V5 and V6
	if (random_edge >= probability_boundaries[0] && random_edge < probability_boundaries[1])
	{
		foundedge = 1;
		x = mitox[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][1] - mitox[mitoid][0]));
		y = mitoy[mitoid][0];
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][4] - mitoz[mitoid][0]));
	}
	// Edge 2: Formed by V2,V4,V6 and V8
	if (random_edge >= probability_boundaries[1] && random_edge < probability_boundaries[2])
	{
		foundedge = 2;
		x = mitox[mitoid][1];
		y = mitoy[mitoid][1] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][3] - mitoy[mitoid][1]));
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][5] - mitoz[mitoid][1]));
	}
	// Edge 3: Formed by V3,V4,V7 and V8
	if (random_edge >= probability_boundaries[2] && random_edge < probability_boundaries[3])
	{
		foundedge = 3;
		x = mitox[mitoid][2] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][3] - mitox[mitoid][2]));
		y = mitoy[mitoid][2];
		z = mitoz[mitoid][2] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][6] - mitoz[mitoid][2]));
	}
	// Edge 4: Formed by V1,V3,V5 and V7
	if (random_edge >= probability_boundaries[3] && random_edge < probability_boundaries[4])
	{
		foundedge = 4;
		x = mitox[mitoid][0];
		y = mitoy[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][2] - mitoy[mitoid][0]));
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][4] - mitoz[mitoid][0]));
	}
	// Edge 5: Formed by V1,V2,V3 and V4
	if (random_edge >= probability_boundaries[4] && random_edge < probability_boundaries[5])
	{
		foundedge = 5;
		x = mitox[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][1] - mitox[mitoid][0]));
		y = mitoy[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][2] - mitoy[mitoid][0]));
		z = mitoz[mitoid][0];
	}
	// Edge 6: Formed by V5,V6,V7 and V8
	if (random_edge >= probability_boundaries[5] && random_edge < probability_boundaries[6])
	{
		foundedge = 6;
		x = mitox[mitoid][4] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][5] - mitox[mitoid][4]));
		y = mitoy[mitoid][4] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][6] - mitoy[mitoid][4]));
		z = mitoz[mitoid][4];
	}

    fprintf(file, "<xagent>\n");
    fprintf(file, "<name>Receptor</name>\n");
    fprintf(file, "<id>%d</id>\n", Receptor_Agents);
    fprintf(file, "<state>53</state>\n");
    fprintf(file, "<Pending>-1</Pending>\n");
    fprintf(file, "<active>1</active>\n");
    fprintf(file, "<x>%f</x>\n", x);
    fprintf(file, "<y>%f</y>\n", y);
    fprintf(file, "<z>%f</z>\n", z);
    fprintf(file, "<Dc>%f</Dc>\n", Dc_state_53);
	fprintf(file, "<collision_counter_ligand>%f</collision_counter_ligand>\n", collision_counter_ligand);
	fprintf(file, "<check_collision>%d</check_collision>\n", check_collision);
	fprintf(file, "<foundedge>%d</foundedge>\n", foundedge);
	fprintf(file, "<mitoid>%d</mitoid>\n", mitoid);
    fprintf(file, "</xagent>\n");
state53_counter++;
Receptor_Agents++;
}

/* Create tBidMcl1 Receptor Agents */
printf ("Creating tBidMcl1 Receptor Agents\n");
int state54_counter = 0;
while (state54_counter < number_state_54)
{
	double collision_counter_ligand = 0;
	int check_collision = 0;
	int foundedge = 0;
	
    // creates random seed on mitochondria
	int mitoid = mito_number + 1;
	while (mitoid > mito_number || mitoid == mito_number)
	{
	mitoid = (int) rand()/(double)(RAND_MAX)*mito_number; // Rounds down a random double value between 0 and mito_number
	}
	
	double random_edge = rand()/(double)(RAND_MAX);
	// Edge 1: Formed by V1,V2,V5 and V6
	if (random_edge >= probability_boundaries[0] && random_edge < probability_boundaries[1])
	{
		foundedge = 1;
		x = mitox[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][1] - mitox[mitoid][0]));
		y = mitoy[mitoid][0];
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][4] - mitoz[mitoid][0]));
	}
	// Edge 2: Formed by V2,V4,V6 and V8
	if (random_edge >= probability_boundaries[1] && random_edge < probability_boundaries[2])
	{
		foundedge = 2;
		x = mitox[mitoid][1];
		y = mitoy[mitoid][1] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][3] - mitoy[mitoid][1]));
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][5] - mitoz[mitoid][1]));
	}
	// Edge 3: Formed by V3,V4,V7 and V8
	if (random_edge >= probability_boundaries[2] && random_edge < probability_boundaries[3])
	{
		foundedge = 3;
		x = mitox[mitoid][2] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][3] - mitox[mitoid][2]));
		y = mitoy[mitoid][2];
		z = mitoz[mitoid][2] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][6] - mitoz[mitoid][2]));
	}
	// Edge 4: Formed by V1,V3,V5 and V7
	if (random_edge >= probability_boundaries[3] && random_edge < probability_boundaries[4])
	{
		foundedge = 4;
		x = mitox[mitoid][0];
		y = mitoy[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][2] - mitoy[mitoid][0]));
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][4] - mitoz[mitoid][0]));
	}
	// Edge 5: Formed by V1,V2,V3 and V4
	if (random_edge >= probability_boundaries[4] && random_edge < probability_boundaries[5])
	{
		foundedge = 5;
		x = mitox[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][1] - mitox[mitoid][0]));
		y = mitoy[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][2] - mitoy[mitoid][0]));
		z = mitoz[mitoid][0];
	}
	// Edge 6: Formed by V5,V6,V7 and V8
	if (random_edge >= probability_boundaries[5] && random_edge < probability_boundaries[6])
	{
		foundedge = 6;
		x = mitox[mitoid][4] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][5] - mitox[mitoid][4]));
		y = mitoy[mitoid][4] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][6] - mitoy[mitoid][4]));
		z = mitoz[mitoid][4];
	}

    fprintf(file, "<xagent>\n");
    fprintf(file, "<name>Receptor</name>\n");
    fprintf(file, "<id>%d</id>\n", Receptor_Agents);
    fprintf(file, "<state>54</state>\n");
    fprintf(file, "<Pending>-1</Pending>\n");
    fprintf(file, "<active>1</active>\n");
    fprintf(file, "<x>%f</x>\n", x);
    fprintf(file, "<y>%f</y>\n", y);
    fprintf(file, "<z>%f</z>\n", z);
    fprintf(file, "<Dc>%f</Dc>\n", Dc_state_54);
	fprintf(file, "<collision_counter_ligand>%f</collision_counter_ligand>\n", collision_counter_ligand);
	fprintf(file, "<check_collision>%d</check_collision>\n", check_collision);
	fprintf(file, "<foundedge>%d</foundedge>\n", foundedge);
	fprintf(file, "<mitoid>%d</mitoid>\n", mitoid);
    fprintf(file, "</xagent>\n");
state54_counter++;
Receptor_Agents++;
}

/* Create BclxlaBax Receptor Agents */
printf ("Creating BclxlaBax Receptor Agents\n");
int state55_counter = 0;
while (state55_counter < number_state_55)
{
	double collision_counter_ligand = 0;
	int check_collision = 0;
	int foundedge = 0;
	
	// creates random seed on mitochondria
	int mitoid = mito_number + 1;
	while (mitoid > mito_number || mitoid == mito_number)
	{
	mitoid = (int) rand()/(double)(RAND_MAX)*mito_number; // Rounds down a random double value between 0 and mito_number
	}
	
	double random_edge = rand()/(double)(RAND_MAX);
	// Edge 1: Formed by V1,V2,V5 and V6
	if (random_edge >= probability_boundaries[0] && random_edge < probability_boundaries[1])
	{
		foundedge = 1;
		x = mitox[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][1] - mitox[mitoid][0]));
		y = mitoy[mitoid][0];
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][4] - mitoz[mitoid][0]));
	}
	// Edge 2: Formed by V2,V4,V6 and V8
	if (random_edge >= probability_boundaries[1] && random_edge < probability_boundaries[2])
	{
		foundedge = 2;
		x = mitox[mitoid][1];
		y = mitoy[mitoid][1] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][3] - mitoy[mitoid][1]));
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][5] - mitoz[mitoid][1]));
	}
	// Edge 3: Formed by V3,V4,V7 and V8
	if (random_edge >= probability_boundaries[2] && random_edge < probability_boundaries[3])
	{
		foundedge = 3;
		x = mitox[mitoid][2] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][3] - mitox[mitoid][2]));
		y = mitoy[mitoid][2];
		z = mitoz[mitoid][2] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][6] - mitoz[mitoid][2]));
	}
	// Edge 4: Formed by V1,V3,V5 and V7
	if (random_edge >= probability_boundaries[3] && random_edge < probability_boundaries[4])
	{
		foundedge = 4;
		x = mitox[mitoid][0];
		y = mitoy[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][2] - mitoy[mitoid][0]));
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][4] - mitoz[mitoid][0]));
	}
	// Edge 5: Formed by V1,V2,V3 and V4
	if (random_edge >= probability_boundaries[4] && random_edge < probability_boundaries[5])
	{
		foundedge = 5;
		x = mitox[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][1] - mitox[mitoid][0]));
		y = mitoy[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][2] - mitoy[mitoid][0]));
		z = mitoz[mitoid][0];
	}
	// Edge 6: Formed by V5,V6,V7 and V8
	if (random_edge >= probability_boundaries[5] && random_edge < probability_boundaries[6])
	{
		foundedge = 6;
		x = mitox[mitoid][4] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][5] - mitox[mitoid][4]));
		y = mitoy[mitoid][4] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][6] - mitoy[mitoid][4]));
		z = mitoz[mitoid][4];
	}

    fprintf(file, "<xagent>\n");
    fprintf(file, "<name>Receptor</name>\n");
    fprintf(file, "<id>%d</id>\n", Receptor_Agents);
    fprintf(file, "<state>55</state>\n");
    fprintf(file, "<Pending>-1</Pending>\n");
    fprintf(file, "<active>1</active>\n");
    fprintf(file, "<x>%f</x>\n", x);
    fprintf(file, "<y>%f</y>\n", y);
    fprintf(file, "<z>%f</z>\n", z);
    fprintf(file, "<Dc>%f</Dc>\n", Dc_state_51);
	fprintf(file, "<collision_counter_ligand>%f</collision_counter_ligand>\n", collision_counter_ligand);
	fprintf(file, "<check_collision>%d</check_collision>\n", check_collision);
	fprintf(file, "<foundedge>%d</foundedge>\n", foundedge);
	fprintf(file, "<mitoid>%d</mitoid>\n", mitoid);
    fprintf(file, "</xagent>\n");
state55_counter++;
Receptor_Agents++;
}

/* Create Mcl1aBax Receptor Agents */
printf ("Creating Mcl1aBax Receptor Agents\n");
int state21_counter = 0;
while (state21_counter < number_state_21)
{
	double collision_counter_ligand = 0;
	int check_collision = 0;
	int foundedge = 0;
	
	// creates random seed on mitochondria
	int mitoid = mito_number + 1;
	while (mitoid > mito_number || mitoid == mito_number)
	{
	mitoid = (int) rand()/(double)(RAND_MAX)*mito_number; // Rounds down a random double value between 0 and mito_number
	}
	
	double random_edge = rand()/(double)(RAND_MAX);
	// Edge 1: Formed by V1,V2,V5 and V6
	if (random_edge >= probability_boundaries[0] && random_edge < probability_boundaries[1])
	{
		foundedge = 1;
		x = mitox[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][1] - mitox[mitoid][0]));
		y = mitoy[mitoid][0];
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][4] - mitoz[mitoid][0]));
	}
	// Edge 2: Formed by V2,V4,V6 and V8
	if (random_edge >= probability_boundaries[1] && random_edge < probability_boundaries[2])
	{
		foundedge = 2;
		x = mitox[mitoid][1];
		y = mitoy[mitoid][1] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][3] - mitoy[mitoid][1]));
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][5] - mitoz[mitoid][1]));
	}
	// Edge 3: Formed by V3,V4,V7 and V8
	if (random_edge >= probability_boundaries[2] && random_edge < probability_boundaries[3])
	{
		foundedge = 3;
		x = mitox[mitoid][2] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][3] - mitox[mitoid][2]));
		y = mitoy[mitoid][2];
		z = mitoz[mitoid][2] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][6] - mitoz[mitoid][2]));
	}
	// Edge 4: Formed by V1,V3,V5 and V7
	if (random_edge >= probability_boundaries[3] && random_edge < probability_boundaries[4])
	{
		foundedge = 4;
		x = mitox[mitoid][0];
		y = mitoy[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][2] - mitoy[mitoid][0]));
		z = mitoz[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoz[mitoid][4] - mitoz[mitoid][0]));
	}
	// Edge 5: Formed by V1,V2,V3 and V4
	if (random_edge >= probability_boundaries[4] && random_edge < probability_boundaries[5])
	{
		foundedge = 5;
		x = mitox[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][1] - mitox[mitoid][0]));
		y = mitoy[mitoid][0] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][2] - mitoy[mitoid][0]));
		z = mitoz[mitoid][0];
	}
	// Edge 6: Formed by V5,V6,V7 and V8
	if (random_edge >= probability_boundaries[5] && random_edge < probability_boundaries[6])
	{
		foundedge = 6;
		x = mitox[mitoid][4] + (rand()/(double)(RAND_MAX)* (mitox[mitoid][5] - mitox[mitoid][4]));
		y = mitoy[mitoid][4] + (rand()/(double)(RAND_MAX)* (mitoy[mitoid][6] - mitoy[mitoid][4]));
		z = mitoz[mitoid][4];
	}

    fprintf(file, "<xagent>\n");
    fprintf(file, "<name>Receptor</name>\n");
    fprintf(file, "<id>%d</id>\n", Receptor_Agents);
    fprintf(file, "<state>21</state>\n");
    fprintf(file, "<Pending>-1</Pending>\n");
    fprintf(file, "<active>1</active>\n");
    fprintf(file, "<x>%f</x>\n", x);
    fprintf(file, "<y>%f</y>\n", y);
    fprintf(file, "<z>%f</z>\n", z);
    fprintf(file, "<Dc>%f</Dc>\n", Dc_state_21);
	fprintf(file, "<collision_counter_ligand>%f</collision_counter_ligand>\n", collision_counter_ligand);
	fprintf(file, "<check_collision>%d</check_collision>\n", check_collision);
	fprintf(file, "<foundedge>%d</foundedge>\n", foundedge);
	fprintf(file, "<mitoid>%d</mitoid>\n", mitoid);
    fprintf(file, "</xagent>\n");
state21_counter++;
Receptor_Agents++;
}

/* Finalise 0.xml File */
printf("Finalising 0.xml File\n");
printf("Total Mitochondria: %d\n", Mito_Agents);
printf("Total Agents: %d\n", Receptor_Agents+Ligand_Agents);
printf("Total Ligand Agents: %d\n", Ligand_Agents);
printf("Total Receptor Agents: %d\n", Receptor_Agents);

fputs("</states>\n", file);
fclose(file);
return 0;
}
