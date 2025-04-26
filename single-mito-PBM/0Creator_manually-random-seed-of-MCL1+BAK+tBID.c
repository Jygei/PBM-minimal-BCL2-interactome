/*

Jenny Geiger
Institute of Cell Biology and Immunology (IZI) at the University of Stuttgart

FaST generated Codes were changed and modified to simulate a single mitochondrial environment with protein interactions
of a reduced interactome of the BCL-2 family (MCL-1, BAK, tBID)


Code version information: this code creates the 0.xml file with manually set protein distributions (input data resulted from the whole-cell PBM)
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
srand (1);								/************************************** ATTENTION: change between srand(time(NULL)) -> mitos at random positions or srand(1) -> mitos at the same position **************************************/

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
// particle distributions resulted from whole-cell PBM
int MCL1_mito[] = {6180, 6670, 6130, 6450, 6740, 6410, 6590, 6650, 6400, 6440, 6550, 6780, 6370, 5910, 6490, 6260, 6160, 6080, 5830, 6440, 6400, 6250, 5950, 6350, 6250, 6270, 6640, 6460, 6550, 6160, 6390, 5980, 6440, 6610, 6430, 6310, 6440, 6330, 6090, 6190, 240, 100, 280, 190, 220, 290, 160, 190, 170, 220, 170, 190, 160, 190, 180, 130, 230, 200, 130, 270, 180, 270, 230, 160, 270, 190, 150, 180, 190, 100, 150, 200, 220, 180, 170, 200, 80, 220, 250, 210, 200, 220, 240, 160, 180, 200, 180, 230, 300, 190, 400, 130, 230, 180, 210, 330, 290, 180, 250, 190, 230, 300, 200, 180, 220, 170, 180, 200, 180, 230, 160, 300, 210, 190, 200, 190, 250, 200, 160, 210, 250, 200, 220, 180, 240, 180, 210, 280, 190, 260, 260, 300, 230, 190, 200, 240, 200, 180, 160, 180, 210, 290, 160, 230, 290, 260, 200, 270, 250, 170, 210, 250, 230, 240, 190, 260, 220, 190, 220, 200, 240, 220, 200, 180, 210, 170, 140, 230, 180, 220, 220, 230, 250, 180, 130, 210, 210, 190, 260, 240, 170, 310, 200, 220, 250, 290, 240, 160, 270, 170, 240, 200, 210, 180, 230, 180, 240, 230};
int state20_counter = 0;
int m_number = 0;
int mitosize = 2;
while (m_number < mito_number)
{
  	int particles_counter = 1;
	while (particles_counter <= MCL1_mito[m_number])
	{
		double collision_counter_ligand = 0;
		int check_collision = 0;
		int foundedge = 0;

	    int mitoid = m_number;

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

		if (mitowidth[mitoid] == 25428.1)		{ mitosize = 1; }
		else if (mitowidth[mitoid] != 25428.1)	{ mitosize = 0; }

	    fprintf(file, "<xagent>\n");
	    fprintf(file, "<name>Receptor</name>\n");
	    fprintf(file, "<id>%d</id>\n", particles_counter);
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
		fprintf(file, "<mitosize>%d</mitosize>\n", mitosize);
	    fprintf(file, "</xagent>\n");

		//state4_counter++;
		Receptor_Agents++;
        particles_counter++;
	}
	m_number++;
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
printf ("Creating Bax Receptor Agents\n");
// particle distributions resulted from whole-cell PBM
int BAK_mito[] = {52068, 52349, 51909, 52209, 52258, 52459, 52116, 52382, 52337, 52409, 51606, 52065, 52294, 52015, 52208, 51914, 52001, 52245, 52226, 52646, 51993, 51767, 52202, 51910, 52459, 52252, 52210, 52276, 52085, 52006, 52147, 52068, 52570, 52520, 51946, 51967, 52204, 51909, 52306, 52190, 1807, 1816, 1717, 1718, 1713, 1791, 1713, 1780, 1807, 1725, 1682, 1796, 1776, 1792, 1774, 1808, 1804, 1799, 1872, 1770, 1742, 1756, 1663, 1702, 1822, 1748, 1765, 1805, 1751, 1687, 1778, 1729, 1769, 1813, 1756, 1774, 1740, 1733, 1701, 1752, 1734, 1788, 1768, 1772, 1746, 1761, 1745, 1800, 1740, 1747, 1684, 1792, 1747, 1712, 1807, 1775, 1745, 1752, 1774, 1818, 1787, 1765, 1688, 1815, 1698, 1790, 1758, 1773, 1850, 1818, 1798, 1792, 1857, 1707, 1761, 1701, 1797, 1719, 1805, 1753, 1691, 1793, 1799, 1720, 1802, 1706, 1715, 1892, 1828, 1797, 1793, 1813, 1739, 1724, 1718, 1693, 1750, 1861, 1788, 1786, 1716, 1797, 1846, 1760, 1758, 1791, 1799, 1776, 1700, 1811, 1735, 1714, 1761, 1776, 1708, 1757, 1679, 1732, 1786, 1771, 1744, 1731, 1733, 1740, 1680, 1731, 1874, 1698, 1775, 1791, 1790, 1750, 1851, 1753, 1757, 1861, 1699, 1771, 1758, 1736, 1864, 1750, 1704, 1731, 1787, 1789, 1789, 1760, 1828, 1689, 1760, 1727, 1778, 1837, 1766, 1760, 1803, 1784};
int state4_counter = 0;
//int m_number = 0;
//int mitosize = 2;
while (m_number < mito_number)
{
  	int particles_counter = 1;
	while (particles_counter <= BAK_mito[m_number])
	{
		double collision_counter_ligand = 0;
		int check_collision = 0;
		int foundedge = 0;

	    int mitoid = m_number;

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

		if (mitowidth[mitoid] == 25428.1)		{ mitosize = 1; }
		else if (mitowidth[mitoid] != 25428.1)	{ mitosize = 0; }

	    fprintf(file, "<xagent>\n");
	    fprintf(file, "<name>Receptor</name>\n");
	    fprintf(file, "<id>%d</id>\n", particles_counter);
	    fprintf(file, "<state>40</state>\n");
	    fprintf(file, "<Pending>-1</Pending>\n");
	    fprintf(file, "<active>1</active>\n");
	    fprintf(file, "<x>%f</x>\n", x);
	    fprintf(file, "<y>%f</y>\n", y);
	    fprintf(file, "<z>%f</z>\n", z);
	    fprintf(file, "<Dc>%f</Dc>\n", Dc_state_40);
		fprintf(file, "<collision_counter_ligand>%f</collision_counter_ligand>\n", collision_counter_ligand);
		fprintf(file, "<check_collision>%d</check_collision>\n", check_collision);
		fprintf(file, "<foundedge>%d</foundedge>\n", foundedge);
		fprintf(file, "<mitoid>%d</mitoid>\n", mitoid);
		fprintf(file, "<mitosize>%d</mitosize>\n", mitosize);
	    fprintf(file, "</xagent>\n");

		//state4_counter++;
		Receptor_Agents++;
        particles_counter++;
	}
	m_number++;
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
// particle distributions resulted from whole-cell PBM
int tBID_mito[] = {521, 483, 510, 552, 563, 500, 542, 502, 525, 520, 523, 525, 535, 528, 546, 520, 519, 558, 516, 491, 486, 512, 531, 501, 527, 526, 495, 535, 501, 522, 571, 520, 511, 491, 539, 522, 547, 462, 533, 520, 13, 22, 13, 18, 8, 22, 16, 26, 25, 18, 21, 16, 19, 17, 20, 14, 15, 30, 13, 15, 22, 18, 21, 17, 17, 17, 23, 16, 14, 17, 22, 22, 15, 11, 19, 22, 17, 11, 19, 18, 17, 27, 15, 18, 25, 12, 15, 16, 19, 25, 22, 19, 22, 16, 13, 15, 21, 13, 17, 13, 22, 14, 15, 13, 23, 18, 18, 15, 18, 17, 12, 21, 6, 17, 22, 17, 11, 19, 15, 18, 20, 16, 13, 17, 22, 16, 22, 24, 18, 18, 17, 20, 22, 11, 16, 11, 14, 18, 17, 22, 15, 22, 16, 14, 18, 17, 18, 25, 17, 14, 14, 22, 23, 10, 15, 21, 12, 18, 24, 21, 16, 15, 19, 13, 17, 22, 23, 11, 14, 18, 16, 10, 9, 16, 16, 19, 18, 17, 13, 25, 22, 14, 14, 19, 15, 14, 24, 22, 11, 10, 18, 23, 11, 17, 16, 20, 22, 21};
int state56_counter = 0;
while (m_number < mito_number)
{
  	int particles_counter = 1;
	while (particles_counter <= tBID_mito[m_number])
	{
		double collision_counter_ligand = 0;
		int check_collision = 0;
		int foundedge = 0;
	    int mitoid = m_number;

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

		if (mitowidth[mitoid] == 25428.1)		{ mitosize = 1; }
		else if (mitowidth[mitoid] != 25428.1)	{ mitosize = 0; }

	    fprintf(file, "<xagent>\n");
	    fprintf(file, "<name>Receptor</name>\n");
	    fprintf(file, "<id>%d</id>\n", particles_counter);
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
		fprintf(file, "<mitosize>%d</mitosize>\n", mitosize);
	    fprintf(file, "</xagent>\n");

		//state4_counter++;
		Receptor_Agents++;
        particles_counter++;
	}
	m_number++;
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
