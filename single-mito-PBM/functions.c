/*

Jenny Geiger
Institute of Cell Biology and Immunology (IZI) at the University of Stuttgart

FaST generated Codes were changed and modified to simulate a single mitochondrial environment with protein interactions
of a reduced interactome of the BCL-2 family (MCL-1, BAK, tBID)



Code version information:
- this code is for the single-mitochondria particle-based model; main differences to whole cell particle-based model:
	- this code is running without any ligands functions to save memory
	- unbinding of receptors are turned off to simulate only the reactions on the mitochondrial surface area
- new FLAME functions were used: use the newest FLAME v1.5 version!
- periodic boundary conditions of the cell boundaries, otherwise artifacts will occur


problems in the used FLAME v1 (resolved in FLAME v2):
Increasing the buffersize in XMLModelFile.xml changed the results because the random seed is based on the declared
buffersize variable. Normally the RNG seed was defined by the maximum buffersize, so its definition has been redefined
to a buffersize of 1045876 (always used in old simulations).
Workaround: The template of the header file (in "C:\FLAME-GPU\FLAMEGPU\templates") has been changed to 1045876.
Therefore, changes to the buffersize now have no effect on the RNG seed. The redefinition in the templates is marked
with "added by Jenny".

*/


#ifndef _FLAMEGPU_FUNCTIONS
#define _FLAMEGPU_FUNCTIONS
#include <header.h>
#include "AgentVariables.h"
#include "ReactionVariables.h"
#include "globals.h"

#include <stdio.h>
#include <math.h>

#include <stdlib.h>
#include <time.h>


/***********************************************************************************/
/****************************** Mitochondria Functions *****************************/
/***********************************************************************************/

/*********************************/
/****** Mitochondria Output ******/
/*********************************/

__FLAME_GPU_FUNC__ int Mito_Output(xmachine_memory_Mito* xmemory, xmachine_message_Mito_Location_list* Mito_Location_messages, RNG_rand48* rand48){

	// Mito Location data
	int id = xmemory->id;
	double x1 = xmemory->x1;
	double y1 = xmemory->y1;
	double z1 = xmemory->z1;
	double x2 = xmemory->x2;
	double y2 = xmemory->y2;
	double z2 = xmemory->z2;
	double x3 = xmemory->x3;
	double y3 = xmemory->y3;
	double z3 = xmemory->z3;
	double x4 = xmemory->x4;
	double y4 = xmemory->y4;
	double z4 = xmemory->z4;
	double x5 = xmemory->x5;
	double y5 = xmemory->y5;
	double z5 = xmemory->z5;
	double x6 = xmemory->x6;
	double y6 = xmemory->y6;
	double z6 = xmemory->z6;
	double x7 = xmemory->x7;
	double y7 = xmemory->y7;
	double z7 = xmemory->z7;
	double x8 = xmemory->x8;
	double y8 = xmemory->y8;
	double z8 = xmemory->z8;
	double mitocentrex = xmemory->mitocentrex;
	double mitocentrey = xmemory->mitocentrey;
	double mitocentrez = xmemory->mitocentrez;
	double mitolength = xmemory->mitolength;
	double mitoheight = xmemory->mitoheight;
	double mitowidth = xmemory->mitowidth;
	double radius = xmemory->radius;
	double radius2 = xmemory->radius2;
	double mito_collision_counter = xmemory->mito_collision_counter;
	double mito_energetics = xmemory->mito_energetics;
	int current_mitosize = xmemory->current_mitosize;
	
	add_Mito_Location_message(Mito_Location_messages, xmemory->id, xmemory->x1, xmemory->y1, xmemory->z1, xmemory->x2, xmemory->y2, xmemory->z2, xmemory->x3, xmemory->y3, xmemory->z3, xmemory->x4, xmemory->y4, xmemory->z4, xmemory->x5, xmemory->y5, xmemory->z5, xmemory->x6, xmemory->y6, xmemory->z6, xmemory->x7, xmemory->y7, xmemory->z7, xmemory->x8, xmemory->y8, xmemory->z8, xmemory->mitocentrex, xmemory->mitocentrey, xmemory->mitocentrez, xmemory->mitolength, xmemory->mitoheight, xmemory->mitowidth, xmemory->radius, xmemory->radius2, xmemory->mito_collision_counter, xmemory->mito_energetics, xmemory->current_mitosize);
	return 0;
}


/*****************************************************************************/
/***************************** Mito_CC function ****************************/
/*****************************************************************************/
/* function to successfully load mitochondrial environment; function makes nothing else */

__FLAME_GPU_FUNC__ int Mito_CC(xmachine_memory_Mito* xmemory, RNG_rand48* rand48)
{

    return 0;
}
	
/*****************************************************************************/
/****************************** Ligand Functions *****************************/
/*****************************************************************************/
/* no ligand functions needed in single-mitochondria PBM */

/***************************/
/****** Ligand Output ******/
/***************************/



/*****************************************************************************/
/***************************** Receptor Functions ****************************/
/*****************************************************************************/

/***************************/
/***** Receptor Output *****/
/***************************/

__FLAME_GPU_FUNC__ int Receptor_Output(xmachine_memory_Receptor* xmemory, xmachine_message_Mito_Location_list* Mito_Location_messages, xmachine_message_Receptor_Location_list* Receptor_Location_messages, RNG_rand48* rand48)
{
	int mito_counter, mito_counter2;
	mito_counter = xmemory->mitoid;			// Counter keeps the mito_number with which it collided in the Ligand_output fct
	mito_counter2 = 0;						// Counter to go through all mitos
	
	/* cornerpoints of the Mito planes */
	double x1[mito_number], y1[mito_number], z1[mito_number], x2[mito_number], y2[mito_number], z2[mito_number], x3[mito_number], y3[mito_number], z3[mito_number], x4[mito_number], y4[mito_number], z4[mito_number], x5[mito_number], y5[mito_number], z5[mito_number], x6[mito_number], y6[mito_number], z6[mito_number], x7[mito_number], y7[mito_number], z7[mito_number], x8[mito_number], y8[mito_number], z8[mito_number];
	double mito_energetics[mito_number];
	int current_mitosize[mito_number];

	/* Cycles through all Mito Locations and safe the Mito-cornerpoints */
	xmachine_message_Mito_Location* current_message = get_first_Mito_Location_message(Mito_Location_messages);

	while (current_message)
	{
		x1[current_message->id] = current_message->x1;
		y1[current_message->id] = current_message->y1;
		z1[current_message->id] = current_message->z1;
		x2[current_message->id] = current_message->x2;
		y2[current_message->id] = current_message->y2;
		z2[current_message->id] = current_message->z2;
		x3[current_message->id] = current_message->x3;
		y3[current_message->id] = current_message->y3;
		z3[current_message->id] = current_message->z3;
		x4[current_message->id] = current_message->x4;
		y4[current_message->id] = current_message->y4;
		z4[current_message->id] = current_message->z4;
		x5[current_message->id] = current_message->x5;
		y5[current_message->id] = current_message->y5;
		z5[current_message->id] = current_message->z5;
		x6[current_message->id] = current_message->x6;
		y6[current_message->id] = current_message->y6;
		z6[current_message->id] = current_message->z6;
		x7[current_message->id] = current_message->x7;
		y7[current_message->id] = current_message->y7;
		z7[current_message->id] = current_message->z7;
		x8[current_message->id] = current_message->x8;
		y8[current_message->id] = current_message->y8;
		z8[current_message->id] = current_message->z8;
		mito_energetics[current_message->id] = current_message->mito_energetics;
		current_mitosize[current_message->id] = current_message->current_mitosize;
		
		current_message = get_next_Mito_Location_message(current_message, Mito_Location_messages);
		
		mito_counter2++;
	}
	
	if (xmemory->foundedge != 0)
	{

		/* Diffusion of Mcl-1 on the mitochondria surface */
		/*********************************************************************************************************/
		/*************************************** Diffusion on Mito surface ***************************************/
		/*********************************************************************************************************/

			/* first locations of Receptor */
			double x_old0, y_old0, z_old0;
			x_old0 = xmemory->x;
			y_old0 = xmemory->y;
			z_old0 = xmemory->z;
			
			int counter = 0;
			while (counter < brownianloop)
			{
				/* previous locations of Receptor */
				double x_old, y_old, z_old;
				x_old = xmemory->x;
				y_old = xmemory->y;
				z_old = xmemory->z;

				/* Assign deltax, deltay and deltaz (in nm) */
				double n1, n2, n3;
				double a1, b1, r1;
				
				do {
					a1 = 2.0*rnd<CONTINUOUS>(rand48) -1;
					b1 = 2.0*rnd<CONTINUOUS>(rand48) -1;
					r1 = a1*a1 + b1*b1;
				} while (r1 == 0.0 || r1 > 1.0); /* Check that x^2 + y^2 is less than 1 */
				{
					double d = sqrt(-2.0*log(r1) / r1);
					n1 = a1*d;
					n2 = b1*d;
				}
				double a2, b2, r2;
				do {
					a2 = 2.0*rnd<CONTINUOUS>(rand48) -1;
					b2 = 2.0*rnd<CONTINUOUS>(rand48) -1;
					r2 = a2*a2 + b2*b2;
				} while (r2 == 0.0 || r2 > 1.0); /* Check that x^2 + y^2 is less than 1 */
				{
					double d = sqrt(-2.0*log(r2) / r2);
					n3 = a2*d;
				}
				double a3, b3, r3;
				do {
					a3 = 2.0*rnd<CONTINUOUS>(rand48) -1;
					b3 = 2.0*rnd<CONTINUOUS>(rand48) -1;
					r3 = a3*a3 + b3*b3;
				} while (r3 == 0.0 || r3 > 1.0); /* Check that x^2 + y^2 is less than 1 */
				{
					double d = sqrt(-2.0*log(r3) / r3);
					n2 = a3*d;
				}
				
				/* walk distance */
				double deltax = xmemory->Dc * n1;
				double deltay = xmemory->Dc * n2;
				double deltaz = xmemory->Dc * n3;
				
				/* Update x, y and z Variables */
				xmemory->x += deltax;
				xmemory->y += deltay;
				xmemory->z += deltaz;
				
				/* calculation of the walk distances until its crossing a border/plane */
				double fractionx1 = (x_old - x1[mito_counter])/(x_old - xmemory->x);
				double fractionx2 = (x_old - x2[mito_counter])/(x_old - xmemory->x);
				double fractiony1 = (y_old - y1[mito_counter])/(y_old - xmemory->y);
				double fractiony3 = (y_old - y3[mito_counter])/(y_old - xmemory->y);
				double fractionz1 = (z_old - z1[mito_counter])/(z_old - xmemory->z);
				double fractionz5 = (z_old - z5[mito_counter])/(z_old - xmemory->z);
				
				/* calculation of the walking distances into positive numbers */
				double deltax_pos = deltax;
				double deltay_pos = deltay;
				double deltaz_pos = deltaz;
				if (deltax < 0)	{	deltax_pos = fabs(deltax);	}
				if (deltay < 0)	{	deltay_pos = fabs(deltay);	}
				if (deltaz < 0)	{	deltaz_pos = fabs(deltaz);	}
				
				/* correction of the constant variable for each plane. If Mcl-1 is on a plane it keeps one variable constant */
				if (xmemory->foundedge == 1)		{	xmemory->y = y_old0;		}
				if (xmemory->foundedge == 3)		{	xmemory->y = y_old0;		}
				if (xmemory->foundedge == 2)		{	xmemory->x = x_old0;		}
				if (xmemory->foundedge == 4)		{	xmemory->x = x_old0;		}
				if (xmemory->foundedge == 5)		{	xmemory->z = z_old0;		}
				if (xmemory->foundedge == 6)		{	xmemory->z = z_old0;		}
				
				/* counts how many borders of the mito are crossed */
				int borders_crossed = 0;
				
				/* Check on which plane is the ligand -> then checking if its walking outside of the plane -> then checking into which plane its walking and updating the position */
				// Mcl-1 on plane 1 or 3
				/* Mcl-1 on plane 1 or 3		->		Mcl-1 stays in the borders of plane 1 or 3 	*/
				if ((borders_crossed == 0 && (xmemory->foundedge == 1)) || (borders_crossed == 0 && (xmemory->foundedge == 3)))
				{ 
					if ((xmemory->x >= x1[mito_counter]) && (xmemory->x <= x2[mito_counter]) && (xmemory->z >= z1[mito_counter]) && (xmemory->z <= z5[mito_counter]))
					{
						if (xmemory->foundedge == 1)		{	xmemory->y = y1[mito_counter];		}
						if (xmemory->foundedge == 3)		{	xmemory->y = y3[mito_counter];		}
						borders_crossed++;
					}
					else
					{
						/* Mcl-1 on plane 1 or 3		->		walks on to plane 4	*/
						if (borders_crossed == 0 && (xmemory->x < x1[mito_counter]) && (xmemory->z >= z1[mito_counter]) && (xmemory->z <= z5[mito_counter]))
						{
							if (xmemory->foundedge == 1)		{	xmemory->y = y1[mito_counter] + fractionx1 * deltax_pos;		}
							if (xmemory->foundedge == 3)		{	xmemory->y = y3[mito_counter] - fractionx1 * deltax_pos;		}
							xmemory->foundedge = 4;
							xmemory->x = x1[mito_counter];
							borders_crossed++;
						}
						/* Mcl-1 on plane 1 or 3		->		walks on to plane 2	*/
						if (borders_crossed == 0 && (xmemory->x > x2[mito_counter]) && (xmemory->z >= z1[mito_counter]) && (xmemory->z <= z5[mito_counter]))
						{
							if (xmemory->foundedge == 1)		{	xmemory->y = y1[mito_counter] + fractionx2 * deltax_pos;		}
							if (xmemory->foundedge == 3)		{	xmemory->y = y3[mito_counter] - fractionx2 * deltax_pos;		}
							xmemory->foundedge = 2;
							xmemory->x = x2[mito_counter];
							borders_crossed++;
						}
						/* Mcl-1 on plane 1 or 3		->		walks on to plane 5	*/
						if (borders_crossed == 0 && (xmemory->z < z1[mito_counter]) && (xmemory->x >= x1[mito_counter]) && (xmemory->x <= x2[mito_counter]))
						{
							if (xmemory->foundedge == 1)		 {	xmemory->y = y1[mito_counter] + fractionz1 * deltaz_pos;		}
							if (xmemory->foundedge == 3)		 {	xmemory->y = y3[mito_counter] - fractionz1 * deltaz_pos;		}
							xmemory->foundedge = 5;
							xmemory->z = z1[mito_counter];
							borders_crossed++;
						}
						/* Mcl-1 on plane 1 or 3		->		walks on to plane 6	*/
						if (borders_crossed == 0 && (xmemory->z > z5[mito_counter]) && (xmemory->x >= x1[mito_counter]) && (xmemory->x <= x2[mito_counter]))
						{
							if (xmemory->foundedge == 1)		{	xmemory->y = y1[mito_counter] + fractionz5 * deltaz_pos;		}
							if (xmemory->foundedge == 3)		{	xmemory->y = y3[mito_counter] - fractionz5 * deltaz_pos;		}
							xmemory->foundedge = 6;
							xmemory->z = z5[mito_counter];
							borders_crossed++;
						}
						/* Mcl-1 on plane 1 or 3		->		walks on to plane 2/4/5/6	*/
						if (borders_crossed == 0 && (xmemory->x < x1[mito_counter]) && (xmemory->z < z1[mito_counter]))
						{
							if (fractionx1 < fractionz1)
							{
								if (xmemory->foundedge == 1)	{	xmemory->y = y1[mito_counter] + fractionz1 * deltaz_pos;	}
								if (xmemory->foundedge == 3)	{	xmemory->y = y3[mito_counter] - fractionz1 * deltaz_pos;	}
								xmemory->x = x1[mito_counter] + fractionx1 * deltax_pos;
								xmemory->z = z1[mito_counter];
								xmemory->foundedge = 5;
								borders_crossed++;
							}
							if (fractionx1 >= fractionz1)
							{
								if (xmemory->foundedge == 1)	{	xmemory->y = y1[mito_counter] + fractionx1 * deltax_pos;	}
								if (xmemory->foundedge == 3)	{	xmemory->y = y3[mito_counter] - fractionx1 * deltax_pos;	}
								xmemory->x = x1[mito_counter];
								xmemory->z = z1[mito_counter] + fractionz1 * deltaz_pos;
								xmemory->foundedge = 4;
								borders_crossed++;
							}
						}
						if (borders_crossed == 0 && (xmemory->x < x1[mito_counter]) && (xmemory->z > z5[mito_counter]))
						{
							if (fractionx1 < fractionz5)
							{
								if (xmemory->foundedge == 1)	{	xmemory->y = y1[mito_counter] + fractionz5 * deltaz_pos;	}
								if (xmemory->foundedge == 3)	{	xmemory->y = y3[mito_counter] - fractionz5 * deltaz_pos;	}
								xmemory->x = x1[mito_counter] + fractionx1 * deltax_pos;
								xmemory->z = z5[mito_counter];
								xmemory->foundedge = 6;
								borders_crossed++;
							}
							if (fractionx1 >= fractionz5)
							{
								if (xmemory->foundedge == 1)	{	xmemory->y = y1[mito_counter] + fractionx1 * deltax_pos;	}
								if (xmemory->foundedge == 3)	{	xmemory->y = y3[mito_counter] - fractionx1 * deltax_pos;	}
								xmemory->x = x1[mito_counter];
								xmemory->z = z5[mito_counter] - fractionz5 * deltaz_pos;
								xmemory->foundedge = 4;
								borders_crossed++;
							}
						}
						if (borders_crossed == 0 && (xmemory->x > x2[mito_counter]) && (xmemory->z < z1[mito_counter]))
						{
							if (fractionx2 < fractionz1)
							{
								if (xmemory->foundedge == 1)	{	xmemory->y = y1[mito_counter] + fractionz1 * deltaz_pos;	}
								if (xmemory->foundedge == 3)	{	xmemory->y = y3[mito_counter] - fractionz1 * deltaz_pos;	}
								xmemory->x = x2[mito_counter] - fractionx2 * deltax_pos;
								xmemory->z = z1[mito_counter];
								xmemory->foundedge = 5;
								borders_crossed++;
							}
							if (fractionx2 >= fractionz1)
							{
								if (xmemory->foundedge == 1)	{	xmemory->y = y1[mito_counter] + fractionx2 * deltax_pos;	}
								if (xmemory->foundedge == 3)	{	xmemory->y = y3[mito_counter] - fractionx2 * deltax_pos;	}
								xmemory->x = x2[mito_counter];
								xmemory->z = z1[mito_counter] + fractionz1 * deltaz_pos;
								xmemory->foundedge = 2;
								borders_crossed++;
							}
						}
						if (borders_crossed == 0 && (xmemory->x > x2[mito_counter]) && (xmemory->z > z5[mito_counter]))
						{
							if (fractionx2 < fractionz5)
							{
								if (xmemory->foundedge == 1)	{	xmemory->y = y1[mito_counter] + fractionz5 * deltaz_pos;	}
								if (xmemory->foundedge == 3)	{	xmemory->y = y3[mito_counter] - fractionz5 * deltaz_pos;	}
								xmemory->x = x2[mito_counter] - fractionx2 * deltax_pos;
								xmemory->z = z5[mito_counter];
								xmemory->foundedge = 6;
								borders_crossed++;
							}
							if (fractionx2 >= fractionz5)
							{
								if (xmemory->foundedge == 1)	{	xmemory->y = y1[mito_counter] + fractionx2 * deltax_pos;	}
								if (xmemory->foundedge == 3)	{	xmemory->y = y3[mito_counter] - fractionx2 * deltax_pos;	}
								xmemory->x = x2[mito_counter];
								xmemory->z = z5[mito_counter] - fractionz5 * deltaz_pos;
								xmemory->foundedge = 2;
								borders_crossed++;
							}
						}
					}
				}
				
				//	Mcl-1 on plane 2 or 4
				/* Mcl-1 on plane 2 or 4		->		Mcl-1 stays in the borders of plane 2 or 4 	*/
				if ((borders_crossed == 0 && (xmemory->foundedge == 2)) || (borders_crossed == 0 && (xmemory->foundedge == 4)))
				{
					if ((xmemory->y >= y1[mito_counter]) && (xmemory->y <= y3[mito_counter]) && (xmemory->z >= z1[mito_counter]) && (xmemory->z <= z5[mito_counter]))
					{
						if (xmemory->foundedge == 2)		{	xmemory->x = x2[mito_counter];		}
						if (xmemory->foundedge == 4)		{	xmemory->x = x1[mito_counter];		}
						borders_crossed++;
					}
					else
					{
						/* Mcl-1 on plane 2 or 4		->		walks on to plane 1	*/
						if (borders_crossed == 0 && (xmemory->y < y1[mito_counter]) && (xmemory->z >= z1[mito_counter]) && (xmemory->z <= z5[mito_counter]))
						{
							if (xmemory->foundedge == 2)		{	xmemory->x = x2[mito_counter] - fractiony1 * deltay_pos;		}
							if (xmemory->foundedge == 4)		{	xmemory->x = x1[mito_counter] + fractiony1 * deltay_pos;		}
							xmemory->foundedge = 1;
							xmemory->y = y1[mito_counter];
							borders_crossed++;
						}
						/* Mcl-1 on plane 2 or 4		->		walks on to plane 3	*/
						if (borders_crossed == 0 && (xmemory->y > y3[mito_counter]) && (xmemory->z >= z1[mito_counter]) && (xmemory->z <= z5[mito_counter]))
						{
							if (xmemory->foundedge == 2)		{	xmemory->x = x2[mito_counter] - fractiony3 * deltay_pos;		}
							if (xmemory->foundedge == 4)		{	xmemory->x = x1[mito_counter] + fractiony3 * deltay_pos;		}
							xmemory->foundedge = 3;
							xmemory->y = y3[mito_counter];
							borders_crossed++;
						}
						/* Mcl-1 on plane 2 or 4		->		walks on to plane 5	*/
						if (borders_crossed == 0 && (xmemory->z < z1[mito_counter]) && (xmemory->y >= y1[mito_counter]) && (xmemory->y <= y3[mito_counter]))
						{
							if (xmemory->foundedge == 2)		{	xmemory->x = x2[mito_counter] - fractionz1 * deltaz_pos;		}
							if (xmemory->foundedge == 4)		{	xmemory->x = x1[mito_counter] + fractionz1 * deltaz_pos;		}
							xmemory->foundedge = 5;
							xmemory->z = z1[mito_counter];
							borders_crossed++;
						}
						/* Mcl-1 on plane 2 or 4		->		walks on to plane 6	*/
						if (borders_crossed == 0 && (xmemory->z > z5[mito_counter]) && (xmemory->y >= y1[mito_counter]) && (xmemory->y <= y3[mito_counter]))
						{
							if (xmemory->foundedge == 2)		{	xmemory->x = x2[mito_counter] - fractionz5 * deltaz_pos;		}
							if (xmemory->foundedge == 4)		{	xmemory->x = x1[mito_counter] + fractionz5 * deltaz_pos;		}
							xmemory->foundedge = 6;
							xmemory->z = z5[mito_counter];
							borders_crossed++;
						}
						
						/* Mcl-1 on plane 2 or 4		->		walks on to plane 1/3/5/6	*/
						if (borders_crossed == 0 && (xmemory->y < y1[mito_counter]) && (xmemory->z < z1[mito_counter]))
						{
							if (fractiony1 < fractionz1)
							{
								if (xmemory->foundedge == 2)	{	xmemory->x = x1[mito_counter] + fractionz1 * deltaz_pos;	}
								if (xmemory->foundedge == 4)	{	xmemory->x = x2[mito_counter] - fractionz1 * deltaz_pos;	}
								xmemory->y = y1[mito_counter] + fractiony1 * deltay_pos;
								xmemory->z = z1[mito_counter];
								xmemory->foundedge = 5;
								borders_crossed++;
							}
							if (fractiony1 >= fractionz1)
							{
								if (xmemory->foundedge == 2)	{	xmemory->x = x1[mito_counter] + fractiony1 * deltay_pos;	}
								if (xmemory->foundedge == 4)	{	xmemory->x = x2[mito_counter] - fractiony1 * deltay_pos;	}
								xmemory->y = y1[mito_counter];
								xmemory->z = z1[mito_counter] + fractionz1 * deltaz_pos;
								xmemory->foundedge = 1;
								borders_crossed++;
							}
						}
						if (borders_crossed == 0 && (xmemory->y < y1[mito_counter]) && (xmemory->z > z5[mito_counter]))
						{
							if (fractiony1 < fractionz5)
							{
								if (xmemory->foundedge == 2)	{	xmemory->x = x1[mito_counter] + fractionz5 * deltaz_pos;	}
								if (xmemory->foundedge == 4)	{	xmemory->x = x2[mito_counter] - fractionz5 * deltaz_pos;	}
								xmemory->y = y1[mito_counter] + fractiony1 * deltay_pos;
								xmemory->z = z5[mito_counter];
								xmemory->foundedge = 6;
								borders_crossed++;
							}
							if (fractiony1 >= fractionz5)
							{
								if (xmemory->foundedge == 2)	{	xmemory->x = x1[mito_counter] + fractiony1 * deltay_pos;	}
								if (xmemory->foundedge == 4)	{	xmemory->x = x2[mito_counter] - fractiony1 * deltay_pos;	}
								xmemory->y = y1[mito_counter];
								xmemory->z = z5[mito_counter] - fractionz5 * deltaz_pos;
								xmemory->foundedge = 1;
								borders_crossed++;
							}
						}
						if (borders_crossed == 0 && (xmemory->y > y3[mito_counter]) && (xmemory->z < z1[mito_counter]))
						{
							if (fractiony3 < fractionz1)
							{
								if (xmemory->foundedge == 2)	{	xmemory->x = x1[mito_counter] + fractionz1 * deltaz_pos;	}
								if (xmemory->foundedge == 4)	{	xmemory->x = x2[mito_counter] - fractionz1 * deltaz_pos;	}
								xmemory->y = y3[mito_counter] - fractiony3 * deltay_pos;
								xmemory->z = z1[mito_counter];
								xmemory->foundedge = 5;
								borders_crossed++;
							}
							if (fractiony3 >= fractionz1)
							{
								if (xmemory->foundedge == 2)	{	xmemory->x = x1[mito_counter] + fractiony3 * deltay_pos;	}
								if (xmemory->foundedge == 4)	{	xmemory->x = x2[mito_counter] - fractiony3 * deltay_pos;	}
								xmemory->y = y3[mito_counter];
								xmemory->z = z1[mito_counter] + fractionz1 * deltaz_pos;
								xmemory->foundedge = 3;
								borders_crossed++;
							}
						}
						if (borders_crossed == 0 && (xmemory->y > y3[mito_counter]) && (xmemory->z > z5[mito_counter]))
						{
							if (fractiony3 < fractionz5)
							{
								if (xmemory->foundedge == 2)	{	xmemory->x = x1[mito_counter] + fractionz5 * deltaz_pos;	}
								if (xmemory->foundedge == 4)	{	xmemory->x = x2[mito_counter] - fractionz5 * deltaz_pos;	}
								xmemory->y = y3[mito_counter] - fractiony3 * deltay_pos;
								xmemory->z = z5[mito_counter];
								xmemory->foundedge = 6;
							}
							if (fractiony3 >= fractionz5)
							{
								if (xmemory->foundedge == 2)	{	xmemory->x = x1[mito_counter] + fractiony3 * deltay_pos;	}
								if (xmemory->foundedge == 4)	{	xmemory->x = x2[mito_counter] - fractiony3 * deltay_pos;	}
								xmemory->y = y3[mito_counter];
								xmemory->z = z5[mito_counter] - fractionz5 * deltaz_pos;
								xmemory->foundedge = 3;
								borders_crossed++;
							}
						}
					}
				}
				
				// Mcl-1 on plane 5 or 6	
				/* Mcl-1 on plane 5 or 6		->		Mcl-1 stays in the borders of plane 5 or 5 	*/
				if ((borders_crossed == 0 && (xmemory->foundedge == 5)) || (borders_crossed == 0 && (xmemory->foundedge == 6)))
				{
					if ((xmemory->x >= x1[mito_counter]) && (xmemory->x <= x2[mito_counter]) && (xmemory->y >= y1[mito_counter]) && (xmemory->y <= y3[mito_counter]))
					{
						if (xmemory->foundedge == 5)		{	xmemory->z = z1[mito_counter];		}
						if (xmemory->foundedge == 6)		{	xmemory->z = z5[mito_counter];		}
						borders_crossed++;
					}
					
					else
					{
						/* Mcl-1 on plane 5 or 6		->		walks on to plane 4	*/
						if (borders_crossed == 0 && (xmemory->x < x1[mito_counter]) && (xmemory->y >= y1[mito_counter]) && (xmemory->y <= y3[mito_counter]))
						{
							if (xmemory->foundedge == 5)		{	xmemory->z = z1[mito_counter] + fractionx1 * deltax_pos;		}
							if (xmemory->foundedge == 6)		{	xmemory->z = z5[mito_counter] - fractionx1 * deltax_pos;		}
							xmemory->foundedge = 4;
							xmemory->x = x1[mito_counter];
							borders_crossed++;
						}
						/* Mcl-1 on plane 5 or 6		->		walks on to plane 2	*/
						if (borders_crossed == 0 && (xmemory->x > x2[mito_counter]) && (xmemory->y >= y1[mito_counter]) && (xmemory->y <= y3[mito_counter]))
						{
							if (xmemory->foundedge == 5)		{	xmemory->z = z1[mito_counter] + fractionx2 * deltax_pos;		}
							if (xmemory->foundedge == 6)		{	xmemory->z = z5[mito_counter] - fractionx2 * deltax_pos;		}
							xmemory->foundedge = 2;
							xmemory->x = x2[mito_counter];
							borders_crossed++;
						}
						/* Mcl-1 on plane 5 or 6		->		walks on to plane 1	*/
						if (borders_crossed == 0 && (xmemory->y < y1[mito_counter]) && (xmemory->x >= x1[mito_counter]) && (xmemory->x <= x2[mito_counter]))
						{
							if (xmemory->foundedge == 5)		{	xmemory->z = z1[mito_counter] + fractiony1 * deltay_pos;		}
							if (xmemory->foundedge == 6)		{	xmemory->z = z5[mito_counter] - fractiony1 * deltay_pos;		}
							xmemory->foundedge = 1;
							xmemory->y = y1[mito_counter];
							borders_crossed++;
						}
						/* Mcl-1 on plane 5 or 6		->		walks on to plane 3	*/
						if (borders_crossed == 0 && (xmemory->y > y3[mito_counter]) && (xmemory->x >= x1[mito_counter]) && (xmemory->x <= x2[mito_counter]))
						{
							if (xmemory->foundedge == 5)		{	xmemory->z = z1[mito_counter] + fractiony3 * deltay_pos;		}
							if (xmemory->foundedge == 6)		{	xmemory->z = z5[mito_counter] - fractiony3 * deltay_pos;		}
							xmemory->foundedge = 3;
							xmemory->y = y3[mito_counter];
							borders_crossed++;
						}
						
						/* Mcl-1 on plane 5 or 6		->		walks on to plane 1/2/3/4	*/
						if (borders_crossed == 0 && (xmemory->x < x1[mito_counter]) && (xmemory->y < y1[mito_counter]))
						{
							if (fractionx1 < fractiony1)
							{
								if (xmemory->foundedge == 5)	{	xmemory->z = z1[mito_counter] + fractiony1 * deltay_pos;	}
								if (xmemory->foundedge == 6)	{	xmemory->z = z5[mito_counter] - fractiony1 * deltay_pos;	}
								xmemory->x = x1[mito_counter] + fractionx1 * deltax_pos;
								xmemory->y = y1[mito_counter];
								xmemory->foundedge = 1;
								borders_crossed++;
							}
							if (fractionx1 >= fractiony1)
							{
								if (xmemory->foundedge == 5)	{	xmemory->z = z1[mito_counter] + fractionx1 * deltax_pos;	}
								if (xmemory->foundedge == 6)	{	xmemory->z = z5[mito_counter] - fractionx1 * deltax_pos;	}
								xmemory->x = x1[mito_counter];
								xmemory->y = y1[mito_counter] + fractiony1 * deltay_pos;
								xmemory->foundedge = 4;
								borders_crossed++;
							}
						}
						if (borders_crossed == 0 && (xmemory->x < x1[mito_counter]) && (xmemory->y > y3[mito_counter]))
						{
							if (fractionx1 < fractiony3)
							{
								if (xmemory->foundedge == 5)	{	xmemory->z = z1[mito_counter] + fractiony3 * deltay_pos;	}
								if (xmemory->foundedge == 6)	{	xmemory->z = z5[mito_counter] - fractiony3 * deltay_pos;	}
								xmemory->x = x1[mito_counter] + fractionx1 * deltax_pos;
								xmemory->y = y3[mito_counter];
								xmemory->foundedge = 3;
								borders_crossed++;
							}
							if (fractionx1 >= fractiony3)
							{
								if (xmemory->foundedge == 5)	{	xmemory->z = z1[mito_counter] + fractionx1 * deltax_pos;	}
								if (xmemory->foundedge == 6)	{	xmemory->z = z5[mito_counter] - fractionx1 * deltax_pos;	}
								xmemory->x = x1[mito_counter];
								xmemory->y = y3[mito_counter] - fractiony3 * deltay_pos;
								xmemory->foundedge = 4;
								borders_crossed++;
							}
						}
						if (borders_crossed == 0 && (xmemory->x > x2[mito_counter]) && (xmemory->y < y1[mito_counter]))
						{
							if (fractionx2 < fractiony1)
							{
								if (xmemory->foundedge == 5)	{	xmemory->z = z1[mito_counter] + fractiony1 * deltay_pos;	}
								if (xmemory->foundedge == 6)	{	xmemory->z = z5[mito_counter] - fractiony1 * deltay_pos;	}
								xmemory->x = x2[mito_counter] - fractionx2 * deltax_pos;
								xmemory->y = y1[mito_counter];
								xmemory->foundedge = 1;
								borders_crossed++;
							}
							if (fractionx2 >= fractiony1)
							{
								if (xmemory->foundedge == 5)	{	xmemory->z = z1[mito_counter] + fractionx2 * deltax_pos;	}
								if (xmemory->foundedge == 6)	{	xmemory->z = z5[mito_counter] - fractionx2 * deltax_pos;	}
								xmemory->x = x2[mito_counter];
								xmemory->y = y1[mito_counter] + fractiony1 * deltay_pos;
								xmemory->foundedge = 2;
								borders_crossed++;
							}
						}
						if (borders_crossed == 0 && (xmemory->x > x2[mito_counter]) && (xmemory->y > y3[mito_counter]))
						{
							if (fractionx2 < fractiony3)
							{
								if (xmemory->foundedge == 5)	{	xmemory->z = z1[mito_counter] + fractiony3 * deltay_pos;	}
								if (xmemory->foundedge == 6)	{	xmemory->z = z5[mito_counter] - fractiony3 * deltay_pos;	}
								xmemory->x = x2[mito_counter] - fractionx2 * deltax_pos;
								xmemory->y = y3[mito_counter];
								xmemory->foundedge = 3;
								borders_crossed++;
							}
							if (fractionx2 >= fractiony3)
							{
								if (xmemory->foundedge == 5)	{	xmemory->z = z1[mito_counter] + fractionx2 * deltax_pos;	}
								if (xmemory->foundedge == 6)	{	xmemory->z = z5[mito_counter] - fractionx2 * deltax_pos;	}
								xmemory->x = x2[mito_counter];
								xmemory->y = y3[mito_counter] - fractiony3 * deltay_pos;
								xmemory->foundedge = 2;
								borders_crossed++;
							}
						}
					}
				}
				counter++;
			}
		//}
	}
        xmemory->active = 1;
		/* Add Receptor Location Message */
			add_Receptor_Location_message (Receptor_Location_messages, xmemory->id, xmemory->state, xmemory->x,  xmemory->y,  xmemory->z, xmemory->collision_counter_ligand, xmemory->check_collision, xmemory->foundedge, xmemory->mitoid, xmemory->mitosize);
			
return 0;
}


/**********************************/
/*** Receptor Dissociate Receptor ***/
/**********************************/

__FLAME_GPU_FUNC__ int Receptor_Dissociate_Receptor(xmachine_memory_Receptor* xmemory, xmachine_memory_Receptor_list* Receptor_agents, xmachine_message_Mito_Location_list* Mito_Location_messages, RNG_rand48* rand48) {
	/* general changes in this function (deviating from FaST code):
	- dx, dy, dz were changed that they change accordingly to the mitochondrial plane
	- after dissociation the new Receptor agent gets a new unique ID (FLAME function is only running with the newest FLAME v1.5 version)
	*/

	xmemory->active = 1;
	/* cornerpoints of the Mito planes */
	double x1, x2, x3, x4, x5, x6, x7, x8, y1, y2, y3, y4, y5, y6 ,y7, y8, z1, z2, z3, z4, z5, z6, z7, z8;
	int current_mitosize;

	/* Cycles through all Mito Locations and safe the Mito-cornerpoints */
	xmachine_message_Mito_Location* current_message = get_first_Mito_Location_message(Mito_Location_messages);

	while (current_message)
	{
		x1 = current_message->x1;
		y1 = current_message->y1;
		z1 = current_message->z1;
		x2 = current_message->x2;
		y2 = current_message->y2;
		z2 = current_message->z2;
		x3 = current_message->x3;
		y3 = current_message->y3;
		z3 = current_message->z3;
		x4 = current_message->x4;
		y4 = current_message->y4;
		z4 = current_message->z4;
		x5 = current_message->x5;
		y5 = current_message->y5;
		z5 = current_message->z5;
		x6 = current_message->x6;
		y6 = current_message->y6;
		z6 = current_message->z6;
		x7 = current_message->x7;
		y7 = current_message->y7;
		z7 = current_message->z7;
		x8 = current_message->x8;
		y8 = current_message->y8;
		z8 = current_message->z8;

		//mito_energetics[current_message->id] = current_message->mito_energetics;
		current_mitosize = current_message->current_mitosize;

		current_message = get_next_Mito_Location_message(current_message, Mito_Location_messages);

	}

	int id_new;
	id_new = generate_Receptor_id();

	double dx, dy, dz;
	double theta = 2 * pi * rnd<CONTINUOUS>(rand48);
	double random = rnd<CONTINUOUS>(rand48);
	double x = xmemory->x;
	double y = xmemory->y;
	double z = xmemory->z;
	int new_agent_foundedge = xmemory->foundedge;
	double reaction = 0;
	double new_agent_state, new_agent_Dc_state, ubr_reaction;

	/* Dissociation of aBak2 into Bak and aBak */
	if (xmemory->state == 42 && xmemory->active == 1)
	{
		if (random < kr_reaction_1)
		{
			xmemory->state = 4;
			reaction = 1;
			new_agent_state = 41;
			new_agent_Dc_state = Dc_state_41;
			ubr_reaction = ubr_reaction_1;
			int active = 1;
			xmemory->active = 0;
			xmemory->Dc = Dc_state_4;
		}
	}

	/* Dissociation of aBak2 into aBak and aBak */
	if (xmemory->state == 42 && xmemory->active == 1)
	{
		if (random < kr_reaction_2)
		{
			xmemory->state = 41;
			reaction = 2;
			new_agent_state = 41;
			new_agent_Dc_state = Dc_state_41;
			ubr_reaction = ubr_reaction_2;
			int active = 1;
			xmemory->active = 0;
			xmemory->Dc = Dc_state_41;
		}
	}

	/* Dissociation of aBak4 into aBak2 and aBak2 */
	if (xmemory->state == 44 && xmemory->active == 1)
	{
		if (random < kr_reaction_3)
		{
			xmemory->state = 42;
			reaction = 3;
			new_agent_state = 42;
			new_agent_Dc_state = Dc_state_42;
			ubr_reaction = ubr_reaction_3;
			int active = 1;
			xmemory->active = 0;
			xmemory->Dc = Dc_state_42;
		}
	}

	/* Dissociation of aBak6 into aBak2 and aBak4 */
	if (xmemory->state == 46 && xmemory->active == 1)
	{
		if (random < kr_reaction_4)
		{
			xmemory->state = 42;
			reaction = 4;
			new_agent_state = 44;
			new_agent_Dc_state = Dc_state_44;
			ubr_reaction = ubr_reaction_4;
			int active = 1;
			xmemory->active = 0;
			xmemory->Dc = Dc_state_42;
		}
	}

	/* Dissociation of tBidBak into tBidMito and Bak */
	if (xmemory->state == 51 && xmemory->active == 1)
	{
		if (random < kr_reaction_5)
		{
			xmemory->state = 56;
			reaction = 5;
			new_agent_state = 4;
			new_agent_Dc_state = Dc_state_4;
			ubr_reaction = ubr_reaction_5;
			int active = 1;
			xmemory->active = 0;
			xmemory->Dc = Dc_state_56;
		}
	}

	/* Dissociation of tBidMcl1 into tBidMito and Mcl1Mito */
	if (xmemory->state == 54 && xmemory->active == 1)
	{
		if (random < kr_reaction_6)
		{
			xmemory->state = 56;
			reaction = 6;
			new_agent_state = 20;
			new_agent_Dc_state = Dc_state_20;
			ubr_reaction = ubr_reaction_6;
			int active = 1;
			xmemory->active = 0;
			xmemory->Dc = Dc_state_56;
		}
	}

	/* Dissociation of tBidaBak into tBidMito and aBak */
	if (xmemory->state == 52 && xmemory->active == 1)
	{
		if (random < kr_reaction_7)
		{
			xmemory->state = 56;
			reaction = 7;
			new_agent_state = 41;
			new_agent_Dc_state = Dc_state_41;
			ubr_reaction = ubr_reaction_7;
			int active = 1;
			xmemory->active = 0;
			xmemory->Dc = Dc_state_56;
		}
	}

	/* Dissociation of Mcl1aBak into Mcl1Mito and aBak */
	if (xmemory->state == 21 && xmemory->active == 1)
	{
		if (random < kr_reaction_8)
		{
			xmemory->state = 20;
			reaction = 8;
			new_agent_state = 41;
			new_agent_Dc_state = Dc_state_41;
			ubr_reaction = ubr_reaction_8;
			int active = 1;
			xmemory->active = 0;
			xmemory->Dc = Dc_state_20;
		}
	}

	/* Catalysis of tBidBak */
	if (xmemory->state == 51 && xmemory->active == 1)
	{
		if (random < kf_reaction_9)
		{
			xmemory->state = 52;
			int active = 1;
			xmemory->active = 0;
			xmemory->Dc = Dc_state_51;
		}
	}

	/* Catalysis of tBidaBak */
	if (xmemory->state == 52 && xmemory->active == 1)
	{
		if (random < kf_reaction_10)
		{
			xmemory->state = 51;
			int active = 1;
			xmemory->active = 0;
			xmemory->Dc = Dc_state_52;
		}
	}

	/* Catalysis of Mcl1aBak */
	if (xmemory->state == 21 && xmemory->active == 1)
	{
		if (random < kf_reaction_11)
		{
			xmemory->state = 20;
			reaction = 11;
			new_agent_state = 4;
			new_agent_Dc_state = Dc_state_4;
			ubr_reaction = kf_reaction_11;
			int active = 1;
			xmemory->active = 0;
			xmemory->Dc = Dc_state_20;
		}
	}

	if (reaction != 0)
	{
		int active = 1;

		/* update coordinates of dissociated new particle */
		if (new_agent_foundedge == 1 || new_agent_foundedge == 3)
		{
			dx = ubr_reaction * sin(theta);
			dy = 0;
			dz = ubr_reaction * cos(theta);
		}
		if (new_agent_foundedge == 2 || new_agent_foundedge == 4)
		{
			dx = 0;
			dy = ubr_reaction * sin(theta);
			dz = ubr_reaction * cos(theta);
		}
		if (new_agent_foundedge == 5 || new_agent_foundedge == 6)
		{
			dx = ubr_reaction * sin(theta);
			dy = ubr_reaction * cos(theta);
			dz = 0;
		}

		/* previous locations of Receptor */
	    double x_old, y_old, z_old;
	    x_old = x;
	    y_old = y;
	    z_old = z;

		x += dx;
		y += dy;
		z += dz;

		/* calculation of the walk distances until its crossing a border/plane */
		double fractionx1 = (x_old - x1)/(x_old - x);
		double fractionx2 = (x_old - x2)/(x_old - x);
		double fractiony1 = (y_old - y1)/(y_old - y);
		double fractiony3 = (y_old - y3)/(y_old - y);
		double fractionz1 = (z_old - z1)/(z_old - z);
		double fractionz5 = (z_old - z5)/(z_old - z);

		/* calculation of the walking distances into positive numbers */
		double dx_pos = dx;
		double dy_pos = dy;
		double dz_pos = dz;
		if (dx < 0)	{	dx_pos = fabs(dx);	}
		if (dy < 0)	{	dy_pos = fabs(dy);	}
		if (dz < 0)	{	dz_pos = fabs(dz);	}

		/* counts how many borders of the mito are crossed */
		int borders_crossed = 0;
		int check = 0;

		if ((x >= x1) && (x <= x2) && (y >= y1) && (y <= y3) && (z >= z1) && (x <= z5))
		{ check = 1; }

		if (check == 0)
		{
			/* Check on which plane is the ligand -> then checking if its walking outside of the plane -> then checking into which plane its walking and updating the position */
			// Mcl-1 on plane 1 or 3
			/* Mcl-1 on plane 1 or 3		->		Mcl-1 stays in the borders of plane 1 or 3 	*/
			if ((borders_crossed == 0 && (new_agent_foundedge == 1)) || (borders_crossed == 0 && (new_agent_foundedge == 3)))
			{
				if ((x >= x1) && (x <= x2) && (z >= z1) && (z <= z5))
				{
					if (new_agent_foundedge == 1)		{	y = y1;		}
					if (new_agent_foundedge == 3)		{	y = y3;		}
					borders_crossed++;
				}
				else
				{
					/* Mcl-1 on plane 1 or 3		->		walks on to plane 4	*/
					if (borders_crossed == 0 && (x < x1) && (z >= z1) && (z <= z5))
					{
						if (new_agent_foundedge == 1)		{	y = y1 + fractionx1 * dx_pos;		}
						if (new_agent_foundedge == 3)		{	y = y3 - fractionx1 * dx_pos;		}
						new_agent_foundedge = 4;
						x = x1;
						borders_crossed++;
					}
					/* Mcl-1 on plane 1 or 3		->		walks on to plane 2	*/
					if (borders_crossed == 0 && (x > x2) && (z >= z1) && (z <= z5))
					{
						if (new_agent_foundedge == 1)		{	y = y1 + fractionx2 * dx_pos;		}
						if (new_agent_foundedge == 3)		{	y = y3 - fractionx2 * dx_pos;		}
						new_agent_foundedge = 2;
						x = x2;
						borders_crossed++;
					}
					/* Mcl-1 on plane 1 or 3		->		walks on to plane 5	*/
					if (borders_crossed == 0 && (z < z1) && (x >= x1) && (x <= x2))
					{
						if (new_agent_foundedge == 1)		 {	y = y1 + fractionz1 * dz_pos;		}
						if (new_agent_foundedge == 3)		 {	y = y3 - fractionz1 * dz_pos;		}
						new_agent_foundedge = 5;
						z = z1;
						borders_crossed++;
					}
					/* Mcl-1 on plane 1 or 3		->		walks on to plane 6	*/
					if (borders_crossed == 0 && (z > z5) && (x >= x1) && (x <= x2))
					{
						if (new_agent_foundedge == 1)		{	y = y1 + fractionz5 * dz_pos;		}
						if (new_agent_foundedge == 3)		{	y = y3 - fractionz5 * dz_pos;		}
						new_agent_foundedge = 6;
						z = z5;
						borders_crossed++;
					}
					/* Mcl-1 on plane 1 or 3		->		walks on to plane 2/4/5/6	*/
					if (borders_crossed == 0 && (x < x1) && (z < z1))
					{
						if (fractionx1 < fractionz1)
						{
							if (new_agent_foundedge == 1)	{	y = y1 + fractionz1 * dz_pos;	}
							if (new_agent_foundedge == 3)	{	y = y3 - fractionz1 * dz_pos;	}
							x = x1 + fractionx1 * dx_pos;
							z = z1;
							new_agent_foundedge = 5;
							borders_crossed++;
						}
						if (fractionx1 >= fractionz1)
						{
							if (new_agent_foundedge == 1)	{	y = y1 + fractionx1 * dx_pos;	}
							if (new_agent_foundedge == 3)	{	y = y3 - fractionx1 * dx_pos;	}
							x = x1;
							z = z1 + fractionz1 * dz_pos;
							new_agent_foundedge = 4;
							borders_crossed++;
						}
					}
					if (borders_crossed == 0 && (x < x1) && (z > z5))
					{
						if (fractionx1 < fractionz5)
						{
							if (new_agent_foundedge == 1)	{	y = y1 + fractionz5 * dz_pos;	}
							if (new_agent_foundedge == 3)	{	y = y3 - fractionz5 * dz_pos;	}
							x = x1 + fractionx1 * dx_pos;
							z = z5;
							new_agent_foundedge = 6;
							borders_crossed++;
						}
						if (fractionx1 >= fractionz5)
						{
							if (new_agent_foundedge == 1)	{	y = y1 + fractionx1 * dx_pos;	}
							if (new_agent_foundedge == 3)	{	y = y3 - fractionx1 * dx_pos;	}
							x = x1;
							z = z5 - fractionz5 * dz_pos;
							new_agent_foundedge = 4;
							borders_crossed++;
						}
					}
					if (borders_crossed == 0 && (x > x2) && (z < z1))
					{
						if (fractionx2 < fractionz1)
						{
							if (new_agent_foundedge == 1)	{	y = y1 + fractionz1 * dz_pos;	}
							if (new_agent_foundedge == 3)	{	y = y3 - fractionz1 * dz_pos;	}
							x = x2 - fractionx2 * dx_pos;
							z = z1;
							new_agent_foundedge = 5;
							borders_crossed++;
						}
						if (fractionx2 >= fractionz1)
						{
							if (new_agent_foundedge == 1)	{	y = y1 + fractionx2 * dx_pos;	}
							if (new_agent_foundedge == 3)	{	y = y3 - fractionx2 * dx_pos;	}
							x = x2;
							z = z1 + fractionz1 * dz_pos;
							new_agent_foundedge = 2;
							borders_crossed++;
						}
					}
					if (borders_crossed == 0 && (x > x2) && (z > z5))
					{
						if (fractionx2 < fractionz5)
						{
							if (new_agent_foundedge == 1)	{	y = y1 + fractionz5 * dz_pos;	}
							if (new_agent_foundedge == 3)	{	y = y3 - fractionz5 * dz_pos;	}
							x = x2 - fractionx2 * dx_pos;
							z = z5;
							new_agent_foundedge = 6;
							borders_crossed++;
						}
						if (fractionx2 >= fractionz5)
						{
							if (new_agent_foundedge == 1)	{	y = y1 + fractionx2 * dx_pos;	}
							if (new_agent_foundedge == 3)	{	y = y3 - fractionx2 * dx_pos;	}
							x = x2;
							z = z5 - fractionz5 * dz_pos;
							new_agent_foundedge = 2;
							borders_crossed++;
						}
					}
				}
			}

			//	Mcl-1 on plane 2 or 4
			/* Mcl-1 on plane 2 or 4		->		Mcl-1 stays in the borders of plane 2 or 4 	*/
			if ((borders_crossed == 0 && (new_agent_foundedge == 2)) || (borders_crossed == 0 && (new_agent_foundedge == 4)))
			{
				if ((y >= y1) && (y <= y3) && (z >= z1) && (z <= z5))
				{
					if (new_agent_foundedge == 2)		{	x = x2;		}
					if (new_agent_foundedge == 4)		{	x = x1;		}
					borders_crossed++;
				}
				else
				{
					/* Mcl-1 on plane 2 or 4		->		walks on to plane 1	*/
					if (borders_crossed == 0 && (y < y1) && (z >= z1) && (z <= z5))
					{
						if (new_agent_foundedge == 2)		{	x = x2 - fractiony1 * dy_pos;		}
						if (new_agent_foundedge == 4)		{	x = x1 + fractiony1 * dy_pos;		}
						new_agent_foundedge = 1;
						y = y1;
						borders_crossed++;
					}
					/* Mcl-1 on plane 2 or 4		->		walks on to plane 3	*/
					if (borders_crossed == 0 && (y > y3) && (z >= z1) && (z <= z5))
					{
						if (new_agent_foundedge == 2)		{	x = x2 - fractiony3 * dy_pos;		}
						if (new_agent_foundedge == 4)		{	x = x1 + fractiony3 * dy_pos;		}
						new_agent_foundedge = 3;
						y = y3;
						borders_crossed++;
					}
					/* Mcl-1 on plane 2 or 4		->		walks on to plane 5	*/
					if (borders_crossed == 0 && (z < z1) && (y >= y1) && (y <= y3))
					{
						if (new_agent_foundedge == 2)		{	x = x2 - fractionz1 * dz_pos;		}
						if (new_agent_foundedge == 4)		{	x = x1 + fractionz1 * dz_pos;		}
						new_agent_foundedge = 5;
						z = z1;
						borders_crossed++;
					}
					/* Mcl-1 on plane 2 or 4		->		walks on to plane 6	*/
					if (borders_crossed == 0 && (z > z5) && (y >= y1) && (y <= y3))
					{
						if (new_agent_foundedge == 2)		{	x = x2 - fractionz5 * dz_pos;		}
						if (new_agent_foundedge == 4)		{	x = x1 + fractionz5 * dz_pos;		}
						new_agent_foundedge = 6;
						z = z5;
						borders_crossed++;
					}

					/* Mcl-1 on plane 2 or 4		->		walks on to plane 1/3/5/6	*/
					if (borders_crossed == 0 && (y < y1) && (z < z1))
					{
						if (fractiony1 < fractionz1)
						{
							if (new_agent_foundedge == 2)	{	x = x1 + fractionz1 * dz_pos;	}
							if (new_agent_foundedge == 4)	{	x = x2 - fractionz1 * dz_pos;	}
							y = y1 + fractiony1 * dy_pos;
							z = z1;
							new_agent_foundedge = 5;
							borders_crossed++;
						}
						if (fractiony1 >= fractionz1)
						{
							if (new_agent_foundedge == 2)	{	x = x1 + fractiony1 * dy_pos;	}
							if (new_agent_foundedge == 4)	{	x = x2 - fractiony1 * dy_pos;	}
							y = y1;
							z = z1 + fractionz1 * dz_pos;
							new_agent_foundedge = 1;
							borders_crossed++;
						}
					}
					if (borders_crossed == 0 && (y < y1) && (z > z5))
					{
						if (fractiony1 < fractionz5)
						{
							if (new_agent_foundedge == 2)	{	x = x1 + fractionz5 * dz_pos;	}
							if (new_agent_foundedge == 4)	{	x = x2 - fractionz5 * dz_pos;	}
							y = y1 + fractiony1 * dy_pos;
							z = z5;
							new_agent_foundedge = 6;
							borders_crossed++;
						}
						if (fractiony1 >= fractionz5)
						{
							if (new_agent_foundedge == 2)	{	x = x1 + fractiony1 * dy_pos;	}
							if (new_agent_foundedge == 4)	{	x = x2 - fractiony1 * dy_pos;	}
							y = y1;
							z = z5 - fractionz5 * dz_pos;
							new_agent_foundedge = 1;
							borders_crossed++;
						}
					}
					if (borders_crossed == 0 && (y > y3) && (z < z1))
					{
						if (fractiony3 < fractionz1)
						{
							if (new_agent_foundedge == 2)	{	x = x1 + fractionz1 * dz_pos;	}
							if (new_agent_foundedge == 4)	{	x = x2 - fractionz1 * dz_pos;	}
							y = y3 - fractiony3 * dy_pos;
							z = z1;
							new_agent_foundedge = 5;
							borders_crossed++;
						}
						if (fractiony3 >= fractionz1)
						{
							if (new_agent_foundedge == 2)	{	x = x1 + fractiony3 * dy_pos;	}
							if (new_agent_foundedge == 4)	{	x = x2 - fractiony3 * dy_pos;	}
							y = y3;
							z = z1 + fractionz1 * dz_pos;
							new_agent_foundedge = 3;
							borders_crossed++;
						}
					}
					if (borders_crossed == 0 && (y > y3) && (z > z5))
					{
						if (fractiony3 < fractionz5)
						{
							if (new_agent_foundedge == 2)	{	x = x1 + fractionz5 * dz_pos;	}
							if (new_agent_foundedge == 4)	{	x = x2 - fractionz5 * dz_pos;	}
							y = y3 - fractiony3 * dy_pos;
							z = z5;
							new_agent_foundedge = 6;
						}
						if (fractiony3 >= fractionz5)
						{
							if (new_agent_foundedge == 2)	{	x = x1 + fractiony3 * dy_pos;	}
							if (new_agent_foundedge == 4)	{	x = x2 - fractiony3 * dy_pos;	}
							y = y3;
							z = z5 - fractionz5 * dz_pos;
							new_agent_foundedge = 3;
							borders_crossed++;
						}
					}
				}
			}

			// Mcl-1 on plane 5 or 6
			/* Mcl-1 on plane 5 or 6		->		Mcl-1 stays in the borders of plane 5 or 5 	*/
			if ((borders_crossed == 0 && (new_agent_foundedge == 5)) || (borders_crossed == 0 && (new_agent_foundedge == 6)))
			{
				if ((x >= x1) && (x <= x2) && (y >= y1) && (y <= y3))
				{
					if (new_agent_foundedge == 5)		{	z = z1;		}
					if (new_agent_foundedge == 6)		{	z = z5;		}
					borders_crossed++;
				}

				else
				{
					/* Mcl-1 on plane 5 or 6		->		walks on to plane 4	*/
					if (borders_crossed == 0 && (x < x1) && (y >= y1) && (y <= y3))
					{
						if (new_agent_foundedge == 5)		{	z = z1 + fractionx1 * dx_pos;		}
						if (new_agent_foundedge == 6)		{	z = z5 - fractionx1 * dx_pos;		}
						new_agent_foundedge = 4;
						x = x1;
						borders_crossed++;
					}
					/* Mcl-1 on plane 5 or 6		->		walks on to plane 2	*/
					if (borders_crossed == 0 && (x > x2) && (y >= y1) && (y <= y3))
					{
						if (new_agent_foundedge == 5)		{	z = z1 + fractionx2 * dx_pos;		}
						if (new_agent_foundedge == 6)		{	z = z5 - fractionx2 * dx_pos;		}
						new_agent_foundedge = 2;
						x = x2;
						borders_crossed++;
					}
					/* Mcl-1 on plane 5 or 6		->		walks on to plane 1	*/
					if (borders_crossed == 0 && (y < y1) && (x >= x1) && (x <= x2))
					{
						if (new_agent_foundedge == 5)		{	z = z1 + fractiony1 * dy_pos;		}
						if (new_agent_foundedge == 6)		{	z = z5 - fractiony1 * dy_pos;		}
						new_agent_foundedge = 1;
						y = y1;
						borders_crossed++;
					}
					/* Mcl-1 on plane 5 or 6		->		walks on to plane 3	*/
					if (borders_crossed == 0 && (y > y3) && (x >= x1) && (x <= x2))
					{
						if (new_agent_foundedge == 5)		{	z = z1 + fractiony3 * dy_pos;		}
						if (new_agent_foundedge == 6)		{	z = z5 - fractiony3 * dy_pos;		}
						new_agent_foundedge = 3;
						y = y3;
						borders_crossed++;
					}

					/* Mcl-1 on plane 5 or 6		->		walks on to plane 1/2/3/4	*/
					if (borders_crossed == 0 && (x < x1) && (y < y1))
					{
						if (fractionx1 < fractiony1)
						{
							if (new_agent_foundedge == 5)	{	z = z1 + fractiony1 * dy_pos;	}
							if (new_agent_foundedge == 6)	{	z = z5 - fractiony1 * dy_pos;	}
							x = x1 + fractionx1 * dx_pos;
							y = y1;
							new_agent_foundedge = 1;
							borders_crossed++;
						}
						if (fractionx1 >= fractiony1)
						{
							if (new_agent_foundedge == 5)	{	z = z1 + fractionx1 * dx_pos;	}
							if (new_agent_foundedge == 6)	{	z = z5 - fractionx1 * dx_pos;	}
							x = x1;
							y = y1 + fractiony1 * dy_pos;
							new_agent_foundedge = 4;
							borders_crossed++;
						}
					}
					if (borders_crossed == 0 && (x < x1) && (y > y3))
					{
						if (fractionx1 < fractiony3)
						{
							if (new_agent_foundedge == 5)	{	z = z1 + fractiony3 * dy_pos;	}
							if (new_agent_foundedge == 6)	{	z = z5 - fractiony3 * dy_pos;	}
							x = x1 + fractionx1 * dx_pos;
							y = y3;
							new_agent_foundedge = 3;
							borders_crossed++;
						}
						if (fractionx1 >= fractiony3)
						{
							if (new_agent_foundedge == 5)	{	z = z1 + fractionx1 * dx_pos;	}
							if (new_agent_foundedge == 6)	{	z = z5 - fractionx1 * dx_pos;	}
							x = x1;
							y = y3 - fractiony3 * dy_pos;
							new_agent_foundedge = 4;
							borders_crossed++;
						}
					}
					if (borders_crossed == 0 && (x > x2) && (y < y1))
					{
						if (fractionx2 < fractiony1)
						{
							if (new_agent_foundedge == 5)	{	z = z1 + fractiony1 * dy_pos;	}
							if (new_agent_foundedge == 6)	{	z = z5 - fractiony1 * dy_pos;	}
							x = x2 - fractionx2 * dx_pos;
							y = y1;
							new_agent_foundedge = 1;
							borders_crossed++;
						}
						if (fractionx2 >= fractiony1)
						{
							if (new_agent_foundedge == 5)	{	z = z1 + fractionx2 * dx_pos;	}
							if (new_agent_foundedge == 6)	{	z = z5 - fractionx2 * dx_pos;	}
							x = x2;
							y = y1 + fractiony1 * dy_pos;
							new_agent_foundedge = 2;
							borders_crossed++;
						}
					}
					if (borders_crossed == 0 && (x > x2) && (y > y3))
					{
						if (fractionx2 < fractiony3)
						{
							if (new_agent_foundedge == 5)	{	z = z1 + fractiony3 * dy_pos;	}
							if (new_agent_foundedge == 6)	{	z = z5 - fractiony3 * dy_pos;	}
							x = x2 - fractionx2 * dx_pos;
							y = y3;
							new_agent_foundedge = 3;
							borders_crossed++;
						}
						if (fractionx2 >= fractiony3)
						{
							if (new_agent_foundedge == 5)	{	z = z1 + fractionx2 * dx_pos;	}
							if (new_agent_foundedge == 6)	{	z = z5 - fractionx2 * dx_pos;	}
							x = x2;
							y = y3 - fractiony3 * dy_pos;
							new_agent_foundedge = 2;
							borders_crossed++;
						}
					}
				}
			}

			if ((x >= x1) && (x <= x2) && (y >= y1) && (y <= y3) && (z >= z1) && (x <= z5))
			{ check = 1; }
		}

		// create new receptor agent after reaction (dissociation)
		add_Receptor_agent (Receptor_agents, id_new, new_agent_state, -1, active, x, y, z, new_agent_Dc_state, xmemory->collision_counter_ligand, xmemory->check_collision, new_agent_foundedge, xmemory->mitoid, xmemory->mitosize);
	}

return 0;
}


/***************************/
/** Receptor Bind Receptor */
/***************************/

__FLAME_GPU_FUNC__ int Receptor_Bind_Receptor(xmachine_memory_Receptor* xmemory, xmachine_memory_Receptor_list* Receptor_agents, xmachine_message_Receptor_Location_list* Receptor_Location_messages, xmachine_message_Receptor_Location_PBM* partition_matrix, xmachine_message_Receptor_Bound_list* Receptor_Bound_messages, RNG_rand48* rand48)
{
    double closestdistance = 99999.0;
    int Bound = -1;
    int Boundstate = 0;
    int Reaction = 0;
    double distance = 0.0;
        xmachine_message_Receptor_Location* Receptor_Location_message;
        Receptor_Location_message = get_first_Receptor_Location_message(Receptor_Location_messages, partition_matrix, xmemory->x, xmemory->y, xmemory->z);
        while (Receptor_Location_message)
        {
    	if (xmemory->active == 1) {
    if (xmemory->state == 20)
    {
            if (Receptor_Location_message->state == 56 && Receptor_Location_message->id != xmemory->id && Receptor_Location_message->mitoid == xmemory->mitoid)
            {
                distance = sqrt((xmemory->x - Receptor_Location_message->x)*(xmemory->x - Receptor_Location_message->x)+(xmemory->y - Receptor_Location_message->y)*(xmemory->y - Receptor_Location_message->y)+(xmemory->z - Receptor_Location_message->z)*(xmemory->z - Receptor_Location_message->z));
                if (distance < br_reaction_6 && distance < closestdistance)
                {
                    closestdistance = distance;
                    Bound = Receptor_Location_message->id;
                    Boundstate = 56;
                    Reaction = 6;
                }
            }
            if (Receptor_Location_message->state == 41 && Receptor_Location_message->id != xmemory->id && Receptor_Location_message->mitoid == xmemory->mitoid)
            {
                distance = sqrt((xmemory->x - Receptor_Location_message->x)*(xmemory->x - Receptor_Location_message->x)+(xmemory->y - Receptor_Location_message->y)*(xmemory->y - Receptor_Location_message->y)+(xmemory->z - Receptor_Location_message->z)*(xmemory->z - Receptor_Location_message->z));
                if (distance < br_reaction_8 && distance < closestdistance)
                {
                    closestdistance = distance;
                    Bound = Receptor_Location_message->id;
                    Boundstate = 41;
                    Reaction = 8;
                }
            }
    }

    if (xmemory->state == 4)
    {
            if (Receptor_Location_message->state == 41 && Receptor_Location_message->id != xmemory->id && Receptor_Location_message->mitoid == xmemory->mitoid)
            {
                distance = sqrt((xmemory->x - Receptor_Location_message->x)*(xmemory->x - Receptor_Location_message->x)+(xmemory->y - Receptor_Location_message->y)*(xmemory->y - Receptor_Location_message->y)+(xmemory->z - Receptor_Location_message->z)*(xmemory->z - Receptor_Location_message->z));
                if (distance < br_reaction_1 && distance < closestdistance)
                {
                    closestdistance = distance;
                    Bound = Receptor_Location_message->id;
                    Boundstate = 41;
                    Reaction = 1;
                }
            }
            if (Receptor_Location_message->state == 56 && Receptor_Location_message->id != xmemory->id && Receptor_Location_message->mitoid == xmemory->mitoid)
            {
                distance = sqrt((xmemory->x - Receptor_Location_message->x)*(xmemory->x - Receptor_Location_message->x)+(xmemory->y - Receptor_Location_message->y)*(xmemory->y - Receptor_Location_message->y)+(xmemory->z - Receptor_Location_message->z)*(xmemory->z - Receptor_Location_message->z));
                if (distance < br_reaction_5 && distance < closestdistance)
                {
                    closestdistance = distance;
                    Bound = Receptor_Location_message->id;
                    Boundstate = 56;
                    Reaction = 5;
                }
            }
    }

    if (xmemory->state == 41)
    {
            if (Receptor_Location_message->state == 4 && Receptor_Location_message->id != xmemory->id && Receptor_Location_message->mitoid == xmemory->mitoid)
            {
                distance = sqrt((xmemory->x - Receptor_Location_message->x)*(xmemory->x - Receptor_Location_message->x)+(xmemory->y - Receptor_Location_message->y)*(xmemory->y - Receptor_Location_message->y)+(xmemory->z - Receptor_Location_message->z)*(xmemory->z - Receptor_Location_message->z));
                if (distance < br_reaction_1 && distance < closestdistance)
                {
                    closestdistance = distance;
                    Bound = Receptor_Location_message->id;
                    Boundstate = 4;
                    Reaction = 1;
                }
            }
            if (Receptor_Location_message->state == 41 && Receptor_Location_message->id != xmemory->id && Receptor_Location_message->mitoid == xmemory->mitoid)
            {
                distance = sqrt((xmemory->x - Receptor_Location_message->x)*(xmemory->x - Receptor_Location_message->x)+(xmemory->y - Receptor_Location_message->y)*(xmemory->y - Receptor_Location_message->y)+(xmemory->z - Receptor_Location_message->z)*(xmemory->z - Receptor_Location_message->z));
                if (distance < br_reaction_2 && distance < closestdistance)
                {
                    closestdistance = distance;
                    Bound = Receptor_Location_message->id;
                    Boundstate = 41;
                    Reaction = 2;
                }
            }
            if (Receptor_Location_message->state == 56 && Receptor_Location_message->id != xmemory->id && Receptor_Location_message->mitoid == xmemory->mitoid)
            {
                distance = sqrt((xmemory->x - Receptor_Location_message->x)*(xmemory->x - Receptor_Location_message->x)+(xmemory->y - Receptor_Location_message->y)*(xmemory->y - Receptor_Location_message->y)+(xmemory->z - Receptor_Location_message->z)*(xmemory->z - Receptor_Location_message->z));
                if (distance < br_reaction_7 && distance < closestdistance)
                {
                    closestdistance = distance;
                    Bound = Receptor_Location_message->id;
                    Boundstate = 56;
                    Reaction = 7;
                }
            }
            if (Receptor_Location_message->state == 20 && Receptor_Location_message->id != xmemory->id && Receptor_Location_message->mitoid == xmemory->mitoid)
            {
                distance = sqrt((xmemory->x - Receptor_Location_message->x)*(xmemory->x - Receptor_Location_message->x)+(xmemory->y - Receptor_Location_message->y)*(xmemory->y - Receptor_Location_message->y)+(xmemory->z - Receptor_Location_message->z)*(xmemory->z - Receptor_Location_message->z));
                if (distance < br_reaction_8 && distance < closestdistance)
                {
                    closestdistance = distance;
                    Bound = Receptor_Location_message->id;
                    Boundstate = 20;
                    Reaction = 8;
                }
            }
    }

    if (xmemory->state == 42)
    {
            if (Receptor_Location_message->state == 42 && Receptor_Location_message->id != xmemory->id && Receptor_Location_message->mitoid == xmemory->mitoid)
            {
                distance = sqrt((xmemory->x - Receptor_Location_message->x)*(xmemory->x - Receptor_Location_message->x)+(xmemory->y - Receptor_Location_message->y)*(xmemory->y - Receptor_Location_message->y)+(xmemory->z - Receptor_Location_message->z)*(xmemory->z - Receptor_Location_message->z));
                if (distance < br_reaction_3 && distance < closestdistance)
                {
                    closestdistance = distance;
                    Bound = Receptor_Location_message->id;
                    Boundstate = 42;
                    Reaction = 3;
                }
            }
            if (Receptor_Location_message->state == 44 && Receptor_Location_message->id != xmemory->id && Receptor_Location_message->mitoid == xmemory->mitoid)
            {
                distance = sqrt((xmemory->x - Receptor_Location_message->x)*(xmemory->x - Receptor_Location_message->x)+(xmemory->y - Receptor_Location_message->y)*(xmemory->y - Receptor_Location_message->y)+(xmemory->z - Receptor_Location_message->z)*(xmemory->z - Receptor_Location_message->z));
                if (distance < br_reaction_4 && distance < closestdistance)
                {
                    closestdistance = distance;
                    Bound = Receptor_Location_message->id;
                    Boundstate = 44;
                    Reaction = 4;
                }
            }
    }

    if (xmemory->state == 44)
    {
            if (Receptor_Location_message->state == 42 && Receptor_Location_message->id != xmemory->id && Receptor_Location_message->mitoid == xmemory->mitoid)
            {
                distance = sqrt((xmemory->x - Receptor_Location_message->x)*(xmemory->x - Receptor_Location_message->x)+(xmemory->y - Receptor_Location_message->y)*(xmemory->y - Receptor_Location_message->y)+(xmemory->z - Receptor_Location_message->z)*(xmemory->z - Receptor_Location_message->z));
                if (distance < br_reaction_4 && distance < closestdistance)
                {
                    closestdistance = distance;
                    Bound = Receptor_Location_message->id;
                    Boundstate = 42;
                    Reaction = 4;
                }
            }
    }

    if (xmemory->state == 56)
    {
            if (Receptor_Location_message->state == 4 && Receptor_Location_message->id != xmemory->id && Receptor_Location_message->mitoid == xmemory->mitoid)
            {
                distance = sqrt((xmemory->x - Receptor_Location_message->x)*(xmemory->x - Receptor_Location_message->x)+(xmemory->y - Receptor_Location_message->y)*(xmemory->y - Receptor_Location_message->y)+(xmemory->z - Receptor_Location_message->z)*(xmemory->z - Receptor_Location_message->z));
                if (distance < br_reaction_5 && distance < closestdistance)
                {
                    closestdistance = distance;
                    Bound = Receptor_Location_message->id;
                    Boundstate = 4;
                    Reaction = 5;
                }
            }
            if (Receptor_Location_message->state == 20 && Receptor_Location_message->id != xmemory->id && Receptor_Location_message->mitoid == xmemory->mitoid)
            {
                distance = sqrt((xmemory->x - Receptor_Location_message->x)*(xmemory->x - Receptor_Location_message->x)+(xmemory->y - Receptor_Location_message->y)*(xmemory->y - Receptor_Location_message->y)+(xmemory->z - Receptor_Location_message->z)*(xmemory->z - Receptor_Location_message->z));
                if (distance < br_reaction_6 && distance < closestdistance)
                {
                    closestdistance = distance;
                    Bound = Receptor_Location_message->id;
                    Boundstate = 20;
                    Reaction = 6;
                }
            }
            if (Receptor_Location_message->state == 41 && Receptor_Location_message->id != xmemory->id && Receptor_Location_message->mitoid == xmemory->mitoid)
            {
                distance = sqrt((xmemory->x - Receptor_Location_message->x)*(xmemory->x - Receptor_Location_message->x)+(xmemory->y - Receptor_Location_message->y)*(xmemory->y - Receptor_Location_message->y)+(xmemory->z - Receptor_Location_message->z)*(xmemory->z - Receptor_Location_message->z));
                if (distance < br_reaction_7 && distance < closestdistance)
                {
                    closestdistance = distance;
                    Bound = Receptor_Location_message->id;
                    Boundstate = 41;
                    Reaction = 7;
                }
            }
    }
    }

            Receptor_Location_message = get_next_Receptor_Location_message(Receptor_Location_message, Receptor_Location_messages, partition_matrix);
        }
    if (Bound != -1)
    {
        xmemory->Pending = Bound;
        add_Receptor_Bound_message (Receptor_Bound_messages, Bound, xmemory->id, closestdistance, Reaction);
    }

return 0;
}


/***************************/
/*** Receptor Check Bound **/
/***************************/

__FLAME_GPU_FUNC__ int Receptor_Check_Bound(xmachine_memory_Receptor* xmemory, xmachine_memory_Receptor_list* Receptor_agents, xmachine_message_Receptor_Bound_list* Receptor_Bound_messages, xmachine_message_MReceptor_Confirm_list* MReceptor_Confirm_messages, RNG_rand48* rand48)
{
    int closest_bound = -1;
    double closest_distance = 99999.0;
    int Reaction = 0;
	int die = 0;
    xmachine_message_Receptor_Bound* Receptor_Bound_message;
    Receptor_Bound_message = get_first_Receptor_Bound_message(Receptor_Bound_messages);
    while (Receptor_Bound_message)
    {
        if (xmemory->Pending != -1)
        {
            if (Receptor_Bound_message->id == xmemory->id)
            {
                if (Receptor_Bound_message->idfrom == xmemory->Pending)
                {
	                if (Receptor_Bound_message->Reaction == 1)
	                {
	                    if (xmemory->state == 4)
	                    {
	                        xmemory->state = 42;
	                        xmemory->Pending = -1;
	                    }
	                    if (xmemory->state == 41)
	                    {
	                        die = 1;
	                    }
	                }
	                if (Receptor_Bound_message->Reaction == 2)
	                {
	                    if (xmemory->state == 41 && Receptor_Bound_message->idfrom < xmemory->id)
	                    {
	                        xmemory->state = 42;
	                        xmemory->Pending = -1;
	                    }
	                    if (xmemory->state == 41 && Receptor_Bound_message->idfrom > xmemory->id)
	                    {
	                        die = 1;
	                    }
	                }
	                if (Receptor_Bound_message->Reaction == 3)
	                {
	                    if (xmemory->state == 42 && Receptor_Bound_message->idfrom < xmemory->id)
	                    {
	                        xmemory->state = 44;
	                        xmemory->Pending = -1;
	                    }
	                    if (xmemory->state == 42 && Receptor_Bound_message->idfrom > xmemory->id)
	                    {
	                        die = 1;
	                    }
	                }
	                if (Receptor_Bound_message->Reaction == 4)
	                {
	                    if (xmemory->state == 42)
	                    {
	                        xmemory->state = 46;
	                        xmemory->Pending = -1;
	                    }
	                    if (xmemory->state == 44)
	                    {
	                        die = 1;
	                    }
	                }
	                if (Receptor_Bound_message->Reaction == 5)
	                {
	                    if (xmemory->state == 56)
	                    {
	                        xmemory->state = 51;
	                        xmemory->Pending = -1;
	                    }
	                    if (xmemory->state == 4)
	                    {
	                        die = 1;
	                    }
	                }
	                if (Receptor_Bound_message->Reaction == 6)
	                {
	                    if (xmemory->state == 56)
	                    {
	                        xmemory->state = 54;
	                        xmemory->Pending = -1;
	                    }
	                    if (xmemory->state == 20)
	                    {
	                        die = 1;
	                    }
	                }
	                if (Receptor_Bound_message->Reaction == 7)
	                {
	                    if (xmemory->state == 56)
	                    {
	                        xmemory->state = 52;
	                        xmemory->Pending = -1;
	                    }
	                    if (xmemory->state == 41)
	                    {
	                        die = 1;
	                    }
	                }
	                if (Receptor_Bound_message->Reaction == 8)
	                {
	                    if (xmemory->state == 20)
	                    {
	                        xmemory->state = 21;
	                        xmemory->Pending = -1;
	                    }
	                    if (xmemory->state == 41)
	                    {
	                        die = 1;
	                    }
	                }
                }
				else
                {
                    if (Receptor_Bound_message->distance < closest_distance)
                    {
                        closest_distance = Receptor_Bound_message->distance;
                        closest_bound = Receptor_Bound_message->idfrom;
                        Reaction = Receptor_Bound_message->Reaction;
                    }
                }
            }
        }
        Receptor_Bound_message = get_next_Receptor_Bound_message(Receptor_Bound_message, Receptor_Bound_messages);
    }
	if (die == 1)
	{
		return 1;
	}
    
    
        if (closest_bound != -1)
        {
        	add_MReceptor_Confirm_message(MReceptor_Confirm_messages, closest_bound, Reaction);
        }

return 0;
}


/**********************************/
/***** Receptor Confirm Ligand ****/
/**********************************/


/************************************/
/***** Receptor Confirm Receptor ****/
/************************************/

__FLAME_GPU_FUNC__ int Receptor_Confirm_Receptor(xmachine_memory_Receptor* xmemory, xmachine_memory_Receptor_list* Receptor_agents, xmachine_message_MReceptor_Confirm_list* MReceptor_Confirm_messages, RNG_rand48* rand48)
{
xmemory->Pending = -1;
xmemory->active = 1;
return 0;
}
#endif //_FLAMEGPU_FUNCTIONS
