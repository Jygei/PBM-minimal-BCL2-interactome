/*

Jenny Geiger
Institute of Cell Biology and Immunology (IZI) at the University of Stuttgart

FaST generated Codes were changed and modified to simulate a whole cell environment to simulate subcellular localization of MCL-1






Code version information:
- this code is for the whole-cell particle-based model: particles diffuse and potentially bind to the mitochondrial surface where they are able to diffuse laterally
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
__FLAME_GPU_FUNC__ int Mito_CC(xmachine_memory_Mito* xmemory, /*xmachine_message_Ligand_Location_list* Ligand_Location_messages, /*xmachine_message_Mito_Location_list* Mito_Location_messages, */ RNG_rand48* rand48)
{
    return 0;
}
	
/*****************************************************************************/
/****************************** Ligand Functions *****************************/
/*****************************************************************************/

/***************************/
/****** Ligand Output ******/
/***************************/

__FLAME_GPU_FUNC__ int Ligand_Output(xmachine_memory_Ligand* xmemory, xmachine_memory_Receptor_list* Receptor_agents, xmachine_message_Mito_Location_list* Mito_Location_messages, xmachine_message_Ligand_Location_list* Ligand_Location_messages, RNG_rand48* rand48)
{
	/* variables to check MCL1 collisions witch Mitos */
	int check_bound = 0;
	double collision_counter_ligand = xmemory->collision_counter_ligand;
	xmemory->check_collision = 0;										// collision-counter per ligand per iteration (always =0 in the beginning)
	
	/* probability: MCL1 binding/reflection */
	double probability;
	probability = 0.00000310367; /* = 0.00249982/175.6623168; after mito_number calculation with median */

	/* Readout of the Mito_Location messages*/
	int mito_counter = 0;
	int mito_counter2 = 0;
	
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
	}
	
	/* Ligands diffusing and checking collision with the boundaries of the cell and afterwards of the mitos */
	int counter = 0;
	while (counter < brownianloop)
	{		
		/* previous locations of Ligand */
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
        
        double deltax = xmemory->Dc * n1;
        double deltay = xmemory->Dc * n2;
        double deltaz = xmemory->Dc * n3;

        
		/* Update x, y and z Variables */
        xmemory->x += deltax;
        xmemory->y += deltay;
        xmemory->z += deltaz;
		
        
		/* periodic boundary conditions */	 
		if ((xmemory->state == 2) || (xmemory->state == 4) || (xmemory->state == 5));
        {
            if (xmemory->x > maxintracellular_xboundary)
            {
                xmemory->x = minintracellular_xboundary + (xmemory->x - maxintracellular_xboundary);
            }
            if (xmemory->x < minintracellular_xboundary)
            {
                xmemory->x = maxintracellular_xboundary - (minintracellular_xboundary - xmemory->x);
            }
            if (xmemory->y > maxintracellular_yboundary)
            {
                xmemory->y = minintracellular_yboundary + (xmemory->y - maxintracellular_yboundary);
            }
            if (xmemory->y < minintracellular_yboundary)
            {
                xmemory->y = maxintracellular_yboundary - (minintracellular_yboundary - xmemory->y);
            }
            if (xmemory->z > maxintracellular_zboundary)
            {
                xmemory->z = minintracellular_zboundary + (xmemory->z - maxintracellular_zboundary);
            }
            if (xmemory->z < minintracellular_zboundary)
            {
                xmemory->z = maxintracellular_zboundary - (minintracellular_zboundary - xmemory->z);
            }
		}
	

		/* Check if Ligand collides with Mitos, reflect off Mito boundaries and count collisions */
		check_bound = 0;
		mito_counter = 0;		//counter which safes the mito ID with which it collides
		mito_counter2 = 0;		//counter for the general collision check with a mito
		int foundedge = 0;
		int edges_found = 0;
		double intersectx, intersecty, intersectz;
		double saveintersectx, saveintersecty, saveintersectz;
		
		while (mito_counter2 < mito_number)
		{				
			if (((xmemory->x > x1[mito_counter2]) && (xmemory->x < x2[mito_counter2])) && ((xmemory->y > y1[mito_counter2]) && (xmemory->y < y3[mito_counter2])) && ((xmemory->z > z1[mito_counter2]) && (xmemory->z < z5[mito_counter2])))
			{
				xmemory->collision_counter_ligand += 1;				
				xmemory->check_collision += 1;
				edges_found = 0;
				xmemory->mitosize = current_mitosize[mito_counter2];
				// Check Edge 1: Defined by vertices 1, 2, 5 & 6	C
				if (y_old <= y1[mito_counter2] && y1[mito_counter2] <= xmemory->y)
				{
					// Crosses the plane: Check point of intersection
					double fraction = (y1[mito_counter2]-y_old)/(xmemory->y - y_old);
					double intersectx = fraction * (xmemory->x - x_old) + x_old;
					double intersecty = y1[mito_counter2];
					double intersectz = fraction * (xmemory->z - z_old) + z_old;
					// If intersection point is on mito
					if (intersectx > x1[mito_counter2] && intersectx < x2[mito_counter2] && intersectz > z1[mito_counter2] && intersectz < z5[mito_counter2])
					{
						foundedge = 1;
						edges_found++;
						saveintersectx = intersectx;
						saveintersecty = intersecty;
						saveintersectz = intersectz;
					}
				}
				// Check Edge 2: Defined by vertices 2, 4, 6 & 8	F
				if (x_old >= x2[mito_counter2] && x2[mito_counter2] >= xmemory->x)
				{
					// Crosses the plane: Check point of intersection
					double fraction = (x2[mito_counter2]-x_old)/(xmemory->x - x_old);
					double intersectx = x2[mito_counter2];
					double intersecty = fraction * (xmemory->y - y_old) + y_old;
					double intersectz = fraction * (xmemory->z - z_old) + z_old;
					// If intersection point is on mito
					if (intersecty > y2[mito_counter2] && intersecty < y4[mito_counter2] && intersectz > z2[mito_counter2] && intersectz < z6[mito_counter2])
					{
						foundedge = 2;
						edges_found++;
						saveintersectx = intersectx;
						saveintersecty = intersecty;
						saveintersectz = intersectz;
					}
				}
				// Check Edge 3: Defined by vertices 3, 4, 7 & 8	D
				if (y_old >= y3[mito_counter2] && y3[mito_counter2] >= xmemory->y)
				{
					// Crosses the plane: Check point of intersection
					double fraction = (y3[mito_counter2]-y_old)/(xmemory->y - y_old);
					double intersectx = fraction * (xmemory->x - x_old) + x_old;
					double intersecty = y3[mito_counter2];
					double intersectz = fraction * (xmemory->z - z_old) + z_old;
					// If intersection point is on mito
					if (intersectx > x3[mito_counter2] && intersectx < x4[mito_counter2] && intersectz > z3[mito_counter2] && intersectz < z7[mito_counter2])
					{
						foundedge = 3;
						edges_found++;
						saveintersectx = intersectx;
						saveintersecty = intersecty;
						saveintersectz = intersectz;
					}
				}
				// Check Edge 4: Defined by vertices 1, 3, 5 & 7	E
				if (x_old <= x1[mito_counter2] && x1[mito_counter2] <= xmemory->x)
				{
					// Crosses the plane: Check point of intersection
					double fraction = (x1[mito_counter2]-x_old)/(xmemory->x - x_old);
					double intersectx = x1[mito_counter2];
					double intersecty = fraction * (xmemory->y - y_old) + y_old;
					double intersectz = fraction * (xmemory->z - z_old) + z_old;
					// If intersection point is on mito
					if (intersecty > y1[mito_counter2] && intersecty < y3[mito_counter2] && intersectz > z1[mito_counter2] && intersectz < z5[mito_counter2])
					{
						foundedge = 4;
						edges_found++;
						saveintersectx = intersectx;
						saveintersecty = intersecty;
						saveintersectz = intersectz;
					}
				}
				// Check Edge 5: Defined by vertices 1, 2, 3 & 4	A
				if (z_old <= z1[mito_counter2] && z1[mito_counter2] <= xmemory->z)
				{
					// Crosses the plane: Check point of intersection
					double fraction = (z1[mito_counter2]-z_old)/(xmemory->z - z_old);
					double intersectx = fraction * (xmemory->x - x_old) + x_old;
					double intersecty = fraction * (xmemory->y - y_old) + y_old;
					double intersectz = z1[mito_counter2];
					// If intersection point is on mito
					if (intersectx > x1[mito_counter2] && intersectx < x2[mito_counter2] && intersecty > y1[mito_counter2] && intersecty < y3[mito_counter2])
					{
						foundedge = 5;
						edges_found++;
						saveintersectx = intersectx;
						saveintersecty = intersecty;
						saveintersectz = intersectz;
					}
				}
				// Check Edge 6: Defined by vertices 5, 6, 7 & 8	B
				if (z_old >= z5[mito_counter2] && z5[mito_counter2] >= xmemory->z)
				{
					// Crosses the plane: Check point of intersection
					double fraction = (z5[mito_counter2]-z_old)/(xmemory->z - z_old);
					double intersectx = fraction * (xmemory->x - x_old) + x_old;
					double intersecty = fraction * (xmemory->y - y_old) + y_old;
					double intersectz = z5[mito_counter2];
					// If intersection point is on mito
					if (intersectx > x5[mito_counter2] && intersectx < x6[mito_counter2] && intersecty > y5[mito_counter2] && intersecty < y7[mito_counter2])
					{
						foundedge = 6;
						edges_found++;
						saveintersectx = intersectx;
						saveintersecty = intersecty;
						saveintersectz = intersectz;
					}
				}
				
				// Reflect according
				if (foundedge == 1) //C
				{
					xmemory->y = y1[mito_counter2] - (xmemory->y - y1[mito_counter2]);
				}
				if (foundedge == 2)	//F
				{
					xmemory->x = x2[mito_counter2] + (x2[mito_counter2] - xmemory->x);
				}
				if (foundedge == 3)	//D
				{
					xmemory->y = y3[mito_counter2] + (y3[mito_counter2] - xmemory->y);
				}
				if (foundedge == 4)	//E
				{
					xmemory->x = x1[mito_counter2] - (xmemory->x - x1[mito_counter2]);
				}
				if (foundedge == 5)	//A
				{
					xmemory->z = z1[mito_counter2] - (xmemory->z - z1[mito_counter2]);
				}
				if (foundedge == 6)	//B
				{
					xmemory->z = z5[mito_counter2] + (z5[mito_counter2] - xmemory->z);
				}		
				
			}
			
			if (foundedge == 0)
			{
				mito_counter++;
				mito_counter2++;
			}
			else { mito_counter2 = mito_number; }
		}	
		
		/* periodic boundary conditions */	 
        if (xmemory->x > maxintracellular_xboundary)
        {
            xmemory->x = minintracellular_xboundary + (xmemory->x - maxintracellular_xboundary);
        }
        if (xmemory->x < minintracellular_xboundary)
        {
            xmemory->x = maxintracellular_xboundary - (minintracellular_xboundary - xmemory->x);
        }
        if (xmemory->y > maxintracellular_yboundary)
        {
            xmemory->y = minintracellular_yboundary + (xmemory->y - maxintracellular_yboundary);
        }
        if (xmemory->y < minintracellular_yboundary)
        {
            xmemory->y = maxintracellular_yboundary - (minintracellular_yboundary - xmemory->y);
        }
        if (xmemory->z > maxintracellular_zboundary)
        {
            xmemory->z = minintracellular_zboundary + (xmemory->z - maxintracellular_zboundary);
        }
        if (xmemory->z < minintracellular_zboundary)
        {
            xmemory->z = maxintracellular_zboundary - (minintracellular_zboundary - xmemory->z);
        }
		
		counter++;

		/* Check if Ligand (MCL-1) binds to the Mito or is reflected: check with binding probability */
		if ((foundedge != 0) && (xmemory->state == 2))
		{
			double random = rnd<CONTINUOUS>(rand48);
			// ligand is binding at the intersection point and will be safed as receptor agent
			if (random < probability)				
			{
				xmemory->x = saveintersectx;
				xmemory->y = saveintersecty;
				xmemory->z = saveintersectz;
				xmemory->foundedge = foundedge;
				xmemory->mitoid = mito_counter;
				xmemory->state = 20;
				xmemory->active = 0;
				xmemory->Dc = Dc_state_20;
				add_Receptor_agent(Receptor_agents,  xmemory->id, xmemory->state, xmemory->Pending, xmemory->active,  xmemory->x,  xmemory->y,  xmemory->z, xmemory->Dc, xmemory->collision_counter_ligand, xmemory->check_collision, xmemory->foundedge, xmemory->mitoid, xmemory->mitosize);
				counter = brownianloop;
				return 1; 	//Kill agent
			}
		}
		
		/* Check if Ligand (BAX) binds to the Mito or is reflected: check with binding probability */
		if ((foundedge != 0) && (xmemory->state == 4))
		{
			double random = rnd<CONTINUOUS>(rand48);
			// ligand is binding at the intersection point and will be safed as receptor agent
			if (random < probability)				
			{
				xmemory->x = saveintersectx;
				xmemory->y = saveintersecty;
				xmemory->z = saveintersectz;
				xmemory->foundedge = foundedge;
				xmemory->mitoid = mito_counter;
				xmemory->state = 41;
				xmemory->active = 0;
				xmemory->Dc = Dc_state_41;
				add_Receptor_agent(Receptor_agents,  xmemory->id, xmemory->state, xmemory->Pending, xmemory->active,  xmemory->x,  xmemory->y,  xmemory->z, xmemory->Dc, xmemory->collision_counter_ligand, xmemory->check_collision, xmemory->foundedge, xmemory->mitoid, xmemory->mitosize);
				counter = brownianloop;
				return 1; 	//Kill agent
			}
		}
    }
	
	/* Add Ligand Location Message */
	add_Ligand_Location_message (Ligand_Location_messages, xmemory->id, xmemory->state, xmemory->active, xmemory->x, xmemory->y, xmemory->z, xmemory->Dc, xmemory->collision_counter_ligand, xmemory->check_collision, xmemory->foundedge, xmemory->mitoid, xmemory->mitosize);
return 0;
}


/*****************************************************************************/
/***************************** Receptor Functions ****************************/
/*****************************************************************************/

/***************************/
/***** Receptor Output *****/
/***************************/

__FLAME_GPU_FUNC__ int Receptor_Output(xmachine_memory_Receptor* xmemory, xmachine_memory_Ligand_list* Ligand_agents, xmachine_message_Mito_Location_list* Mito_Location_messages, xmachine_message_Receptor_Location_list* Receptor_Location_messages, RNG_rand48* rand48)
{
	/* the following loop is only to delete the dummy receptors, which was created as a work-around to solve the problem with wrong results when increasing the buffersize. Yes it is nosense that it makes problems,
	but this is necessary to get closer to the correct result */
	if (xmemory->state == 1)	{	return 1;	}
	
	
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
		/* Check if Ligand is going to be a Ligand again or stays on the Mito membrane as receptor agent: check with unbinding probability*/
		double probability;
		probability = 0.00019779;			/* = 1-e^(-0.00199994*0.05) after fitting with FRAP data */
		
		/* Mito energy status is checked -> if Mito is damaged than higher energised and mito_energetics != 0.0; if Mito unchanged/intact than mito_energetics == 0.0; -> Mcl-1 stays on Mito */
		if (mito_energetics[mito_counter] != 0.0)
		{
			probability = 1.0;
		}
		
		/* Mito is normal intact and will check with normal unbinding probability and if random <= unbinding probability -> Mcl-1 becomes a ligand */
		double random = rnd<CONTINUOUS>(rand48);
		if (random <= probability && xmemory->state == 20)
		{
			xmemory->state = 2;
			xmemory->active = 1;
			xmemory->Dc = Dc_state_2;
			add_Ligand_agent(Ligand_agents,  xmemory->id, xmemory->state,  xmemory->Pending, xmemory->active,  xmemory->x,  xmemory->y,  xmemory->z, xmemory->Dc, xmemory->collision_counter_ligand, xmemory->check_collision, xmemory->foundedge, xmemory->mitoid, xmemory->mitosize);
			return 1; 	//Kill agent
		}
		
		/* Mito is normal intact and will check with normal unbinding probability and if random > unbinding probability -> Mcl-1 stays as a receptor and will diffuse in the mito surface */
		/* Diffusion of Mcl-1 on the mitochondria surface if it remains binded */
		/*********************************************************************************************************/
		/*************************************** Diffusion on Mito surface ***************************************/
		/*********************************************************************************************************/
		if ((random > probability && xmemory->state == 20) || xmemory->state != 20)
		{
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
		}
	}
		/* Add Receptor Location Message */
			add_Receptor_Location_message (Receptor_Location_messages, xmemory->id, xmemory->state, xmemory->x,  xmemory->y,  xmemory->z, xmemory->collision_counter_ligand, xmemory->check_collision, xmemory->foundedge, xmemory->mitoid, xmemory->mitosize);
  			
return 0;
}



/***************************/
/**** Ligand Dissociate ****/
/***************************/
/* remaining of FaST code; function makes nothing in this model scenario */
__FLAME_GPU_FUNC__ int Ligand_Dissociate(xmachine_memory_Ligand* xmemory, xmachine_memory_Ligand_list* Ligand_agents, RNG_rand48* rand48)
{
return 0;
}


/***************************/
/*** Ligand Bind Ligand ****/
/***************************/
/* remaining of FaST code; function makes nothing in this model scenario */
__FLAME_GPU_FUNC__ int Ligand_Bind_Ligand(xmachine_memory_Ligand* xmemory, xmachine_memory_Ligand_list* Ligand_agents, xmachine_message_Ligand_Location_list* Ligand_Location_messages, xmachine_message_Ligand_Location_PBM* partition_matrix, xmachine_message_Ligand_byLigand_Bound_list* Ligand_byLigand_Bound_messages, RNG_rand48* rand48)
{
return 0;
}


/*********************************/
/** Ligand Check Bound Receptor **/
/*********************************/

__FLAME_GPU_FUNC__ int Ligand_Check_Bound_Receptor(xmachine_memory_Ligand* xmemory, xmachine_memory_Ligand_list* Ligand_agents, xmachine_message_Ligand_byReceptor_Bound_list* Ligand_byReceptor_Bound_messages, xmachine_message_MReceptor_Confirm_L_list* MReceptor_Confirm_L_messages, RNG_rand48* rand48)
{
    int closest_bound = -1;
    double closest_distance = 99999.0;
    int Reaction = 0;
    xmachine_message_Ligand_byReceptor_Bound* Ligand_byReceptor_Bound_message;
    Ligand_byReceptor_Bound_message = get_first_Ligand_byReceptor_Bound_message(Ligand_byReceptor_Bound_messages);
    while (Ligand_byReceptor_Bound_message)
    {
        if (Ligand_byReceptor_Bound_message->id == xmemory->id)
        {
            if (Ligand_byReceptor_Bound_message->distance < closest_distance)
            {
                 closest_distance = Ligand_byReceptor_Bound_message->distance;
                 closest_bound = Ligand_byReceptor_Bound_message->idfrom;
                 Reaction = Ligand_byReceptor_Bound_message->Reaction;
            }
        }
    Ligand_byReceptor_Bound_message = get_next_Ligand_byReceptor_Bound_message(Ligand_byReceptor_Bound_message, Ligand_byReceptor_Bound_messages);
    }
    
    
    if (closest_bound != -1)
    {
        add_MReceptor_Confirm_L_message(MReceptor_Confirm_L_messages, closest_bound, Reaction);
        return 1;
    }

return 0;
}


/*******************************/
/** Ligand Check Bound Ligand **/
/*******************************/
/* remaining of FaST code; function makes nothing in this model scenario */
__FLAME_GPU_FUNC__ int Ligand_Check_Bound_Ligand(xmachine_memory_Ligand* xmemory, xmachine_memory_Ligand_list* Ligand_agents, xmachine_message_Ligand_byLigand_Bound_list* Ligand_byLigand_Bound_messages, xmachine_message_MLigand_Confirm_list* MLigand_Confirm_messages, RNG_rand48* rand48)
{
return 0;
}


/***************************/
/****** Ligand Confirm *****/
/***************************/
/* remaining of FaST code; function makes nothing in this model scenario */
__FLAME_GPU_FUNC__ int Ligand_Confirm(xmachine_memory_Ligand* xmemory, xmachine_memory_Ligand_list* Ligand_agents, xmachine_message_MLigand_Confirm_list* MLigand_Confirm_messages, RNG_rand48* rand48)
{
return 0;
}


/**********************************/
/*** Receptor Dissociate Ligand ***/
/**********************************/

__FLAME_GPU_FUNC__ int Receptor_Dissociate_Ligand(xmachine_memory_Receptor* xmemory, xmachine_memory_Ligand_list* Ligand_agents, RNG_rand48* rand48)
{
    /* Dissociation of aBax2 into Bax and aBax */
    if (xmemory->state == 42 && xmemory->active == 1)
    {
        double random = rnd<CONTINUOUS>(rand48);
        if (random < kr_reaction_1)
        {
            xmemory->state = 41;
            double x = xmemory->x;
            double y = xmemory->y;
            double z = xmemory->z;
            double theta = 2 * pi * rnd<CONTINUOUS>(rand48);
            double dx = ubr_reaction_1 * sin(theta);
            double dy = 0;
            double dz = ubr_reaction_1 * cos(theta);
            int active = 1;
            xmemory->active = 0;
            xmemory->Dc = Dc_state_41;
            add_Ligand_agent (Ligand_agents, xmemory->id + 25000, 4, -1, active, x+dx, y+dy, z+dz, Dc_state_4, xmemory->collision_counter_ligand, xmemory->check_collision, xmemory->foundedge, xmemory->mitoid, xmemory->mitosize);
        }
    }


    /* Dissociation of Mcl1aBax into Mcl1 and aBax */
     if (xmemory->state == 21 && xmemory->active == 1)
    {
        double random = rnd<CONTINUOUS>(rand48);
        if (random < kr_reaction_10)
        {
            xmemory->state = 41;
            double x = xmemory->x;
            double y = xmemory->y;
            double z = xmemory->z;
            double theta = 2 * pi * rnd<CONTINUOUS>(rand48);
            double dx = ubr_reaction_10 * sin(theta);
            double dy = 0;
            double dz = ubr_reaction_10 * cos(theta);
            int active = 1;
            xmemory->active = 0;
            xmemory->Dc = Dc_state_41;
            add_Ligand_agent (Ligand_agents, xmemory->id + 25000, 2, -1, active, x+dx, y+dy, z+dz, Dc_state_2, xmemory->collision_counter_ligand, xmemory->check_collision, xmemory->foundedge, xmemory->mitoid, xmemory->mitosize);
        }
    } 

return 0;
}


/**********************************/
/*** Receptor Dissociate Receptor ***/
/**********************************/
/* remaining of FaST code; function makes nothing in this model scenario */
__FLAME_GPU_FUNC__ int Receptor_Dissociate_Receptor(xmachine_memory_Receptor* xmemory, xmachine_memory_Receptor_list* Receptor_agents, RNG_rand48* rand48)
{
return 0;
}


/***************************/
/*** Receptor Bind Ligand **/
/***************************/

__FLAME_GPU_FUNC__ int Receptor_Bind_Ligand(xmachine_memory_Receptor* xmemory, xmachine_message_Ligand_Location_list* Ligand_Location_messages, xmachine_message_Ligand_Location_PBM* partition_matrix, xmachine_message_Ligand_byReceptor_Bound_list* Ligand_byReceptor_Bound_messages, RNG_rand48* rand48)
{
    double closestdistance = 99999.0;
    int Bound = -1;
    int Boundstate = 0;
    int Reaction = 0;

    double distance = 0.0;

        xmachine_message_Ligand_Location* Ligand_Location_message;
        Ligand_Location_message = get_first_Ligand_Location_message(Ligand_Location_messages, partition_matrix, xmemory->x,  xmemory->y,  xmemory->z);
        while (Ligand_Location_message)
        {
    if (xmemory->state == 41)
    {
            if (Ligand_Location_message->state == 2)
            {
                distance = sqrt((xmemory->x - Ligand_Location_message->x)*(xmemory->x - Ligand_Location_message->x)+(xmemory->y - Ligand_Location_message->y)*(xmemory->y - Ligand_Location_message->y)+(xmemory->z - Ligand_Location_message->z)*(xmemory->z - Ligand_Location_message->z));
                if (distance < br_reaction_10 && distance < closestdistance)
                {
                    closestdistance = distance;
                    Bound = Ligand_Location_message->id;
                    Boundstate = 2;
                    Reaction = 10;
                }
            }
    }
            Ligand_Location_message = get_next_Ligand_Location_message(Ligand_Location_message, Ligand_Location_messages, partition_matrix);
        }
    if (Bound != -1)
    {
        xmemory->Pending = Bound;
        add_Ligand_byReceptor_Bound_message (Ligand_byReceptor_Bound_messages, Bound, xmemory->id, closestdistance, Reaction);
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
    if (xmemory->state == 20)
    {
            if (Receptor_Location_message->state == 5 && Receptor_Location_message->id != xmemory->id)
            {
                distance = sqrt((xmemory->x - Receptor_Location_message->x)*(xmemory->x - Receptor_Location_message->x)+(xmemory->y - Receptor_Location_message->y)*(xmemory->y - Receptor_Location_message->y)+(xmemory->z - Receptor_Location_message->z)*(xmemory->z - Receptor_Location_message->z));
                if (distance < br_reaction_8 && distance < closestdistance)
                {
                    closestdistance = distance;
                    Bound = Receptor_Location_message->id;
                    Boundstate = 5;
                    Reaction = 8;
                }
            }
            if (Receptor_Location_message->state == 41 && Receptor_Location_message->id != xmemory->id)
            {
                distance = sqrt((xmemory->x - Receptor_Location_message->x)*(xmemory->x - Receptor_Location_message->x)+(xmemory->y - Receptor_Location_message->y)*(xmemory->y - Receptor_Location_message->y)+(xmemory->z - Receptor_Location_message->z)*(xmemory->z - Receptor_Location_message->z));
                if (distance < br_reaction_11 && distance < closestdistance)
                {
                    closestdistance = distance;
                    Bound = Receptor_Location_message->id;
                    Boundstate = 41;
                    Reaction = 11;
                }
            }
    }

    if (xmemory->state == 41)
    {
            if (Receptor_Location_message->state == 41 && Receptor_Location_message->id != xmemory->id)
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
            if (Receptor_Location_message->state == 20 && Receptor_Location_message->id != xmemory->id)
            {
                distance = sqrt((xmemory->x - Receptor_Location_message->x)*(xmemory->x - Receptor_Location_message->x)+(xmemory->y - Receptor_Location_message->y)*(xmemory->y - Receptor_Location_message->y)+(xmemory->z - Receptor_Location_message->z)*(xmemory->z - Receptor_Location_message->z));
                if (distance < br_reaction_11 && distance < closestdistance)
                {
                    closestdistance = distance;
                    Bound = Receptor_Location_message->id;
                    Boundstate = 20;
                    Reaction = 11;
                }
            }
    }

    if (xmemory->state == 42)
    {
            if (Receptor_Location_message->state == 42 && Receptor_Location_message->id != xmemory->id)
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
            if (Receptor_Location_message->state == 44 && Receptor_Location_message->id != xmemory->id)
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
            if (Receptor_Location_message->state == 42 && Receptor_Location_message->id != xmemory->id)
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
                 if (Receptor_Bound_message->Reaction == 11)
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
                if (Reaction == 2)
                {
                    if (xmemory->state == 41 && closest_bound < xmemory->id)
                    {
                        xmemory->state = 42;
                        xmemory->Pending = -1;
                    }
                    if (xmemory->state == 41 && closest_bound > xmemory->id)
                    {
                        return 1;
                    }
                }
                if (Reaction == 3)
                {
                    if (xmemory->state == 42 && closest_bound < xmemory->id)
                    {
                        xmemory->state = 44;
                        xmemory->Pending = -1;
                    }
                    if (xmemory->state == 42 && closest_bound > xmemory->id)
                    {
                        return 1;
                    }
                }
                if (Reaction == 4)
                {
                    if (xmemory->state == 42)
                    {
                        xmemory->state = 46;
                        xmemory->Pending = -1;
                    }
                    if (xmemory->state == 44)
                    {
                        return 1;
                    }
                }
                if (Reaction == 11)
                {
                    if (xmemory->state == 20)
                    {
                        xmemory->state = 21;
                        xmemory->Pending = -1;
                    }
                    if (xmemory->state == 41)
                    {
                        return 1;
                    }
                }
        }

return 0;
}


/**********************************/
/***** Receptor Confirm Ligand ****/
/**********************************/

__FLAME_GPU_FUNC__ int Receptor_Confirm_Ligand(xmachine_memory_Receptor* xmemory, xmachine_memory_Receptor_list* Receptor_agents, xmachine_message_MReceptor_Confirm_L_list* MReceptor_Confirm_L_messages, RNG_rand48* rand48)
{
    xmachine_message_MReceptor_Confirm_L* MReceptor_Confirm_L_message;
    MReceptor_Confirm_L_message = get_first_MReceptor_Confirm_L_message(MReceptor_Confirm_L_messages);
    while (MReceptor_Confirm_L_message)
    {
        if (xmemory->Pending != -1)
        {
            if (MReceptor_Confirm_L_message->id == xmemory->id)
            {
                if (MReceptor_Confirm_L_message->Reaction == 10)
                {
                    if (xmemory->state == 41)
                    {
                        xmemory->state = 21;
                        xmemory->Pending = -1;
                    }
                }
            }
        }
        MReceptor_Confirm_L_message = get_next_MReceptor_Confirm_L_message(MReceptor_Confirm_L_message, MReceptor_Confirm_L_messages);
    }

return 0;
}
/************************************/
/***** Receptor Confirm Receptor ****/
/************************************/

__FLAME_GPU_FUNC__ int Receptor_Confirm_Receptor(xmachine_memory_Receptor* xmemory, xmachine_memory_Receptor_list* Receptor_agents, xmachine_message_MReceptor_Confirm_list* MReceptor_Confirm_messages, RNG_rand48* rand48)
{
    int die = 0;
	xmachine_message_MReceptor_Confirm* MReceptor_Confirm_message;
    MReceptor_Confirm_message = get_first_MReceptor_Confirm_message(MReceptor_Confirm_messages);
    while (MReceptor_Confirm_message)
    {
        if (xmemory->Pending != -1)
        {
            if (MReceptor_Confirm_message->id == xmemory->id)
            {
                if (MReceptor_Confirm_message->Reaction == 2)
                {
                    if (xmemory->state == 41)
                    {
                        xmemory->state = 42;
                        xmemory->Pending = -1;
                    }
                    if (xmemory->state == 41)
                    {
                        die = 1;
                    }
                }
                if (MReceptor_Confirm_message->Reaction == 3)
                {
                    if (xmemory->state == 42)
                    {
                        xmemory->state = 44;
                        xmemory->Pending = -1;
                    }
                    if (xmemory->state == 42)
                    {
                        die = 1;
                    }
                }
                if (MReceptor_Confirm_message->Reaction == 4)
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
                if (MReceptor_Confirm_message->Reaction == 11)
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
        }
        MReceptor_Confirm_message = get_next_MReceptor_Confirm_message(MReceptor_Confirm_message, MReceptor_Confirm_messages);
    }
	if (die == 1)
	{
		return 1;
	}

xmemory->Pending = -1;
xmemory->active = 1;
return 0;
}
#endif //_FLAMEGPU_FUNCTIONS
