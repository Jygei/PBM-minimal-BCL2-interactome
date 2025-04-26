#include <stdio.h>

/**************************/
/*******Ligand States******/
/**************************/

// State 1 MCL-1

/**************************/
/******Receptor States*****/
/**************************/

// State 10 MCL-1

/**************************/
/*********Reactions********/
/**************************/


/**************************/
/****Global Definitions****/
/**************************/

#define pi 3.141592
#define brownianloop 500				// has to be changed if the reaction rate changes

#define fragmentation 43.84132496		// fragmentation degree: how many frag. mitochondria fit into on non-frag. mito.?
#define fragmentation_percent 0.1		// fragmentation amount: how many mitochondria of total mitochondria number are fragmented? 0.1 = 10% 
#define mito_mean_number 45				// total mitochondria number before fragmentation; based on in-house measurements of NCI-H460 
#define mito_number ((unsigned long) (mito_mean_number - fragmentation + mito_mean_number*fragmentation_percent*fragmentation))		// mitochondrial number after fragmentation


/* G1 */	// model cell dimensions based on the in-house measurements of the mean NCI-H460 cell in G1 phase

#define minintracellular_xboundary 0.0
#define minintracellular_yboundary 0.0
#define minintracellular_zboundary 0.0
#define maxintracellular_xboundary 10000
#define maxintracellular_yboundary 4500
#define maxintracellular_zboundary 87053.77778

#define minextracellular_xboundary 0.0
#define minextracellular_yboundary 9000
#define minextracellular_zboundary 0.0
#define maxextracellular_xboundary 10000
#define maxextracellular_yboundary 9000
#define maxextracellular_zboundary 87053.77778

/*	nucleus/blackbox in the middle of the model cell where mitos can't be created -> necessary for simulations of different extreme cases of mitochondrial arrangement	*/
/*
#define minnucleus_xboundary 0		//((unsigned long) (maxintracellular_xboundary / 4))
#define minnucleus_yboundary 0		//((unsigned long) (maxintracellular_yboundary / 4))
#define minnucleus_zboundary 10000		//((unsigned long) (maxintracellular_zboundary / 4))
#define maxnucleus_xboundary 10000		//((unsigned long) (maxintracellular_xboundary / 4 * 3))
#define maxnucleus_yboundary 4500		//((unsigned long) (maxintracellular_yboundary / 4) * 3)
#define maxnucleus_zboundary 77053		//((unsigned long) (maxintracellular_zboundary / 4 * 3))
*/
