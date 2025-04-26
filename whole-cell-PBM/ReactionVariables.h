//remains of FaST code; not used in this simulation setup, but were kept in the code for the sake of completeness

/* Reaction 1: Bax+aBax<-->aBax2 */
/* Reaction 1: Reversible Receptor-Ligand Reaction */
#define br_reaction_1 0.0
#define ubr_reaction_1 0.0
#define kr_reaction_1 0.0

/* Reaction 2: aBax+aBax<-->aBax2 */
/* Reaction 2: Reversible Receptor-Receptor Reaction */
/* #define br_reaction_2 0.000003
#define ubr_reaction_2 0.000010
#define kr_reaction_2 0.000000 */
#define br_reaction_2 0.0
#define ubr_reaction_2 0.0
#define kr_reaction_2 0.000000

/* Reaction 3: aBax2+aBax2<-->aBax4 */
/* Reaction 3: Reversible Receptor-Receptor Reaction */
/* #define br_reaction_3 0.002570
#define ubr_reaction_3 0.010281
#define kr_reaction_3 0.012071 */
#define br_reaction_3 0.0
#define ubr_reaction_3 0.0
#define kr_reaction_3 0.0

/* Reaction 4: aBax2+aBax4<-->aBax6 */
/* Reaction 4: Reversible Receptor-Receptor Reaction */
/* #define br_reaction_4 0.000150
#define ubr_reaction_4 0.000599
#define kr_reaction_4 0.000066 */
#define br_reaction_4 0.0
#define ubr_reaction_4 0
#define kr_reaction_4 0.0

/* Reaction 5: tBid+Bax<-->tBidBax */
/* Reaction 5: Reversible Receptor-Ligand Reaction */
#define br_reaction_5 0.340957
#define ubr_reaction_5 1.363828
#define kr_reaction_5 0.237302

/* Reaction 6: tBid+Bclxl<-->tBidBclxl */
/* Reaction 6: Reversible Receptor-Ligand Reaction */
#define br_reaction_6 0.340957
#define ubr_reaction_6 1.363828
#define kr_reaction_6 0.020127

/* Reaction 7: tBid+Mcl1<-->tBidMcl1 */
/* Reaction 7: Reversible Receptor-Ligand Reaction */
#define br_reaction_7 0.340957
#define ubr_reaction_7 1.363828
#define kr_reaction_7 0.020127

/* Reaction 8: tBid+Mcl1Mito<-->tBidMcl1 */
/* Reaction 8: Reversible Receptor-Receptor Reaction */
#define br_reaction_8 0.002570
#define ubr_reaction_8 0.010281
#define kr_reaction_8 0.020127

/* Reaction 9: Bclxl+aBax<-->BclxlaBax */
/* Reaction 9: Reversible Receptor-Ligand Reaction */
#define br_reaction_9 0.340957
#define ubr_reaction_9 1.363828
#define kr_reaction_9 0.004988

/* Reaction 10: Mcl1+aBax<-->Mcl1aBax */
/* Reaction 10: Reversible Receptor-Ligand Reaction */
#define br_reaction_10 1.5303
#define ubr_reaction_10 0.0
#define kr_reaction_10 0.0

/* Reaction 11: Mcl1Mito+aBax<-->Mcl1aBax */
/* Reaction 11: Reversible Receptor-Receptor Reaction */
#define br_reaction_11 1.5303
#define ubr_reaction_11 0.0
#define kr_reaction_11 0.0

/* Reaction 12: BclxlaBax-->Bclxl+Bax */
/* Reaction 12: Catalysis */
#define kf_reaction_12 0.000288

/* Reaction 13: Mcl1aBax-->Mcl1+Bax */
/* Reaction 13: Catalysis */
#define kf_reaction_13 0.0

/* Reaction 14: Mcl1aBax-->Mcl1Mito+Bax */
/* Reaction 14: Catalysis */
#define kf_reaction_14 0.0

/* Reaction 15: tBidBax-->tBid+aBax */
/* Reaction 15: Catalysis */
#define kf_reaction_15 0.000010