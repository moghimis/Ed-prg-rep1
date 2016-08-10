/* 6 Nov 2012 */ 
/* for local model  including BOU from inner model */
#define ROMS_MODEL
#undef CURVGRID

#define UV_VIS2
#define MIX_S_UV
#define MASKING
#define WET_DRY
#define UV_ADV

#define ANA_INITIAL
#define ANA_SMFLUX

#undef SOLVE3D


/* no analitical tide */
#undef RAMP_TIDES
#undef SSH_TIDES
#undef UV_TIDES

/* if you have only analytical tide */
/*#define ANA_M2OBC    Comment if you like to force boundary condition  */  
/*#define ANA_FSOBC   */ 

/* define open boundary   */
/*you need to have either this one or ana OBC ones    */
#define ADD_FSOBC           /*use to add tidal elevation to processed OBC data      */
#define ADD_M2OBC 

#define WEST_FSOBC
#define EAST_FSOBC
#define SOUTH_FSOBC
#define NORTH_FSOBC

#define WEST_M2OBC
#define EAST_M2OBC
#define SOUTH_M2OBC 
#define NORTH_M2OBC


/*   1st try  */
/*
#define EAST_FSCHAPMAN
#define WEST_FSCHAPMAN
#define NORTH_FSCHAPMAN
#define SOUTH_FSCHAPMAN

#define EAST_M2FLATHER
#define SOUTH_M2FLATHER
#define WEST_M2FLATHER
#define NORTH_M2FLATHER
*/


/*   2ndEAST_M2NUDGING try  */
#define EAST_FSCHAPMAN
#define WEST_FSCHAPMAN
#define NORTH_FSCHAPMAN
#define SOUTH_FSCHAPMAN

#define EAST_M2FLATHER
#define SOUTH_M2FLATHER
#define WEST_M2FLATHER
#define NORTH_M2FLATHER

#define EAST_M2NUDGING
#define SOUTH_M2NUDGING
#define WEST_M2NUDGING
#define NORTH_M2NUDGING

/*  End of open boundary */




#define  UV_QDRAG
#define  DIAGNOSTICS_UV
#undef   AVERAGES   /*  To write out average quantities */

