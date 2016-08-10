
#define ROMS_MODEL
#define CURVGRID

#define UV_VIS2
#define MIX_S_UV
#define MASKING
#undef  WET_DRY
#define UV_ADV

#define ANA_INITIAL
#define ANA_SMFLUX

#undef  SOLVE3D

#define RAMP_TIDES
#define SSH_TIDES
#define UV_TIDES



/* #define ANA_M2OBC         /*  use if analytical 2D momentum boundary conditions     */
/* #define ANA_SSH           /*  use if analytical sea surface height                  */



/* define open boundary   */
/*#define NORTHERN_WALL  */
/*you need to have either this one or ana OBC ones    */
/*#define ADD_FSOBC     */      /*use to add tidal elevation to processed OBC data      */
/*#define ADD_M2OBC

#define WEST_FSOBC
#define EAST_FSOBC
#define SOUTH_FSOBC

#define WEST_M2OBC
#define EAST_M2OBC
#define SOUTH_M2OBC */

/*#define ANA_M2OBC    /*Comment if you like to force boundary condition  */
/*#define ANA_FSOBC  */

/*
#define EAST_FSCHAPMAN
#define EAST_M2FLATHER
#define WEST_FSCHAPMAN
#define WEST_M2FLATHER
#define NORTH_FSCHAPMAN
#define NORTH_M2FLATHER
#define SOUTH_FSCHAPMAN
#define SOUTH_M2FLATHER
*/

#define EAST_FSCLAMPED     /* use if free-surface clamped condition     */
#define EAST_M2CLAMPED     /* use if 2D momentum clamped condition      */*
#define WEST_FSCLAMPED     /* use if free-surface clamped condition     */
#define WEST_M2CLAMPED     /* use if 2D momentum clamped condition      */*
#define NORTH_FSCLAMPED     /* use if free-surface clamped condition     */
#define NORTH_M2CLAMPED     /* use if 2D momentum clamped condition      */*
#define SOUTH_FSCLAMPED     /* use if free-surface clamped condition     */
#define SOUTH_M2CLAMPED     /* use if 2D momentum clamped condition      */*




/*  End of open boundary */

#define  UV_QDRAG

#undef DIAGNOSTICS_UV
#undef AVERAGES

