
#define  ROMS_MODEL
#define  CURVGRID           /*  Orthogonal curvilinear grid. */
#define  MASKING            /*  Land/Sea masking. */
#define  NONLINEAR          /*  Nonlinear Model. */

#define  UV_ADV             /*  Advection of momentum. */
#define  UV_QDRAG           /*  Quadratic bottom stress. */
#define  UV_VIS2            /*  Harmonic mixing of momentum. */


#define  SSH_TIDES          /*  Add tidal elevation to SSH climatology. */
#define  UV_TIDES           /*  Add tidal currents to 2D momentum climatologies. */
#define  RAMP_TIDES         /*  Ramping tidal forcing for one day. */


#define  ANA_FSOBC          /*  Analytical free-surface boundary conditions. */
#define  ANA_M2OBC          /*  Analytical 2D momentum boundary conditions. */

#define  ANA_INITIAL        /*  Analytical initial conditions. */
#define  ANA_SMFLUX         /*  Analytical kinematic surface momentum flux. */

#define  WEST_FSCLAMPED     /*  Western edge, free-surface, Clamped condition. */
#define  WEST_M2CLAMPED    /*   Western edge, 2D momentum, Clamped condition. */
#define  SOUTH_FSCLAMPED    /*  Southern edge, free-surface, Clamped condition. */
#define  SOUTH_M2CLAMPED    /*  Southern edge, 2D momentum, Clamped condition. */
#define  NORTH_FSCLAMPED    /*  Northern edge, free-surface, Clamped condition. */
#define  NORTH_M2CLAMPED    /*  Northern edge, 2D momentum, Clamped condition. */
#define  EAST_FSCLAMPED     /*  Eastern edge, free-surface, Clamped condition. */
#define  EAST_M2CLAMPED     /*  Eastern edge, 2D momentum, Clamped condition. */
