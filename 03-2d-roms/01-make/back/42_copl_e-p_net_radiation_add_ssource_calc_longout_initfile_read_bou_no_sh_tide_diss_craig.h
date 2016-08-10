/*******************************************************************************/
/*                                                                             */
/* What's new:                                                                 */
/*
/* Last Revision : 8-November-2014                                             */
/*******************************************************************************/

#define ROMS_MODEL
#define SWAN_MODEL
#define MCT_LIB

#undef  WEC_MELLOR        /* activate radiation stress terms from Mellor 08    */
#define WEC_VF            /* activ wave-current stresses from Uchiyama et al.  */
#define UV_KIRBY          /* compute depth-avg current based on Hwave that will 
                             sent from ocn to the wave model for coupling      */

#define CRAIG_BANNER      /*  use if Craig and Banner wave breaking surface flux    **/
#define ZOS_HSIG          /*  use if surface roughness from wave amplitude    **/
#define KANTHA_CLAYSON    /*  use if Kantha and Clayson stability function    **/

/** CANUTO_A            use if Canuto A-stability function formulation        **
 ** CANUTO_B            use if Canuto B-stability function formulation        **
 ** CHARNOK             use if Charnok surface roughness from wind stress     **
**#define WDISS_WAVEMOD     /*  activate wave dissipation from a wave model      *
 ** KANTHA_CLAYSON      use if Kantha and Clayson stability function          **
 ** K_C2ADVECTION       use if 2nd-order centered advection                   **
 ** K_C4ADVECTION       use if 4th-order centered advection                   **
 ** N2S2_HORAVG         use if horizontal smoothing of buoyancy/shear         **
 ** #define SURFACE_STREAMING /*  activate wave enhanced surface streaming        **
 ** TKE_WAVEDISS        use if wave breaking surface flux from wave amplitude **/


/* define only one of the following 5 */
#undef UV_LOGDRAG
#define UV_QDRAG
#undef MB_BBL
#undef SG_BBL
#undef SSW_BBL
#ifdef SSW_BBL
# define SSW_CALC_ZNOT
#endif


/* --------------- Options associated with momentum equations ---------------- */
#define UV_ADV            /* turn ON or OFF advection terms                    */ 
#define UV_VIS2           /* turn ON or OFF harmonic horizontal mixing         */
/*#define UV_QDRAG          /* turn ON or OFF quadratic bottom friction          */
#define UV_COR            /* use to activate Coriolis term                     */
#undef  UV_PSOURCE        /* use to activate horizontal point sources/sinks   <<<< is this for river */

/* ---------------- Options associated with tracer equations ----------------- */
#define TS_DIF2           /* use to turn ON and OFF harmonic horizontal mixing */
#define TS_U3HADVECTION   /* use if 3rd order upstream horizontal advection    */
#define TS_C4VADVECTION   /* use if 4th order centered central advection       */
#define SALINITY          /* use if having salinity                            */
#undef TS_PSOURCE        /* use to turn ON and OFF point sources/sinks        */

/* --------------------- Options for Lagrangian drifters --------------------- */
#undef FLOATS            /* use to activate simulated Lagrangian drifters     */

/*------------------ Options for pressure gradient algorithm ------------------*/
#define DJ_GRADPS         /* splines density  Jacobian (Shchepetkin, 2000)     */

/*------------------ Options for atmospheric boundary layer -------------------*/
#define BULK_FLUXES       /* use if bulk fluxes computation                    */

#undef USE_WIND_STRESS    /* if activated, use sustr, svstr from input file instead of Uwind, Vwind */

/* --------------------- Options for model configuration --------------------- */
#define SOLVE3D
#define MASKING           /* use if land/sea masking                           */
#define NONLIN_EOS
/*#undef  STATIONS          /* use if writing out station data                   */

/* --------------- Options for analytical field configuration ---------------- */
#define ANA_BTFLUX        /* use if analytical bottom temperature flux         */
#define ANA_BSFLUX        /* use if analytical bottom salinity flux            */

#define  BULK_FLUXES
# ifdef   BULK_FLUXES
#  undef  LONGWAVE     /* undef to read net longwave, define to use Berliand */
#  define LONGWAVE_OUT /* define to read downward longwave, compute outgoing */
#  undef  DIURNAL_SRFLUX  /* impose shortwave radiation local diurnal cycle */
#  undef  ANA_SRFLUX   /* analytical surface shortwave radiation flux */
#  undef  ALBEDO       /* use albedo equation for shortwave radiation */
#  undef  ANA_CLOUD    /* analytical cloud fraction => zero cloud     */
#  define  SOLAR_SOURCE /* solar shortwave distributed over water column */
#  define EMINUSP      /* turn ON internal calculation of E-P */
#  undef  ANA_RAIN     /* zero rain, with eminusp option can use ncep rain data  */
# else
#  define ANA_SMFLUX /* analytical surface momentum stress */
#  define ANA_STFLUX /* analytical surface temperature flux */
#  define ANA_SSFLUX /*  analytical surface salinity flux */ /* JO: I 
don't think this means anything */
# endif   /* BULK_FLUXES if */

/* ---------------- Options for horizontal mixing of momentum ----------------- */
#define MIX_S_UV          /* use if mixing along constant S-surfaces            */

/* ----------------- Options for horizontal mixing of tracers ----------------- */
#undef MIX_S_TS           /* use if mixing along constant S-surfaces            */
#define MIX_ISO_TS        /* use if mixing along constant epineutral surfaces   */

/* ---------- Options for vertical mixing of momentum and tracers ------------- */
/* Activate only one closure!                                                   */
#define GLS_MIXING 
#ifdef GLS_MIXING
#define KANTHA_CLAYSON    /* Kantha-Clayson stability function                   */
#endif

#undef MY25_MIXING 
#ifdef MY25_MIXING
#undef KANTHA_CLAYSON    /* Kantha-Clayson stability function                   */
#endif

/* --------------- Options for tidal forcing at open boundaries --------------- */
/*#define RAMP_TIDES        /* use if rampin (over one day) tidal forcing         */
#undef SSH_TIDES         /* SE: was "define": read tidal SSH from file         */
#undef UV_TIDES          /* use if imposing tidal currents                     */

/* ----- Options for nearshore stresses and shallow water configurations ------ */
#define WET_DRY


/*you need to have either this one or ana OBC ones    */

/*#undef ANA_M2OBC
#ifdef ANA_M2OBC
# define ANA_FSOBC
#else

# define ADD_FSOBC         /*  use to add tidal elevation to processed OBC data*/      
/*# define ADD_M2OBC 

# define WEST_FSOBC
# define EAST_FSOBC
# define SOUTH_FSOBC

# define WEST_M2OBC
# define EAST_M2OBC
# define SOUTH_M2OBC
#endif 
*/

#undef ANA_INITIAL        /* use if analytical initial conditions       */
#undef NRI_29MAY          /* use if analytical initial conditions       */
/* Lateral boundaries ----------------------------------------------------------*/
/* Begin open boundary condition settings --------------------------------------*/

#define  NORTH_M2FLATHER   /* Northern edge, 2D momentum,  Flather condition  */
#define  SOUTH_M2FLATHER   /* Southern edge, 2D momentum,  Flather condition  */
#define  WEST_M2FLATHER    /* Western edge,  2D momentum,  Flather condition  */
#define  EAST_M2FLATHER    /* EASTern edge,  2D momentum,  Flather condition  */

#define  NORTH_FSCHAPMAN   /* Northern edge, free-surface, Chapman condition  */
#define  SOUTH_FSCHAPMAN   /* Southern edge, free-surface, Chapman condition  */
#define  WEST_FSCHAPMAN    /* Western edge,  free-surface, Chapman condition  */ 
#define  EAST_FSCHAPMAN    /* EASTern edge,  free-surface, Chapman condition  */ 


#undef  NORTH_M2RADIATION 
#undef  SOUTH_M2RADIATION
#undef  WEST_M2RADIATION
#undef  EAST_M2RADIATION


#define  RADIATION_2D       /* Tangential phase speed in radiation conditions  */
#define  NORTH_M3RADIATION  /* Northern edge, 3D momentum, radiation condition */
#define  SOUTH_M3RADIATION  /* Southern edge, 3D momentum, radiation condition */
#define  WEST_M3RADIATION   /* Western edge,  3D momentum, radiation condition */
#define  EAST_M3RADIATION   /* EASTern edge,  3D momentum, radiation condition */


#define  NORTH_TRADIATION   /* Northern edge, tracers, radiation condition     */
#define  SOUTH_TRADIATION   /* Southern edge, tracers, radiation condition     */
#define  WEST_TRADIATION    /* Western edge,  tracers, radiation condition     */
#define  EAST_TRADIATION    /* EASTern edge,  tracers, radiation condition     */

#define  NORTH_TNUDGING     /* Northern edge, tracers, passive/active term */
#define  EAST_TNUDGING      /* Western edge, tracers, passive/active term */
#define  SOUTH_TNUDGING     /* Western edge, tracers, passive/active term */
#define  WEST_TNUDGING      /* Western edge, tracers, passive/active term */


/* Other options -------------------------------------------------------------*/
#define DOUBLE_PRECISION   /* double precision arithmetic                     */
#define UV_U3HADVECTION    /* 3rd-order upstream horizontal advection of 3D mom */
#define VAR_RHO_2D         /* variable density barotropic mode                */
#define K_GSCHEME          /* 3rd-order upstream advection of TKE fields      */


/* DIAG  -------------------------------------------------------------*/
#define  DIAGNOSTICS_UV     /*    use if writing out momentum diagnostics               ***/
#define  DIAGNOSTICS_TS     /*   use if writing out tracer diagnostics */




