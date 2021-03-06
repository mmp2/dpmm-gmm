/* Defines */

#ifdef __GNUC__
#define inline __inline__
#else
#define inline __inline
#endif

/* Number of items between c resample print notifications */
#define C_NOTIFICATION_PERIOD 1000000

/* Number of items between sigma resample print notifications */
#define SIGMA_NOTIFICATION_PERIOD 10000

/* Maximum number of times to run slice sampler before aborting */
#define SLICE_SAMPLE_RETRY_LIMIT 3

/* Macro for converting positive double to int */
#define saferound(x) (unsigned int)((x) + 0.5)

/* Model types */
#define BETA_GIBBS   0
#define SLICE_GIBBS  1

/* Locations of parameters */

/* inv_pi:     top-t orderings, uint32(max_t, N) matrix */
#define R_INV_PI          0

/* t:          Length of each ordering, uint32(N) vector */
#define R_T               1

/* n:          Maximum observed item being ranked, uint32 scalar */
#define R_N               2

/* c:          Cluster assignments, uint32(N) vector */
#define R_C               3

/* sigma:      Cluster centroids in rank (not ordering) format, uint32(n, N) matrix */
/* inv_sigma holds the sigma's in ordering format */
#define R_SIGMA           4

/* rho:        Cluster dispersion parameters, double(max_t, N) matrix */
#define R_RHO             5

/* nc:         Cache of number of occurrences of each cluster at each rank, uint32(max_t, N) matrix */
/*             Pass in empty matrix to construct initially */
#define R_NC              6

/* S:          Sum of lower triangles of sufficient stats for each cluster, double(max_t, N) matrix */
/*             Pass in empty matrix to construct initially */
#define R_S               7

/* r_0:        Prior hyperparameter r, double(max_t) vector */
#define R_R_0             8

/* nu_0:       Prior hyperparameter nu (strength), double scalar */
#define R_NU_0            9

/* alpha_0:    Prior hyperparameter alpha (concentration parameter), double scalar */
#define R_ALPHA_0         10

/* R:          Sufficient statistics of each permutation; cell(max_t, N) of uint32(x) vectors where x can be any value */
/*             Vector entries are linear indexes in the actual R matrix that should be 1 */
/*             Must be constructed externally in MATLAB */
#define R_PI_R            11

/* beta_table: Precomputed table of probabilities for sampling sigma from one permutation; double(n, max_t) matrix */
/*             Must be constructed externally in MATLAB */
#define R_BETA_TABLE      12

/* Temperature for sampling discrete variables */
#define R_TEMPERATURE     13

/* Number of iterations to perform for sampling a new sigma for clusters larger than one */
#define R_SIGMA_GIBBS     14

/* Number of iterations to perform for slice sampling a new rho */
#define R_RHO_SLICES      15

/* Total number of iterations to perform */
#define R_ITERATIONS      16

/* Iteration number to start on (used only for display) */
#define R_STARTITERATIONS 17

/* Model specification parameter, see above */
#define R_MODEL           18

/* Cluster assignments for each permutation, see R_C */
#define L_C               0

/* Cluster centroid permutations, see R_SIGMA */
#define L_SIGMA           1

/* Cluster dispersions, see R_RHO */
#define L_RHO             2

/* Cluster counts, see R_NC */
#define L_NC              3

/* Sum of sufficient stats, see R_S */
#define L_S               4

/* Log probability associated with the last iteration */
#define L_LOGPROB         5
