#ifndef CPNX_H
#define CPNX_H

#ifdef _WIN32
    #ifdef BUILDING_CPNX
        #define CPNX_API __declspec(dllexport)
    #else
        #define CPNX_API __declspec(dllimport)
    #endif
#else
    #define CPNX_API
#endif

#include <math.h>
#include <stddef.h>

#define NO_ERROR          0
#define VBAP_INIT_ERROR   1
#define VBAP2D_INIT_ERROR 2
#define VBAP3D_INIT_ERROR 3
#define DBAP_INIT_ERROR   4

#define PI (double)acos(-1.0)
#define TO_RAD(x) ((x) * PI / 180.0)
#define TO_DEG(x) ((x) * 180.0 / PI)

typedef enum {
    PLANAR = 0,
    INTERLEAVED = 1,
} CHANNEL_MODE;

typedef enum {
    DEGREE = 0,
    RADIANS = 1
} ANGLE_KIND;

typedef enum {
    PVBAP2D,
    PVBAP3D,
    PDBAP,
} PANNIX_DIMENSIONS;

typedef struct {
    double x;
    double y;
    double z;
} CartesianPoint;

typedef struct {
    double rho;
    double phi;
    double theta;
} PolarPoint;

typedef struct {
    CartesianPoint p;
    int index;
} HullPoint;

typedef struct {
    int p0;
    int p1;
    double inverse_matrix[4];
} Pair;

typedef struct {
    int p0;
    int p1;
    int p2;
    double inverse_matrix[9];
} Triplet;

typedef struct {
    double *cur_gains;
    double *prev_gains;
    double *temp_gains;
} Gains;

typedef struct {
    double *degs;
    double *rads;
    CartesianPoint *lpos;
    HullPoint *hull;
    Pair *lpairs;
    Gains *lgains;
    int n;
    int n_hull;
    int spread_points;
    CartesianPoint *spread_cloud;
    size_t n_samples;
    double *out_frame;
} VBAP2D;

typedef struct {
    PolarPoint *degs;
    PolarPoint *rads;
    CartesianPoint *lpos;
    Triplet *ltriplets;
    Gains *lgains;
    int n;
    int n_hull;
    int spread_points;
    double (*spread_angles)[2];
    CartesianPoint *spread_cloud;
    size_t n_samples;
    double *out_frame;
} VBAP3D;

typedef struct {
    double *distances;
    double *sorted_distances;
    double max_distance;
    double mean_distance;
} TempDistances;

typedef struct {
    PolarPoint *degs;
    PolarPoint *rads;
    CartesianPoint *lpos;
    Gains *lgains;
    double *weights;
    CartesianPoint center;
    TempDistances *temp_distances;
    double *temp_b;
    double *temp_u;
    double spatial_blur;
    double eta;
    double a;
    int n;
    size_t n_samples;
    double *out_frame;
} DBAP;

typedef struct {
    PANNIX_DIMENSIONS kind;
    union {
        VBAP2D *vbap2d;
        VBAP3D *vbap3d;
        DBAP *dbap;
    };
} PANNIX;

typedef struct {
    CartesianPoint start;
    CartesianPoint end;
} Position;

static inline double wrap_angle(double x, ANGLE_KIND akind);
static inline double clamp_elevation(double x, ANGLE_KIND akind);

// --- PANNIX (2D-3D) ---

/**
 * Allocates and initializes a new PANNIX instance.
 *
 * @param kind  The spatialization kind to use. Must be one of:
 *              - PVBAP2D : 2D Vector Base Amplitude Panning
 *              - PVBAP3D : 3D Vector Base Amplitude Panning
 *              - PDBAP   : Distance-Based Amplitude Panning
 *
 * @return Pointer to a newly allocated PANNIX structure on success.
 *         Returns NULL if allocation fails.
 *
 * @note The returned pointer must be freed using pannix_dealloc()
 *       to avoid memory leaks.
 *
 * @warning The internal structures of the PANNIX instance are
 *          not initialized beyond memory allocation; you must
 *          call the appropriate initialize_vbap() or initialize_dbap()
 *          before using the instance for audio processing.
 */
PANNIX* pannix_alloc(PANNIX_DIMENSIONS kind);

/**
 * Frees all memory associated with a PANNIX instance.
 *
 * @param p  Pointer to a PANNIX instance previously allocated with pannix_alloc().
 *
 * @note After calling this function, the pointer becomes invalid and must not be used.
 */
void pannix_dealloc(PANNIX *p);

/**
 * Initializes a VBAP (Vector Base Amplitude Panning) spatializer.
 *
 * @param pannix        Pointer to a PANNIX instance (must have kind = PVBAP2D or PVBAP3D).
 * @param frame_size    Audio frame size for processing.
 * @param pos_in_degs   Pointer to an array of loudspeaker positions in degrees (PolarPoint (PVBAP3D) or double[] (PVBAP2D)).
 * @param n             Number of loudspeakers (size of pos_in_degs array).
 * @param spread_points Optional spread parameter (0 for standard VBAP, >0 to diffuse the source over multiple directions).
 *
 * @return 0 on success, non-zero on error.
 *
 * @note The PANNIX instance must be allocated with pannix_alloc() before calling.
 */
int initialize_vbap(PANNIX *pannix, size_t frame_size, void *pos_in_degs, int n, int spread_points);

/**
 * Initializes a DBAP (Distance-Based Amplitude Panning) spatializer.
 *
 * @param pannix      Pointer to a PANNIX instance (must have kind = PDBAP).
 * @param frame_size  Audio frame size for processing.
 * @param pos_in_degs Pointer to an array of loudspeaker positions in degrees (PolarPoint or double[] depending on implementation).
 * @param n           Number of loudspeakers (size of pos_in_degs array).
 * @param rolloff     Rolloff exponent for distance attenuation.
 * @param weights     Optional array of loudspeaker weights (length = n), can be NULL for uniform weighting.
 *
 * @return 0 on success, non-zero on error.
 *
 * @note The PANNIX instance must be allocated with pannix_alloc() before calling.
 */
int initialize_dbap(PANNIX *pannix, size_t frame_size, void *pos_in_degs, int n, double rolloff, double *weights);

/**
 * Computes the gain vector for a given source position.
 *
 * @param pannix Pointer to an initialized PANNIX instance.
 * @param source Pointer to a CartesianPoint representing the source position.
 * @param spread Optional spread parameter for DBAP (smooth exponential attenuation).
 *
 * @note Updates internal gain buffers; call gain_vector_interpolation() to apply gains to audio frames.
 */
void solve_gain_vector(PANNIX *pannix, CartesianPoint *source, double spread);

/**
 * Applies gain vector interpolation to an audio frame.
 *
 * @param out          Pointer to the output frame buffer (length = nchnls * frame_size).
 * @param input_frame  Pointer to the input audio frame buffer (length = frame_size).
 * @param pannix       Pointer to an initialized PANNIX instance with solved gain vector.
 * @param ch_mode      Channel mode (PLANAR or INTERLEAVED).
 *
 * @note The output buffer will contain the spatialized frame.
 */
void gain_vector_interpolation(double *out, double *input_frame, PANNIX *pannix, CHANNEL_MODE ch_mode);

PolarPoint car_to_pol(const CartesianPoint *p);
CartesianPoint pol_to_car(const PolarPoint *p, ANGLE_KIND akind);
void generate_spread_cloud(CartesianPoint *cloud, PANNIX_DIMENSIONS pdim, CartesianPoint *source, double (*angles_table)[2], double spread, int n_points);
// ---

// --- PANNIX VBAP 2D ---
void convex_hull_2d(CartesianPoint *points, int n, HullPoint *hull, Pair *pairs, int *n_hull);
int initialize_vbap2d_from_loc(VBAP2D *vbap2d, size_t frame_size, double *pos_in_degs, int n, int spread_points);
int find_active_2d_arc(Pair *pair, VBAP2D *vbap2d, Position *source, double *angle);
void solve_vbap2d_internal_gain_vector(VBAP2D *vbap2d, CartesianPoint *source, double spread);
void calculate_matrix_2d_inv(Pair *pairs, CartesianPoint *pos, int n);
// ---

// --- PANNIX VBAP 3D ---
int initialize_vbap3d_from_loc(VBAP3D *vbap3d, size_t frame_size, PolarPoint *pos_in_degs, int n, int spread_points);
int find_active_3d_arc(Triplet *triplet, double *gains, VBAP3D *vbap3d, CartesianPoint *source, PolarPoint *angle);
void solve_vbap3d_internal_gain_vector(VBAP3D *vbap3d, CartesianPoint *source, double spread);
void calculate_matrix_3d_inv(Triplet *pairs, CartesianPoint *pos, int n);
// ---

// --- PANNIX DBAP 2D ---
int initialize_dbap_from_loc(DBAP *dbap, size_t frame_size, PolarPoint *pos_in_degs, int n, double rolloff, double *weights);
void solve_dbap_internal_gain_vector(DBAP *dbap, CartesianPoint *source, double spread); // spread means selectivity: 0 -> classic DBAP, spread > 0 more selective (smooth decay)
// ---

// --- INLINE IMPL ---
static inline double wrap_angle(double x, ANGLE_KIND akind) {
    double limit = (akind == RADIANS) ? (2 * PI) : 360.0;
    x = fmod(x, limit);
    if (x < 0.0) x += limit;
    return x;
}

static inline double clamp_elevation(double x, ANGLE_KIND akind) {
    double limit = (akind == RADIANS) ? (PI / 2.0) : 90.0;
    if (x > limit) return limit;
    if (x < -limit) return -limit;
    return x;
}
// --- END INLINE IMPL ---

#endif
