#ifndef CPNX_H
#define CPNX_H

#ifdef _WIN_32
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

#define NO_ERROR 0
#define VBAP_INIT_ERROR 1
#define VBAP2D_INIT_ERROR 2
#define VBAP3D_INIT_ERROR 3
#define DBAP_INIT_ERROR 4

#define PI (double) acos(-1.0)
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
PANNIX* pannix_alloc(PANNIX_DIMENSIONS kind);
void pannix_dealloc(PANNIX *p);
int initialize_vbap(PANNIX *pannix, size_t frame_size, void *pos_in_degs, int n, int spread_points);
int initialize_dbap(PANNIX *pannix, size_t frame_size, void *pos_in_degs, int n, double rolloff, double *weights);
void solve_gain_vector(PANNIX *pannix, CartesianPoint *source, double spread);
void gain_vector_interpolation(double *out, double *input_frame, PANNIX *pannix, CHANNEL_MODE ch_mode);

PolarPoint car_to_pol(const CartesianPoint *p);
CartesianPoint pol_to_car(const PolarPoint *p, ANGLE_KIND ANGLE_KIND);
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
