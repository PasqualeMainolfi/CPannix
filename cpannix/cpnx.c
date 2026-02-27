#include "cpnx.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CONVHULL_3D_ENABLE
#define CONVHULL_3D_IMPLEMENTATION
#include "convhull_3d.h"

#define SPREAD_ATTENUATION 0.5
#define DB_RATIO 6.02059991327962

static inline double get_distance(CartesianPoint a, CartesianPoint b, double spatial_blur) {
    spatial_blur = spatial_blur <= 0.0 ? 0.0 : spatial_blur * spatial_blur;
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    double dz = a.z - b.z;
    double dist = dx * dx + dy * dy + dz * dz;
    return sqrt(dist + spatial_blur);
}

static inline void swap(double *a, double *b) {
    double temp = *a;
    *a = *b;
    *b = temp;
}

static int32_t partition(double *arr, int left, int right) {
    int index = left;
    double pivot = arr[right];
    for (int i = left; i < right; i++) {
        if (arr[i] <= pivot) {
            swap(&arr[index], &arr[i]);
            index++;
        }
    }
    swap(&arr[index], &arr[right]);
    return index;
}

static double quickselect(double *arr, int n, int k) {
    int left = 0;
    int right = n - 1;

    while (left <= right) {
        if (left == right) return arr[left];

        int pivot_index = partition(arr, left, right);

        if (k == pivot_index) return arr[pivot_index];
        if (k > pivot_index) {
            left = pivot_index + 1;
        } else {
            right = pivot_index - 1;
        }
    }
    return arr[k];
}

PANNIX *pannix_alloc(PANNIX_DIMENSIONS kind) {
    PANNIX *pannix = malloc(sizeof(PANNIX));

    if (!pannix) {
        return NULL;
    };

    pannix->kind = kind;
    switch (kind) {
        case PVBAP2D:
            pannix->vbap2d = malloc(sizeof(VBAP2D));
            if (!pannix->vbap2d) {
                return NULL;
            }
            return pannix;
        case PVBAP3D:
            pannix->vbap3d = malloc(sizeof(VBAP3D));
            if (!pannix->vbap3d) {
                return NULL;
            }
            return pannix;
        case PDBAP:
            pannix->dbap = malloc(sizeof(DBAP));
            if (!pannix->dbap) {
                return NULL;
            }
            return pannix;
        default:
            return NULL;
    }
}

void pannix_dealloc(PANNIX *p) {
    if (!p) return;
    switch (p->kind) {
        case PVBAP2D:
            if (p->vbap2d){
                free(p->vbap2d->degs);
                free(p->vbap2d->rads);
                free(p->vbap2d->lpos);
                free(p->vbap2d->hull);
                free(p->vbap2d->lpairs);
                if (p->vbap2d->lgains) {
                    if (p->vbap2d->lgains->cur_gains) free(p->vbap2d->lgains->cur_gains);
                    if (p->vbap2d->lgains->prev_gains) free(p->vbap2d->lgains->prev_gains);
                    if (p->vbap2d->lgains->temp_gains) free(p->vbap2d->lgains->temp_gains);
                    free(p->vbap2d->lgains);
                }
                free(p->vbap2d->spread_cloud);
                if (p->vbap2d->out_frame) free(p->vbap2d->out_frame);
                free(p->vbap2d);
            }
            free(p);
            return;
        case PVBAP3D:
            if (p->vbap3d) {
                free(p->vbap3d->degs);
                free(p->vbap3d->rads);
                free(p->vbap3d->lpos);
                free(p->vbap3d->ltriplets);
                if (p->vbap3d->lgains) {
                    if (p->vbap3d->lgains->cur_gains) free(p->vbap3d->lgains->cur_gains);
                    if (p->vbap3d->lgains->prev_gains) free(p->vbap3d->lgains->prev_gains);
                    if (p->vbap3d->lgains->temp_gains) free(p->vbap3d->lgains->temp_gains);
                    free(p->vbap3d->lgains);
                }
                if (p->vbap3d->spread_angles) free(p->vbap3d->spread_angles);
                if (p->vbap3d->spread_cloud) free(p->vbap3d->spread_cloud);
                if (p->vbap3d->out_frame) free(p->vbap3d->out_frame);
                free(p->vbap3d);
            }
            free(p);
            return;
        case PDBAP:
            if (p->dbap){
                free(p->dbap->degs);
                free(p->dbap->rads);
                free(p->dbap->lpos);
                if (p->dbap->lgains) {
                    if (p->dbap->lgains->cur_gains) free(p->dbap->lgains->cur_gains);
                    if (p->dbap->lgains->prev_gains) free(p->dbap->lgains->prev_gains);
                    if (p->dbap->lgains->temp_gains) free(p->dbap->lgains->temp_gains);
                    free(p->dbap->lgains);
                }
                if (p->dbap->weights) free(p->dbap->weights);
                if (p->dbap->temp_distances) {
                    if (p->dbap->temp_distances->distances) free(p->dbap->temp_distances->distances);
                    if (p->dbap->temp_distances->sorted_distances) free(p->dbap->temp_distances->sorted_distances);
                    free(p->dbap->temp_distances);
                }
                if (p->dbap->temp_b) free(p->dbap->temp_b);
                if (p->dbap->temp_u) free(p->dbap->temp_u);
                if (p->dbap->out_frame) free(p->dbap->out_frame);
                free(p->dbap);
            }
            free(p);
            return;
        default:
            return;
    }
}

int cmp_points(const void *a, const void *b) {
    CartesianPoint *p = (CartesianPoint*) a;
    CartesianPoint *q = (CartesianPoint*) b;
    double a_angle = atan2(p->y, p->x);
    double b_angle = atan2(q->y, q->x);
    return a_angle < b_angle ? -1 : a_angle > b_angle;
}

double point_cross(CartesianPoint o, CartesianPoint a, CartesianPoint b) {
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
}

void convex_hull_2d(CartesianPoint *points, int n, HullPoint *hull, Pair *pairs, int *n_hull) {
    if (n < 3) {
        hull[0].p = points[0];
        hull[1].p = points[1];
        hull[0].index = 0;
        hull[1].index = 1;
        pairs[0].p0 = 0;
        pairs[0].p1 = 1;
        *n_hull = n;
        return;
    };

    CartesianPoint *temp_points = malloc(sizeof(CartesianPoint) * n);
    memcpy(temp_points, points, sizeof(CartesianPoint) * n);
    qsort(temp_points, n, sizeof(CartesianPoint), cmp_points);

    // LOWER HULL
    int k = 0;
    for (int i = 0; i < n; ++i) {
        while (k >= 2 && point_cross(hull[k - 2].p, hull[k - 1].p, temp_points[i]) <= 0) {
            k--;
        }
        int index = k++;
        hull[index].p = temp_points[i];
        hull[index].index = i;
    }

    // UPPER HULL
    int t = k + 1;
    for (int i = n - 2; i >= 0; --i) {
        while (k >= t && point_cross(hull[k - 2].p, hull[k - 1].p, temp_points[i]) <= 0) {
            k--;
        }
        int index = k++;
        hull[index].p = temp_points[i];
        hull[index].index = i;
    }

    free(temp_points);
    *n_hull = k - 1;

    for (int i = 0; i < *n_hull; i++) {
        pairs[i].p0 = hull[i].index;
        pairs[i].p1 = hull[(i + 1) % *n_hull].index;
    }
}

PolarPoint car_to_pol(const CartesianPoint *p) {
    PolarPoint pol;
    pol.rho = sqrt(p->x * p->x + p->y * p->y + p->z * p->z);
    pol.phi = atan2(p->y, p->x);
    pol.theta = asin(p->z / pol.rho);
    return pol;
};

CartesianPoint pol_to_car(const PolarPoint *p, ANGLE_KIND deg_mode) {
    CartesianPoint car;

    double phi = p->phi;
    double theta = p->theta;
    if (deg_mode == DEGREE) {
        phi = TO_RAD(wrap_angle(p->phi, DEGREE));
        theta = TO_RAD(clamp_elevation(p->theta, DEGREE));
    }

    car.x = p->rho * cos(theta) * cos(phi);
    car.y = p->rho * cos(theta) * sin(phi);
    car.z = p->rho * sin(theta);

    return car;
};

void set_2d_loudspeaker_position_from_num(double *angles, int n) {
    double offset = 180.0 / (double)n;
    double step = 360.0 / (double)n;
    for (int i = 0; i < n; i++) {
        angles[i] = wrap_angle(step * i + offset, DEGREE);
    }
};

void set_2d_loudspeaker_position(CartesianPoint *lpos, double *angles, int n) {
    for (int i = 0; i < n; i++) {
        lpos[i].x = cos(angles[i]);
        lpos[i].y = sin(angles[i]);
    }
};

void set_3d_loudspeaker_position(CartesianPoint *lpos, PolarPoint *angles, int n) {
    for (int i = 0; i < n; i++) {
        lpos[i] = pol_to_car(&angles[i], RADIANS);
    }
};

// --- VBAP ---

int init_vbap2d(VBAP2D *p, int n, int spread_points, size_t frame_size) {
    if (n < 2) {
        return VBAP2D_INIT_ERROR;
    }

    p->n = n;
    p->degs = malloc(sizeof(double) * n);
    p->rads = malloc(sizeof(double) * n);
    p->lpos = malloc(sizeof(CartesianPoint) * n);
    p->lpairs = malloc(sizeof(Pair) * 2 * n);
    p->hull = malloc(sizeof(HullPoint) * 2 * n);

    p->lgains = malloc(sizeof(Gains));
    if (!p->lgains) return VBAP2D_INIT_ERROR;

    p->lgains->cur_gains = calloc(n, sizeof(double));
    p->lgains->prev_gains = calloc(n, sizeof(double));
    p->lgains->temp_gains = calloc(n, sizeof(double));
    p->out_frame = calloc(frame_size * n, sizeof(double));
    p->n_hull = 0;
    p->spread_points = spread_points;
    p->spread_cloud = spread_points > 0 ? malloc(sizeof(CartesianPoint) * spread_points) : NULL;

    if (
        !p->degs ||
        !p->rads ||
        !p->lpos ||
        !p->hull ||
        !p->lgains->cur_gains  ||
        !p->lgains->prev_gains ||
        !p->lgains->temp_gains ||
        !p->out_frame
    ) {
        return VBAP2D_INIT_ERROR;
    }

    return NO_ERROR;
}

int init_vbap3d(VBAP3D *p, int n, int spread_points, size_t frame_size) {
    if (n < 3) {
        return VBAP3D_INIT_ERROR;
    }

    p->n = n;
    p->degs = malloc(sizeof(PolarPoint) * n);
    p->rads = malloc(sizeof(PolarPoint) * n);
    p->lpos = malloc(sizeof(CartesianPoint) * n);
    p->ltriplets = NULL;

    p->lgains = malloc(sizeof(Gains));
    if (!p->lgains) return VBAP3D_INIT_ERROR;

    p->lgains->cur_gains = calloc(n, sizeof(double));
    p->lgains->prev_gains = calloc(n, sizeof(double));
    p->lgains->temp_gains = calloc(n, sizeof(double));
    p->out_frame = calloc(frame_size * n, sizeof(double));
    p->n_hull = 0;
    p->spread_points = spread_points;
    p->spread_angles = NULL;
    p->spread_cloud = NULL;

    if (
        !p->degs ||
        !p->rads ||
        !p->lpos ||
        !p->lgains->cur_gains  ||
        !p->lgains->prev_gains ||
        !p->lgains->temp_gains ||
        !p->out_frame
    ) {
        return VBAP3D_INIT_ERROR;
    }

    if (spread_points > 0) {
        p->spread_cloud = malloc(sizeof(CartesianPoint) * spread_points);
        p->spread_angles = malloc(sizeof(double[2]) * spread_points);
        for (int i = 0; i < spread_points; i++) {
            double a = (2 * PI * i) / spread_points;
            p->spread_angles[i][0] = cos(a);
            p->spread_angles[i][1] = sin(a);
        }
    }
    return NO_ERROR;
}

void complete_vbap2d(VBAP2D *p) {
    for (int i = 0; i < p->n; i++) {
        p->rads[i] = TO_RAD(wrap_angle(p->degs[i], DEGREE));
    }

    set_2d_loudspeaker_position(p->lpos, p->rads, p->n);
    convex_hull_2d(p->lpos, p->n, p->hull, p->lpairs, &p->n_hull);
    calculate_matrix_2d_inv(p->lpairs, p->lpos, (*p).n_hull);
}

void fix_triplets_orientation(VBAP3D *p) {
    for (int k = 0; k < p->n_hull; k++) {
        int i0 = p->ltriplets[k].p0;
        int i1 = p->ltriplets[k].p1;
        int i2 = p->ltriplets[k].p2;

        CartesianPoint v0 = p->lpos[i0];
        CartesianPoint v1 = p->lpos[i1];
        CartesianPoint v2 = p->lpos[i2];

        double det =
            v0.x * (v1.y * v2.z - v1.z * v2.y) -
            v0.y * (v1.x * v2.z - v1.z * v2.x) +
            v0.z * (v1.x * v2.y - v1.y * v2.x);

        if (det < 0) {
            int tmp = p->ltriplets[k].p1;
            p->ltriplets[k].p1 = p->ltriplets[k].p2;
            p->ltriplets[k].p2 = tmp;
        }
    }
}

int complete_vbap3d(VBAP3D *p) {
    for (int i = 0; i < p->n; i++) {
        p->rads[i].phi = TO_RAD(wrap_angle(p->degs[i].phi, DEGREE));
        p->rads[i].theta = TO_RAD(clamp_elevation(p->degs[i].theta, DEGREE));
        p->rads[i].rho = p->degs[i].rho;
    }

    set_3d_loudspeaker_position(p->lpos, p->rads, p->n);

    ch_vertex *vertices = malloc(sizeof(ch_vertex) * p->n);
    for (int i = 0; i < p->n; i++) {
        vertices[i].x = p->lpos[i].x;
        vertices[i].y = p->lpos[i].y;
        vertices[i].z = p->lpos[i].z;
    };

    int *faces = NULL;
    convhull_3d_build(vertices, p->n, &faces, &p->n_hull);

    p->ltriplets = malloc(sizeof(Triplet) * p->n_hull);
    if (!p->ltriplets) {
        return VBAP3D_INIT_ERROR;
    }

    for (int i = 0; i < p->n_hull; i++) {
        p->ltriplets[i].p0 = faces[i * 3];
        p->ltriplets[i].p1 = faces[i * 3 + 1];
        p->ltriplets[i].p2 = faces[i * 3 + 2];
    }
    free(faces);

    fix_triplets_orientation(p);
    calculate_matrix_3d_inv(p->ltriplets, p->lpos, (*p).n_hull);
    return NO_ERROR;
}

void calculate_matrix_2d_inv(Pair *p, CartesianPoint *pos, int n) {
    for (int i = 0; i < n; i++) {
        int p0 = p[i].p0;
        int p1 = p[i].p1;

        double a = pos[p0].x;
        double b = pos[p1].x;
        double c = pos[p0].y;
        double d = pos[p1].y;

        double det = a * d - b * c;
        if (fabs(det) > 1e-15) {
            double inv_det = 1.0 / det;
            p[i].inverse_matrix[0] =  d * inv_det;
            p[i].inverse_matrix[1] = -b * inv_det;
            p[i].inverse_matrix[2] = -c * inv_det;
            p[i].inverse_matrix[3] =  a * inv_det;
        }
    }
}

void calculate_matrix_3d_inv(Triplet *p, CartesianPoint *pos, int n) {
    for (int k = 0; k < n; k++) {
        CartesianPoint a = pos[p[k].p0];
        CartesianPoint b = pos[p[k].p1];
        CartesianPoint c = pos[p[k].p2];

        double cp_x = b.y * c.z - b.z * c.y;
        double cp_y = b.z * c.x - b.x * c.z;
        double cp_z = b.x * c.y - b.y * c.x;
        double det = a.x * cp_x + a.y * cp_y + a.z * cp_z;

        if (fabs(det) > 1e-12) {
            double inv_det = 1.0 / det;
            p[k].inverse_matrix[0] = cp_x * inv_det;
            p[k].inverse_matrix[1] = cp_y * inv_det;
            p[k].inverse_matrix[2] = cp_z * inv_det;
            p[k].inverse_matrix[3] = (c.y * a.z - c.z * a.y) * inv_det;
            p[k].inverse_matrix[4] = (c.z * a.x - c.x * a.z) * inv_det;
            p[k].inverse_matrix[5] = (c.x * a.y - c.y * a.x) * inv_det;
            p[k].inverse_matrix[6] = (a.y * b.z - a.z * b.y) * inv_det;
            p[k].inverse_matrix[7] = (a.z * b.x - a.x * b.z) * inv_det;
            p[k].inverse_matrix[8] = (a.x * b.y - a.y * b.x) * inv_det;
        }
    }
}

void cross_matrix2d_prod(double *res, const double *inv_mat, const CartesianPoint *source) {
    res[0] = inv_mat[0] * source->x + inv_mat[1] * source->y;
    res[1] = inv_mat[2] * source->x + inv_mat[3] * source->y;
}

void cross_matrix3d_prod(double *res, const double *inv_mat, const CartesianPoint *source) {
    res[0] = inv_mat[0] * source->x + inv_mat[1] * source->y + inv_mat[2] * source->z;
    res[1] = inv_mat[3] * source->x + inv_mat[4] * source->y + inv_mat[5] * source->z;
    res[2] = inv_mat[6] * source->x + inv_mat[7] * source->y + inv_mat[8] * source->z;
}

int initialize_vbap2d_from_loc(VBAP2D *vbap2d, size_t frame_size, double *pos_in_degs, int n, int spread_points) {
    int err = init_vbap2d(vbap2d, n, spread_points, frame_size);
    if (err != NO_ERROR) {
        return err;
    }
    memcpy(vbap2d->degs, pos_in_degs, sizeof(double) * n);
    complete_vbap2d(vbap2d);
    return NO_ERROR;
}

int initialize_vbap3d_from_loc(VBAP3D *vbap3d, size_t frame_size, PolarPoint *pos_in_degs, int n, int spread_points) {
    int err = init_vbap3d(vbap3d, n, spread_points, frame_size);
    if (err != NO_ERROR) {
        return err;
    }
    memcpy(vbap3d->degs, pos_in_degs, sizeof(PolarPoint) * n);
    err = complete_vbap3d(vbap3d);
    if (err != NO_ERROR) {
        return err;
    }
    return NO_ERROR;
}

int get_intersection(Position *source, Position *lbase) {
    double sp1x = source->start.x;
    double sp1y = source->start.y;
    double sp2x = source->end.x;
    double sp2y = source->end.y;

    double bp1x = lbase->start.x;
    double bp1y = lbase->start.y;
    double bp2x = lbase->end.x;
    double bp2y = lbase->end.y;

    double den = (bp1x - bp2x) * (sp1y - sp2y) - (bp1y - bp2y) * (sp1x - sp2x);

    if (den == 0.0) return 0;

    double t = ((bp1x - sp1x) * (sp1y - sp2y) - (bp1y - sp1y) * (sp1x - sp2x)) / den;
    double u = ((bp1x - sp1x) *(bp1y - bp2y) -(bp1y - sp1y) *(bp1x - bp2x)) / den;

    if ((t >= 0.0 && t <= 1.0) && (u > 0.0)) {
        // double xi = base_p1x + t * (base_p2x - base_p1x);
        // double yi = base_p1y + t * (base_p2y - base_p1y);
        return 1;
    }

    return 0;
}

int find_active_2d_arc(Pair *pair, VBAP2D *vbap2d, Position *source, double *angle) {
    for (int i = 0; i < vbap2d->n; i++) {
        if (vbap2d->rads[i] == *angle) {
            pair->p0 = i;
            pair->p1 = -1;
            return 1;
        }
    }

    Position p;
    for (int i = 0; i < vbap2d->n_hull; i++) {
        p.start = vbap2d->lpos[vbap2d->lpairs[i].p0];
        p.end = vbap2d->lpos[vbap2d->lpairs[i].p1];
        int intersect = get_intersection(source, &p);
        if (intersect) {
            *pair = vbap2d->lpairs[i];
            return 1;
        }
    }
    return 0;
}

int find_active_3d_arc(Triplet *triplet, double *gains, VBAP3D *vbap3d, CartesianPoint *source, PolarPoint *angle) {
    for (int i = 0; i < vbap3d->n; i++) {
        if (vbap3d->rads[i].phi == angle->phi && vbap3d->rads[i].theta == angle->theta) {
            triplet->p0 = i;
            triplet->p1 = -1;
            triplet->p2 = -1;
            return 1;
        }
    }

    double best_max_min_gain = -1e15;
    int best_triplets = -1;
    double best_res[3];

    for (int i = 0; i < vbap3d->n_hull; i++) {
        Triplet t = vbap3d->ltriplets[i];
        double res[3];
        cross_matrix3d_prod(res, t.inverse_matrix, source);

        double min_gain = res[0];
        if (res[1] < min_gain) min_gain = res[1];
        if (res[2] < min_gain) min_gain = res[2];

        if (min_gain > -0.0001) {
            *triplet = vbap3d->ltriplets[i];
            memcpy(gains, res, sizeof(double) * 3);
            return 1;
        }

        if (min_gain > best_max_min_gain) {
            best_max_min_gain = min_gain;
            best_triplets = i;
            memcpy(best_res, res, sizeof(double) * 3);
        }

    }

    if (best_triplets != -1) {
        *triplet = vbap3d->ltriplets[best_triplets];
        memcpy(gains, best_res, sizeof(double) * 3);

        if (gains[0] < 0.0) gains[0] = 0.0;
        if (gains[1] < 0.0) gains[1] = 0.0;
        if (gains[2] < 0.0) gains[2] = 0.0;
        return 1;
    }

    return 0;
}

void _compute_vbap2d_loudspeaker_gains(double *current_gains, VBAP2D *vbap2d, CartesianPoint *source) {
    memset(current_gains, 0, sizeof(double) * vbap2d->n);

    Position source_position;
    source_position.start.x = 0.0;
    source_position.start.y = 0.0;
    source_position.end.x = source->x;
    source_position.end.y = source->y;
    double source_angle = atan2(source->y, source->x);

    int is_dim = 0;
    Position local_base;
    Pair local_pair;

    int has_arc = find_active_2d_arc(&local_pair, vbap2d, &source_position, &source_angle);

    if (has_arc) {
        int start_index = local_pair.p0;
        int end_index = local_pair.p1;
        local_base.start.x = vbap2d->lpos[start_index].x;
        local_base.start.y = vbap2d->lpos[start_index].y;

        if (end_index != -1) {
            local_base.end.x = vbap2d->lpos[local_pair.p1].x;
            local_base.end.y = vbap2d->lpos[local_pair.p1].y;
            is_dim = 1;
        }
    }

    if (is_dim) {
        double gains[2];
        cross_matrix2d_prod(gains, local_pair.inverse_matrix, source);
        current_gains[local_pair.p0] = gains[0];
        current_gains[local_pair.p1] = gains[1];
    } else {
        current_gains[local_pair.p0] = 1.0;
    }
}

void solve_vbap_gain_vector_helper(
    void *vbap,
    PANNIX_DIMENSIONS kind,
    CartesianPoint *spread_cloud,
    int spread_points,
    int n,
    Gains *lgains,
    CartesianPoint *source,
    double spread,
    double (*spread_angles)[2],
    void (*get_gain_vector2d)(double*, VBAP2D*, CartesianPoint*),
    void (*get_gain_vector3d)(double*, VBAP3D*, CartesianPoint*)
) {

    if (spread > 0.0 && spread_points > 0 && spread_cloud) {
        generate_spread_cloud(spread_cloud, kind, source, spread_angles, spread, spread_points);

        double weight = SPREAD_ATTENUATION / spread_points;
        for (int i = 0; i < spread_points; i++) {

            if (get_gain_vector2d) {
                get_gain_vector2d(lgains->temp_gains, (VBAP2D *)vbap, &spread_cloud[i]);
            }
            if (get_gain_vector3d) {
                get_gain_vector3d(lgains->temp_gains, (VBAP3D *)vbap, &spread_cloud[i]);
            }

            for (int j = 0; j < n; j++) {
                lgains->cur_gains[j] += lgains->temp_gains[j] * weight;
            }
        }
    }

    double norm = 0.0;
    for (int i = 0; i < n; i++) {
        double g = lgains->cur_gains[i];
        norm += g * g;
    }

    if (norm > 0.0) {
        norm = sqrt(norm);
        for (int i = 0; i < n; i++) {
            lgains->cur_gains[i] /= norm;
        }
    } else {
        memset(lgains->cur_gains, 0, sizeof(double) * n);
    }
}

void solve_vbap2d_internal_gain_vector(VBAP2D *vbap2d, CartesianPoint *source, double spread) {
    _compute_vbap2d_loudspeaker_gains(vbap2d->lgains->cur_gains, vbap2d, source);
    solve_vbap_gain_vector_helper(
        vbap2d,
        PVBAP2D,
        vbap2d->spread_cloud,
        vbap2d->spread_points,
        vbap2d->n,
        vbap2d->lgains,
        source,
        spread,
        NULL,
        _compute_vbap2d_loudspeaker_gains,
        NULL
    );
}

void _compute_vbap3d_loudspeaker_gains(double *current_gains, VBAP3D *vbap3d, CartesianPoint *source) {
    memset(current_gains, 0, sizeof(double) * vbap3d->n);

    PolarPoint polar_source = car_to_pol(source);

    Triplet local_triplet;
    double gains[3];
    int has_arc = find_active_3d_arc(&local_triplet, gains, vbap3d, source, &polar_source);

    if (has_arc) {
        int is_dim = local_triplet.p1 != -1 && local_triplet.p2 != -1;
        if (is_dim) {
            current_gains[local_triplet.p0] = gains[0];
            current_gains[local_triplet.p1] = gains[1];
            current_gains[local_triplet.p2] = gains[2];
        } else {
            current_gains[local_triplet.p0] = 1.0;
        }
    }
}

void solve_vbap3d_internal_gain_vector(VBAP3D *vbap3d, CartesianPoint *source, double spread) {
    _compute_vbap3d_loudspeaker_gains(vbap3d->lgains->cur_gains, vbap3d, source);
    solve_vbap_gain_vector_helper(
        vbap3d,
        PVBAP3D,
        vbap3d->spread_cloud,
        vbap3d->spread_points,
        vbap3d->n,
        vbap3d->lgains,
        source,
        spread,
        vbap3d->spread_angles,
        NULL,
        _compute_vbap3d_loudspeaker_gains
    );
}

int initialize_vbap(PANNIX *pannix, size_t frame_size, void *pos_in_degs, int n, int spread_points) {
    int err = NO_ERROR;
    switch (pannix->kind) {
        case PVBAP2D:
            err = initialize_vbap2d_from_loc(pannix->vbap2d, frame_size, (double *)pos_in_degs, n, spread_points);
            break;
        case PVBAP3D:
            err = initialize_vbap3d_from_loc(pannix->vbap3d, frame_size, (PolarPoint *)pos_in_degs, n, spread_points);
            break;
        default:
            return VBAP_INIT_ERROR;
    }
    return err;
}

void solve_gain_vector(PANNIX *pannix, CartesianPoint *source, double spread) {
    switch (pannix->kind) {
        case PVBAP2D:
            solve_vbap2d_internal_gain_vector(pannix->vbap2d, source, spread);
            return;
        case PVBAP3D:
            solve_vbap3d_internal_gain_vector(pannix->vbap3d, source, spread);
            return;
        case PDBAP:
            solve_dbap_internal_gain_vector(pannix->dbap, source, spread);
            return;
        default:
            return;
    }
}

double get_spread_value(PANNIX_DIMENSIONS pdim, double spread) {
    if (spread <= 0.0) return 0.0;
    if (spread >= 1.0) return pdim == PVBAP2D ? 360.0 : 180.0;
    if (pdim == PVBAP2D) return spread * 360.0;

    double angle = acos(1.0 - 2.0 * spread);
    return TO_DEG(angle);
}

CartesianPoint cart_cross_prod(const CartesianPoint *a, const CartesianPoint *b) {
    CartesianPoint res;
    res.x = a->y * b->z - a->z * b->y;
    res.y = a->z * b->x - a->x * b->z;
    res.z = a->x * b->y - a->y * b->x;
    return res;
}

CartesianPoint get_orthogonal(const CartesianPoint *p) {
    if (fabs(p->x) < 0.707) {
        CartesianPoint vx = { 1.0, 0.0, 0.0 };
        return cart_cross_prod(p, &vx);
    } else {
        CartesianPoint vz = { 0.0, 0.0, 1.0 };
        return cart_cross_prod(p, &vz);
    }
}

CartesianPoint normalize(const CartesianPoint *p) {
    CartesianPoint normalized;
    double norm = sqrt(p->x * p->x + p->y * p->y + p->z * p->z);

    if (norm > 0) {
        normalized.x = p->x / norm;
        normalized.y = p->y / norm;
        normalized.z = p->z / norm;
    } else {
        normalized.x = 1.0;
    }
    return normalized;
}

void generate_spread_cloud(CartesianPoint *cloud, PANNIX_DIMENSIONS pdim, CartesianPoint *source, double (*angles_table)[2], double spread, int n_points) {
    double spread_deg = get_spread_value(pdim, spread);
    double spread_rad = TO_RAD(spread_deg / 2.0);

    if (pdim == PVBAP2D) {
        PolarPoint source_polar = car_to_pol(source);
        double offset_step = (n_points > 1) ? (2.0 * spread_rad) / (n_points - 1) : 0;
        for (int i = 0; i < n_points; i++) {
            double alpha_offset = -spread_rad + (i * offset_step);
            PolarPoint virtual = { 1.0, source_polar.phi + alpha_offset, 0.0 };
            cloud[i] = pol_to_car(&virtual, RADIANS);
        }
    } else {
        CartesianPoint ortho = get_orthogonal(source);
        CartesianPoint right = normalize(&ortho);

        CartesianPoint cp = cart_cross_prod(source, &right);
        CartesianPoint up = normalize(&cp);

        for (int i = 0; i < n_points; i++) {

            CartesianPoint offset_dir;
            offset_dir.x = right.x * angles_table[i][0] + up.x * angles_table[i][1];
            offset_dir.y = right.y * angles_table[i][0] + up.y * angles_table[i][1];
            offset_dir.z = right.z * angles_table[i][0] + up.z * angles_table[i][1];

            cloud[i].x = source->x * cos(spread_rad) + offset_dir.x * sin(spread_rad);
            cloud[i].y = source->y * cos(spread_rad) + offset_dir.y * sin(spread_rad);
            cloud[i].z = source->z * cos(spread_rad) + offset_dir.z * sin(spread_rad);

            cloud[i] = normalize(&cloud[i]);
        }
    }
}
// --- END VBAP ---

// --- DBAP ---

int init_dbap(DBAP *dbap, int n, double rolloff, double *weights, size_t frame_size) {
    if (n < 2) {
        return DBAP_INIT_ERROR;
    }

    dbap->n = n;
    dbap->degs = malloc(sizeof(PolarPoint) * n);
    dbap->rads = malloc(sizeof(PolarPoint) * n);
    dbap->lpos = malloc(sizeof(CartesianPoint) * n);

    dbap->lgains = malloc(sizeof(Gains));
    dbap->temp_distances = malloc(sizeof(TempDistances));
    if (!dbap->lgains || !dbap->temp_distances) return DBAP_INIT_ERROR;

    dbap->lgains->cur_gains = calloc(n, sizeof(double));
    dbap->lgains->prev_gains = calloc(n, sizeof(double));
    dbap->lgains->temp_gains = calloc(n, sizeof(double));
    dbap->temp_distances->distances = calloc(n, sizeof(double));
    dbap->temp_distances->sorted_distances = calloc(n, sizeof(double));
    dbap->out_frame = calloc(n * frame_size, sizeof(double));
    dbap->temp_b = calloc(n, sizeof(double));
    dbap->temp_u = calloc(n, sizeof(double));
    dbap->weights = malloc(sizeof(double) * n);
    dbap->a = rolloff / DB_RATIO;
    dbap->spatial_blur = 0.0;
    dbap->eta = 0.0;

    if (
        !dbap->degs ||
        !dbap->rads ||
        !dbap->lpos ||
        !dbap->lgains->cur_gains  ||
        !dbap->lgains->prev_gains ||
        !dbap->lgains->temp_gains ||
        !dbap->temp_distances->distances        ||
        !dbap->temp_distances->sorted_distances ||
        !dbap->temp_b ||
        !dbap->temp_u ||
        !dbap->out_frame
    ) {
        return DBAP_INIT_ERROR;
    }

    if (weights) {
        memcpy(dbap->weights, weights, sizeof(double) * n);
    } else {
        for (int i = 0; i < n; i++)  dbap->weights[i] = 1.0;
    }
    return NO_ERROR;
}

void update_loudspeaker_distances(DBAP *dbap, const CartesianPoint *source, double spatial_blur) {
    dbap->temp_distances->max_distance = 0.0;
    dbap->temp_distances->mean_distance = 0.0;
    double max_dist = 0.0;
    double mean_dist = 0.0;
    for (int i = 0; i < dbap->n; i++) {
        double d = get_distance(dbap->lpos[i], *source, spatial_blur);
        dbap->temp_distances->distances[i] = d;
        max_dist = d > max_dist ? d : max_dist;
        mean_dist += d;
    }
    mean_dist /= (double)dbap->n;
    dbap->temp_distances->max_distance = max_dist;
    dbap->temp_distances->mean_distance = mean_dist;
}

double get_spatial_blur(DBAP *dbap) {
    update_loudspeaker_distances(dbap, &dbap->center, 0.0);
    double sb = 0.0;
    for (int i = 0; i < dbap->n; i++) {
        sb += dbap->temp_distances->distances[i];
    }
    return (sb / dbap->n) + 0.2;
}

void set_loudspeaker_position(DBAP *dbap, PolarPoint *input_coords) {
    for (int i = 0; i < dbap->n; i++) {
        dbap->lpos[i] = pol_to_car(&input_coords[i], DEGREE);
    }
};

void complete_dbap(DBAP *dbap) {
    for (int i = 0; i < dbap->n; i++) {
        dbap->rads[i].rho = dbap->degs[i].rho;
        dbap->rads[i].phi = TO_RAD(wrap_angle(dbap->degs[i].phi, DEGREE));
        dbap->rads[i].theta = TO_RAD(clamp_elevation(dbap->degs[i].theta, DEGREE));
    }

    set_3d_loudspeaker_position(dbap->lpos, dbap->rads, dbap->n);

    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    for (int i = 0; i < dbap->n; i++) {
        x += dbap->lpos[i].x;
        y += dbap->lpos[i].y;
        z += dbap->lpos[i].z;
    }

    dbap->center.x = x / dbap->n;
    dbap->center.y = y / dbap->n;
    dbap->center.z = z / dbap->n;

    dbap->spatial_blur = get_spatial_blur(dbap);
    dbap->eta = dbap->spatial_blur / dbap->n;
    memset(dbap->temp_distances->distances, 0, sizeof(double) * dbap->n);
}

int initialize_dbap_from_loc(DBAP *dbap, size_t frame_size, PolarPoint *pos_in_degs, int n, double rolloff, double *weights) {
    int err = init_dbap(dbap, n, rolloff, weights, frame_size);
    if (err != NO_ERROR) {
        return err;
    }
    memcpy(dbap->degs, pos_in_degs, sizeof(PolarPoint) * n);
    complete_dbap(dbap);
    return NO_ERROR;
}

double get_p(DBAP *dbap, const CartesianPoint *source) {
    double dist_source_center = get_distance(dbap->center, *source, dbap->spatial_blur);
    dist_source_center = dist_source_center > 0 ? dist_source_center : 1;
    double q = dbap->temp_distances->max_distance / dist_source_center;
    double p = q < 1 ? q : 1.0;
    return p;
}

void get_b(DBAP *dbap, double p) {
    memset(dbap->temp_b, 0, sizeof(double) * dbap->n);
    memset(dbap->temp_u, 0, sizeof(double) * dbap->n);
    double *u = dbap->temp_u;

    double unorm = 0.0;
    for (int i = 0; i < dbap->n; i++) {
        u[i] = dbap->temp_distances->distances[i] - dbap->temp_distances->max_distance;
        unorm += u[i] * u[i];
    }

    unorm = unorm > 0.001 ? unorm : 1.0;
    for (int i = 0; i < dbap->n; i++) {
        double v = u[i] / unorm;
        u[i] = (v * v) + dbap->eta;
    }

    memcpy(dbap->temp_distances->sorted_distances, dbap->temp_distances->distances, sizeof(double) * dbap->n);
    double upper_mid = quickselect(dbap->temp_distances->sorted_distances, dbap->n, dbap->n / 2);
    double *sorted_distances = dbap->temp_distances->sorted_distances;

    double median;
    if (dbap->n % 2 == 1) {
        median = upper_mid;
    } else {
        double lower_mid = sorted_distances[0];
        for (int32_t i = 0; i < dbap->n / 2; i++) {
            if (sorted_distances[i] > lower_mid) lower_mid = sorted_distances[i];
        }
        median = (lower_mid + upper_mid) * 0.5;
    }

    double *b_values = dbap->temp_b;
    double fac = (1.0 / p) + 1.0;
    for (int i = 0; i < dbap->n; i++) {
        double b = u[i] / (median > 0.0000001 ? median : 0.0000001) * fac;
        b_values[i] = (b * b) + 1.0;
    }
}

double get_k(DBAP *dbap, double p) {
    double pfac = 2 * dbap->a;
    double k_num = pow(p, pfac);
    double k_den = 0.0;
    for (int i = 0; i < dbap->n; i++) {
        double b = dbap->temp_b[i];
        double w = dbap->weights[i];
        double d = dbap->temp_distances->distances[i];
        k_den += (b * b * w * w) / pow(d, pfac);
    }

    double k = k_num / sqrt(k_den);
    return k;
}

void solve_dbap_internal_gain_vector(DBAP *dbap, CartesianPoint *source, double spread) {
    memset(dbap->lgains->cur_gains, 0, sizeof(double) * dbap->n);

    update_loudspeaker_distances(dbap, source, dbap->spatial_blur);
    double p = get_p(dbap, source);
    get_b(dbap, p);
    double k = get_k(dbap, p);
    double mean_distance = dbap->temp_distances->mean_distance;
    double beta = spread / (mean_distance + 0.000001); // soft-limited decay
    double sum_sq = 0.0;
    for (int i = 0; i < dbap->n; i++) {
        double b = dbap->temp_b[i];
        double w = dbap->weights[i];
        double d = dbap->temp_distances->distances[i];
        double gain = (k * w * b) / pow(d, dbap->a);
        gain *= exp(-beta * d);
        dbap->lgains->cur_gains[i] = gain;
        sum_sq += gain * gain;
    }

    double norm = sqrt(sum_sq);
    norm = (norm < 0.000001) ? 1.0 : norm;
    for (int32_t i = 0; i < dbap->n; i++) {
        dbap->lgains->cur_gains[i] /= norm;
    }

}

int initialize_dbap(PANNIX *pannix, size_t frame_size, void *pos_in_degs, int n, double rolloff, double *weights) {
    int err = NO_ERROR;
    switch (pannix->kind) {
        case PDBAP:
            err = initialize_dbap_from_loc(pannix->dbap, frame_size, (PolarPoint *)pos_in_degs, n, rolloff, weights);
            break;
        default:
            return DBAP_INIT_ERROR;
    }
    return err;
}

// --- GENERATE AUDIO OUT ---

void gain_vector_interpolation(double *out, double *input_frame, PANNIX *pannix, CHANNEL_MODE ch_mode) {
    Gains *gains;
    double *out_frame;
    size_t nchnls = 0;
    size_t n_samples = 0;
    switch (pannix->kind) {
        case PVBAP2D:
            gains = pannix->vbap2d->lgains;
            out_frame = pannix->vbap2d->out_frame;
            nchnls = (size_t)pannix->vbap2d->n;
            n_samples = pannix->vbap2d->n_samples;
            break;
        case PVBAP3D:
            gains = pannix->vbap3d->lgains;
            out_frame = pannix->vbap3d->out_frame;
            nchnls = (size_t)pannix->vbap3d->n;
            n_samples = pannix->vbap3d->n_samples;
            break;
        case PDBAP:
            gains = pannix->dbap->lgains;
            out_frame = pannix->dbap->out_frame;
            nchnls = (size_t)pannix->dbap->n;
            n_samples = pannix->dbap->n_samples;
            break;
        default:
            break;
    }

    if (!gains || !out_frame) {
        return;
    }

    for (size_t i = 0; i < nchnls; i++) {
        double g0 = gains->prev_gains[i];
        double g1 = gains->cur_gains[i];

        if (g0 == 0.0 && g1 == 0.0) {
            if (ch_mode == PLANAR) {
                size_t index = (size_t)i *  n_samples;
                memset(out_frame + index, 0, sizeof(double) * n_samples);
            } else {
                for (size_t j = 0; j < n_samples; j++) {
                    size_t index = j * nchnls + i;
                    out_frame[index] = 0.0;
                }
            }
            continue;
        }

        double gdiff = (g1 - g0) / n_samples;
        double g = g0;
        for (size_t j = 0; j < n_samples; j++) {
            size_t index = (ch_mode == PLANAR) ? i * n_samples + j : j * nchnls + i;
            out_frame[index] = input_frame[j] * g;
            g += gdiff;
        }
    }

    memcpy(gains->prev_gains, gains->cur_gains, sizeof(double) * nchnls);
    memcpy(out, out_frame, sizeof(double) * nchnls * n_samples);
}
