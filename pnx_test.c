#include <stdio.h>
#include <stdlib.h>
#include "cpannix/cpnx.h"

#define SPREAD 0.3
#define SPREAD_POINTS 5
#define DBAP_SPREAD 3
#define N_2D_LOUDSPEAKERS 6
#define N_3D_LOUDSPEAKERS 8
#define AUDIO_BLOCK_SIZE 2048

int main() {
    double degs2d[N_2D_LOUDSPEAKERS] = { -30.0, 30.0, 90.0, 150.0, 210.0, 270.0 };
    PolarPoint degs3d[N_3D_LOUDSPEAKERS] = {
        { 1.0,   30.0,  30.0 },
        { 1.0,  150.0,  30.0 },
        { 1.0, -150.0,  30.0 },
        { 1.0,  -30.0,  30.0 },
        { 1.0,   60.0, -30.0 },
        { 1.0,  120.0, -30.0 },
        { 1.0, -120.0, -30.0 },
        { 1.0,  -60.0, -30.0 }
    };

    PolarPoint dbap_degs[N_3D_LOUDSPEAKERS] = {
        { 1.0,   0.0, 0.0 },
        { 1.0,  45.0, 0.0 },
        { 1.0,  90.0, 0.0 },
        { 1.0, 135.0, 0.0 },
        { 1.0, 180.0, 0.0 },
        { 1.0, 225.0, 0.0 },
        { 1.0, 270.0, 0.0 },
        { 1.0, 315.0, 0.0 }
    };

    double source_phi = 90 * PI / 180.0;
    PolarPoint source = { .rho = 1.0, .phi = source_phi, .theta = 0.0 };
    CartesianPoint source_cartesian = pol_to_car(&source, RADIANS);
    printf("Source: (%f, %f, %f)\n", source_cartesian.x, source_cartesian.y, source_cartesian.z);

    PANNIX *vbap2d = pannix_alloc(PVBAP2D);
    int err = initialize_vbap(vbap2d, AUDIO_BLOCK_SIZE, degs2d, N_2D_LOUDSPEAKERS, SPREAD_POINTS);
    if (err != NO_ERROR) {
        printf("[ERROR] Init VBAP2D error!\n");
        return 1;
    }

    printf("[VBAP2D INFO]\n");
    for (int i = 0; i < vbap2d->vbap2d->n_hull; ++i) {
        printf("Hull point: (%f, %f)\n", vbap2d->vbap2d->hull[i].p.x, vbap2d->vbap2d->hull[i].p.y);
        printf("Pairs: (%d, %d)\n", vbap2d->vbap2d->lpairs[i].p0, vbap2d->vbap2d->lpairs[i].p1);
    }

    solve_gain_vector(vbap2d, &source_cartesian, SPREAD);

    for (int i = 0; i < vbap2d->vbap2d->n; ++i) {
        printf("GAINS: %f\n", vbap2d->vbap2d->lgains->cur_gains[i]);
    }

    PolarPoint polar_source3d = { 1.0, -45.0, 10.0 };
    CartesianPoint cartesian_source3d = pol_to_car(&polar_source3d, DEGREE);

    PANNIX *vbap3d = pannix_alloc(PVBAP3D);
    err = initialize_vbap(vbap3d, AUDIO_BLOCK_SIZE, degs3d, N_3D_LOUDSPEAKERS, SPREAD_POINTS);
    if (err != NO_ERROR) {
        printf("[ERROR] Init VBAP3D error!\n");
        return 1;
    }

    printf("[VBAP3D INFO]\n");
    for (int i = 0; i < vbap3d->vbap3d->n_hull; ++i) {
        printf("Triplets: (%d, %d, %d)\n", vbap3d->vbap3d->ltriplets[i].p0, vbap3d->vbap3d->ltriplets[i].p1, vbap3d->vbap3d->ltriplets[i].p2);
    }

    solve_gain_vector(vbap3d, &cartesian_source3d, SPREAD);

    for (int i = 0; i < vbap3d->vbap3d->n; ++i) {
        printf("GAINS: %f\n", vbap3d->vbap3d->lgains->cur_gains[i]);
    }

    PolarPoint dbap_pol_source = { 1.0, -30.0, 0.0 };
    CartesianPoint dbap_cart_source = pol_to_car(&dbap_pol_source, DEGREE);
    PANNIX *dbap = pannix_alloc(PDBAP);
    err = initialize_dbap(dbap, AUDIO_BLOCK_SIZE, dbap_degs, N_3D_LOUDSPEAKERS, 24, NULL);
    if (err != NO_ERROR) {
        printf("[ERROR] Init DBAP error!\n");
        return 1;
    }

    solve_gain_vector(dbap, &dbap_cart_source, DBAP_SPREAD);

    printf("[DBAP INFO]\n");
    for (int i = 0; i < dbap->dbap->n; ++i) {
        printf("GAINS: %f\n", dbap->dbap->lgains->cur_gains[i]);
    }

    pannix_dealloc(vbap2d);
    pannix_dealloc(vbap3d);
    pannix_dealloc(dbap);
    return 0;
}
