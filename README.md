# **CPANNIX (2D-3D PANNER)**

cpannix `cpnx` is a C library for spatial audio rendering, implementing:
- VBAP 2D (Vector Base Amplitude Panning)
- VBAP 3D
- DBAP (Distance Based Amplitude Panning

Pure C implementation, lightweight and dependency-free, designed for realtime processing with odular architecture and easy integration into other projects. 
For convex hull computation (used internally by VBAP), cpannix includes the header-only library [convhull_3d](https://github.com/leomccormack/convhull_3d) by Leo McCormack (already embedded in the source tree and requires no separate installation). 

`cpnx` is designed for realtime audio engines, research applications, DSP development, and language bindings (e.g. Python). 

While VBAP follows its classical formulation, the DBAP implementation in cpannix introduces an optional spatial focus control designed to shape the perceived source width in a smooth and musically natural way. 
In standard DBAP, energy is distributed across loudspeakers purely as a function of distance. This produces a stable and predictable spatial image, but sometimes a sharper focus or a tighter localization is desirable. To achieve this, without introducing discontinuities or abrupt changes, the library optionally applies a smooth exponential attenuation term to the gain computation:

$$g_i = k \frac{w_i}{b_i d_i^{a}} e^{-\beta d_i}$$  

with:  

$$\beta = \frac{\text{spread}}{\text{mean distance}} + \eta$$

What the spread parameter does:
- spread = 0 means: The algorithm behaves exactly like standard DBAP
- increasing spread: The spatial image becomes progressively tighter and more focused. Energy concentrates more strongly around the loudspeakers closest to the virtual source
- smooth behavior: The exponential term ensures that this focusing effect is continuous and free from switching artifacts or spatial discontinuities

By scaling the focusing factor with the mean loudspeaker distance, the perceived behavior remains consistent across small and large arrays. In other words, the same spread value produces comparable perceptual results regardless of the system’s physical size.  

The library provides:
- Gain vector computation
- Audio frame interpolation
- Cartesian and polar coordinate support
- Planar / interleaved channel modes
- Source spread control


## Building

To build the shared library:

```bash
make shared
```

This will generate in the `target` folder:
- macOS → `libcpnx.dylib`
- Linux → `libcpnx.so`
- Windows (MinGW) → `cpnx.dll`

## Compile and run test

```bash
make test
make run_test
```

## Usage

```c
PANNIX *vbap3d = pannix_alloc(PVBAP3D);

// initialize vbap (or dbap)
err = initialize_vbap(vbap3d, AUDIO_BLOCK_SIZE, degs3d, N_3D_LOUDSPEAKERS, SPREAD_POINTS);
if (err != NO_ERROR) {
    printf("[ERROR] Init VBAP3D error!\n");
    return 1;
}

printf("[VBAP3D INFO]\n");
for (int i = 0; i < vbap3d->vbap3d->n_hull; ++i) {
    printf("Triplets: (%d, %d, %d)\n", vbap3d->vbap3d->ltriplets[i].p0, vbap3d->vbap3d->ltriplets[i].p1, vbap3d->vbap3d->ltriplets[i].p2);
}

// get gain vector
solve_gain_vector(vbap3d, &cartesian_source3d, SPREAD);

for (int i = 0; i < vbap3d->vbap3d->n; ++i) {
    printf("GAINS: %f\n", vbap3d->vbap3d->lgains->cur_gains[i]);
}
```
