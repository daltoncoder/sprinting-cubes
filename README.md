# Sprinting Cubes
An optimized marching cube algorithm to get mesh data from 3D volumetric data.

Some optimizations
- Volume Processing: Handles entire 3D volumes instead of single cubes
- Vertex Deduplication: Shared vertices between adjacent cubes are merged using a cache
- Pre-computed lookup tables (That i got from here https://github.com/therealnv6/marching-cubes/blob/main/src/tables.rs )
- Linear Interpolation: Calculates precise surface intersection points rather than using crude approximations, resulting in smooth, accurate surfaces.
- Early Exit: Skips cubes with no surface intersections
