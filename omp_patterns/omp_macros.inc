#ifdef LOOP
#define TEAMS_OUTER_LOOP loop
#define PARALLEL_INNER_LOOP loop
#define COMBINED_LOOP loop
#else
#define TEAMS_OUTER_LOOP distribute
#define PARALLEL_INNER_LOOP parallel do
#define COMBINED_LOOP TEAMS_OUTER_LOOP PARALLEL_INNER_LOOP
#endif

#ifdef __NVCOMPILER_OPENMP_GPU
#define TEAM_SIZE 128
#define WAVEFRONT 32
#define TILE_I 32
#define TILE_J 4
#else
#define TEAM_SIZE 256
#define WAVEFRONT 64
#define TILE_I 64
#define TILE_J 4
#endif
