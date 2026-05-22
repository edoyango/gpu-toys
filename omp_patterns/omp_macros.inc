#ifdef ACC
#define GPU !$acc
#define INIT_TEAMS parallel
#define END_TEAMS end parallel
#define NTEAMS num_gangs
#define THREAD_LIMIT(...) vector_length(__VA_ARGS__)
#define TEAMS_OUTER_LOOP loop
#define PARALLEL_INNER_LOOP loop
#define COMBINED_LOOP loop
#define ENTER_DATA enter data
#define EXIT_DATA exit data
#define MAP_TO(...) copyin(__VA_ARGS__)
#define MAP_FROM(...) copyout(__VA_ARGS__)
#define MAP_TOFROM(...) copy(__VA_ARGS__)
#define MAP_ALLOC(...) create(__VA_ARGS__)
#define MAP_DELETE(...) delete(__VA_ARGS__)
#define SYNC wait
#else
#define GPU !$omp
#define INIT_TEAMS target teams
#define END_TEAMS end target teams
#define ENTER_DATA target enter data
#define EXIT_DATA target exit data
#define NTEAMS num_teams
#define THREAD_LIMIT(...) thread_limit(__VA_ARGS__)
#define MAP_TO(...) map(to: __VA_ARGS__)
#define MAP_FROM(...) map(from: __VA_ARGS__)
#define MAP_TOFROM(...) map(tofrom: __VA_ARGS__)
#define MAP_ALLOC(...) map(alloc: __VA_ARGS__)
#define MAP_DELETE(...) map(delete: __VA_ARGS__)
#define SYNC taskwait
#ifdef LOOP
#define TEAMS_OUTER_LOOP loop
#define PARALLEL_INNER_LOOP loop
#define COMBINED_LOOP loop
#else
#define TEAMS_OUTER_LOOP distribute
#define PARALLEL_INNER_LOOP parallel do
#define COMBINED_LOOP TEAMS_OUTER_LOOP PARALLEL_INNER_LOOP
#endif
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
