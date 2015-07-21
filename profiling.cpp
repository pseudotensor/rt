// Declare global variables only used for profiling the code
struct rusage t_after_geodesics;            // measure time spent in geoint.cpp

time_t t_init=0.;                 // measure time spent in init() [init.cpp]
time_t t_geodesics=0;
time_t t_intensity=0.;            // measure time spent in intensity.cpp
time_t t_solvetrans=0.;           // FAST! measure time spent in solvetrans.cpp
//doub t_solvetrans=0.;           // SLOW! measure time spent in solvetrans.cpp

//time_t t_b4_evalpointzero=0.;     // auxiliary -> local
doub t_evalpointzero=0.;        // measure time spent in evalpointzero.cpp

