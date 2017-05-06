#include "udf.h"
DEFINE_PROFILE(pressure_profile,t,i)
{
real x[ND_ND];
/* this will hold the position vector */
real y;
face_t f;
begin_f_loop(f,t)
{
F_CENTROID(x,f,t);
y = x[1];
F_PROFILE(f,t,i) = 1.1e5 - y*y/(.0745*.0745)*0.1e5;
}
end_f_loop(f,t)
}
