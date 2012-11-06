#include "config.h"
#include <stdlib.h>
#include "gsl_vector.h"

#include "view.h"

#define BASE_DOUBLE
#include "templates_on.h"
#include "view_source.c"
#include "templates_off.h"
#undef  BASE_DOUBLE
