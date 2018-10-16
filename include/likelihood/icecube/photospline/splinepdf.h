#ifndef LOGSPLINEPDF
#define LOGSPLINEPDF

#include <gsl/gsl_rng.h>

#include "helper/photospline/splinetable.h"

#ifdef __cplusplus
extern "C" {
#endif

void logsplinepdf_n_sample(value_t *result, index_t results, index_t burnin,
    value_t *coords, index_t dim, struct splinetable *table, index_t derivatives,
    value_t (* proposal)(void*), value_t (* proposal_pdf)(value_t, value_t, void*),
    void *proposal_info, const gsl_rng *rng);

void splinepdf_n_sample(value_t *result, index_t results, index_t burnin,
    value_t *coords, index_t dim, struct splinetable *table, index_t derivatives,
    value_t (* proposal)(void*), value_t (* proposal_pdf)(value_t, value_t, void*),
    void *proposal_info, const gsl_rng *rng);

#ifdef __cplusplus
}
#endif

#endif
