#ifndef TK_DUSET_EXT_H
#define TK_DUSET_EXT_H

#include <santoku/dvec.h>

static inline void tk_duset_dump (tk_duset_t *s, tk_dvec_t *v)
{
  double x;
  tk_umap_foreach_keys(s, x, ({
    tk_dvec_push(v, x);
  }))
}

#endif