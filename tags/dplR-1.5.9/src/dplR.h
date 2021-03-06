#ifndef DPLR_H
#define DPLR_H

#include <R.h>  /* to include Rconfig.h */
#include <Rversion.h>
#include <Rinternals.h>
size_t dplRlength(SEXP x);
          
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("dplR", String)
#else
#define _(String) (String)
#define dngettext(pkg, String, StringP, N) (N > 1 ? StringP: String)
#endif

#if defined(R_VERSION) && R_VERSION >= R_Version(3, 0, 0)
#define DPLR_RGEQ3
#endif

#endif
