#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: glmBfp package
##        
## Time-stamp: <[configure.ac] by DSB Mon 13/12/2010 15:54 (CET)>
##
## Description:
## Autoconf script for proper handling of OpenMP dependencies. Modified from
## Simon Urbanek's proposal on the R-devel list, 14/9/2010,
## http://tolstoy.newcastle.edu.au/R/e11/devel/10/09/0683.html 
##
## History:
## 13/12/2010	file creation
#####################################################################################


# Initialization:
AC_PREREQ([2.62]) 
AC_INIT(glmBfp, 0.1, daniel.sabanesbove@ifspm.uzh.ch)

# Find R home and set CC/CFLAGS
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
RBIN="${R_HOME}/bin/R"
CC=`"${RBIN}" CMD config CC`;
CFLAGS=`"${RBIN}" CMD config CFLAGS`
LIBS="${PKG_LIBS}" 

# set RCPP_FLAGS according to Windows vs. other OS:
RCPP_FLAGS=`${R_HOME}/bin/Rscript -e "Rcpp:::LdFlags()"`

# Check for programs:
AC_PROG_CC

# Check for OpenMP:
AC_OPENMP 

# Since some systems have broken OMP libraries,
# we also check that the actual package will work:
ac_pkg_openmp=no

if test -n "${OPENMP_CFLAGS}"; then
  AC_MSG_CHECKING([whether OpenMP will work in the package])
  AC_LANG_CONFTEST([AC_LANG_PROGRAM([[#include <omp.h>]], [[ return omp_get_num_threads (); ]])]) 
  MAKEFLAGS="PKG_CFLAGS=${OPENMP_CFLAGS} PKG_LIBS=${OPENMP_CFLAGS}" "$RBIN" CMD SHLIB conftest.c 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && "$RBIN" --vanilla -q -e "dyn.load(paste('conftest',.Platform\$dynlib.ext,sep=''))" 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && ac_pkg_openmp=yes 
  dnl ac_pkg_openmp=yes
  AC_MSG_RESULT([${ac_pkg_openmp}]) 
fi

# If ${ac_pkg_openmp} = "yes" then we have OMP, otherwise it will be "no" 
if test "${ac_pkg_openmp}" = no; then
  OPENMP_CFLAGS=''
  AC_MSG_WARN([NO OpenMP support detected. You should should use gcc >= 4.2])
  # You could also put an AC_MSG_ERROR here if OpenMP is required 
fi

# Output the results:
AC_SUBST(OPENMP_CFLAGS) 
AC_SUBST(RCPP_FLAGS) 
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT 