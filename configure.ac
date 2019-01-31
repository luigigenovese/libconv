AC_PREREQ(2.59)

dnl Init basic data for package
dnl define following variables :
dnl  - PACKAGE_NAME, the name ;
dnl  - PACKAGE_TARNAME, the name used for the tarball ;
dnl  - PACKAGE_VERSION, the version ;
dnl  - PACKAGE_STRING, the name + the version ;
dnl  - PACKAGE_BUGREPORT, the mail.
AC_INIT([LibConv - A library for Convolution Operations],[0.1],[luigi.genovese@cea.fr],[Conv])
AC_PROG_RANLIB()
dnl for configuring the FCFUNC macro
AC_CONFIG_HEADERS(config.h)
AM_INIT_AUTOMAKE(-Wno-portability)
prefix_pwd=$(pwd)
AC_PREFIX_DEFAULT([$prefix_pwd])

dnl search for tar program
AC_PATH_PROG([TAR], [tar])

dnl set the sources to be employed
AC_ARG_WITH(sources, AS_HELP_STRING([--with-sources], [codename of the sources to be compiled (default is 'noopt_f90') Available values are "noopt_f90", "optim_avx".]),
                           [ax_conv_sources=$withval], [ax_conv_sources="noopt_f90"])

dnl set the name of the convolution source directory
AC_SUBST([LIBCONV_SOURCEDIR],[src_${ax_conv_sources}])

dnl set for the directory to untar
UNTAR_CMD="$TAR xJvf $srcdir/${LIBCONV_SOURCEDIR}.tar.xz"
AC_CONFIG_COMMANDS_PRE([$UNTAR_CMD])

dnl If FC environment variable is not already set, it
dnl looks for a modern Fortran compiler, prefering newest
dnl fortran available (e. g. Fortran95). Then FC variable
dnl is set.
AC_PROG_FC([], [Fortran 90])
AC_FC_SRCEXT(f90)
AC_FC_WRAPPERS()

dnl Give the name of file.in to transform to file
AC_CONFIG_FILES([
Makefile
src/Makefile
])

AC_OUTPUT

dnl cat ./${LIBCONV_SOURCEDIR}/Makefile.am ./src/Makefile > toto && mv toto ./src/Makefile
