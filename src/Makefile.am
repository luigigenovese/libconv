include $(top_builddir)/@LIBCONV_SOURCEDIR@/Makefile.am

lib_LIBRARIES = libconv.a

include_HEADERS = ../@LIBCONV_SOURCEDIR@/libconv*.h

#temporary compiling line for gfortran
AM_FCFLAGS = $(LIBCONV_FCFLAGS) $(FCFLAGS)

#temporary compiling line for gcc
AM_CFLAGS = $(LIBCONV_CFLAGS) $(CFLAGS)

libconv_a_SOURCES =  $(LIBCONV_SOURCES)

LIBCONV_OBJECTS_F90 = $(libconv_a_SOURCES:.f90=.o)

am_libconv_a_OBJECTS = $(LIBCONV_OBJECTS_F90:.c=.o)

.f90.o:
	$(FC) $(AM_FCFLAGS) -c $< -o $@

.c.o:
	$(CC) $(AM_CFLAGS) -c $< -o $@


#SUFFIXES = .c .f90 .o
VPATH = ../@LIBCONV_SOURCEDIR@

libconv.a: $(libconv_a_OBJECTS) $(libconv_a_DEPENDENCIES) $(EXTRA_libconv_a_DEPENDENCIES)
	rm -f libconv.a
	$(libconv_a_AR) libconv.a $(libconv_a_OBJECTS) $(libconv_a_LIBADD)
	$(RANLIB) libconv.a
