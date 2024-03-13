check_PROGRAMS = libconv

libconv_libraries = $(top_builddir)/src/libconv.a

#we have to reduce the dict directory such as to move this test only there
AM_LDFLAGS = -L$(top_builddir)/src
AM_FCFLAGS = -I$(top_builddir)/@LIBCONV_SOURCEDIR@

libconv_SOURCES = libconv.f90 generic.f90
libconv_LDADD = $(libconv_libraries) $(FCLIBS) 
