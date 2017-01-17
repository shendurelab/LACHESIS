AC_INIT(lachesis, version-201701)
AC_MSG_NOTICE([A software tool to measure the thread of life.])

AC_LANG(C++)
AC_PROG_CXX
AC_CONFIG_MACRO_DIR([m4])
AC_ISC_POSIX
AC_PROG_CXX
AM_PROG_CC_STDC
AC_HEADER_STDC
AC_PROG_RANLIB
AM_PROG_LIBTOOL
AC_PROG_LIBTOOL
AM_INIT_AUTOMAKE()
AM_CONFIG_HEADER(src/config.h)
AM_MAINTAINER_MODE

AC_CHECK_HEADERS(zlib.h, [], [AC_ERROR([libz is required.])])
AC_SEARCH_LIBS(deflate, z, [], [AC_ERROR([libz is required.])])
AC_DEFINE(HAVE_LIBZ, [], [libz])

AC_CHECK_HEADERS(pthread.h, [], [AC_ERROR([pthread is required.])])
AC_SEARCH_LIBS(pthread_exit, pthread, [], [AC_ERROR([pthread is required.])])
AC_DEFINE(HAVE_PTHREAD, [], [libpthread])

AX_LIB_SAMTOOLS()

BOOST_REQUIRE([1.52.0])
BOOST_FILESYSTEM
BOOST_SYSTEM
BOOST_PROGRAM_OPTIONS
BOOST_REGEX
BOOST_THREAD

AC_CONFIG_FILES([Makefile src/Makefile src/include/Makefile \
 src/include/gtools/Makefile src/include/markov/Makefile])
AC_OUTPUT
