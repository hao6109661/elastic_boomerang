#Name of executable
noinst_PROGRAMS=hao reparametrise_beam_test beam_adapt beam_with_point_load

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
beam_with_point_load_LDADD = -L@libdir@ -lbeam -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#Sources for the executable
beam_with_point_load_SOURCES = beam_with_point_load.cc 



#Sources for the executable
hao_SOURCES = hao.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
hao_LDADD = -L@libdir@ -lbeam -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#Sources for the executable
beam_adapt_SOURCES = beam_adapt.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
beam_adapt_LDADD = -L@libdir@ -lbeam -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#Sources for the executable
reparametrise_beam_test_SOURCES = reparametrise_beam_test.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
reparametrise_beam_test_LDADD = -L@libdir@ -lbeam -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.
AM_CPPFLAGS += -I@includedir@  

