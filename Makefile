CXXFLAGS 	= -Wconversion
LIBS 		= -lfftw3 -lm
BINARY 		= testit
SOURCE		= ffts.cpp testit.cpp

-include ../shared/shared.mak