
include ${FSLCONFDIR}/default.mk

PROJNAME = avwutils

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_ZLIB}

LIBS = -lnewimage -lmiscmaths -lprob -lfslio -lnewmat -lutils -lniftiio -lznz -lm -lz

IOLIBS = -lfslio -lniftiio -lznz -lm -lz
SSLIBS = -lss_32R $(IOLIBS)

XFILES = avwcpgeom avwmerge avwnvols avwsplit \
	 avw2ascii_8UI \
         avw2ascii_16SI avw2ascii_32SI \
         avw2ascii_32R avw2ascii_64R avwhd avwcreatehd \
	 avwmaths_8UI avwmaths_16SI \
         avwmaths_32SI avwmaths_32R avwmaths_64R \
	 avwcc_8UI avwcc_16SI \
         avwcc_32SI avwcc_32R avwcc_64R \
	 avwinterleave_8UI \
         avwinterleave_16SI avwinterleave_32SI \
         avwinterleave_32R avwinterleave_64R avwroi_8UI \
         avwroi_16SI avwroi_32SI avwroi_32R avwroi_64R \
	 avwstats avwcomplex \
         avwstats++ avwfft avwmeants avwcheck avwswapdim avwconv \
         avwslice avworient avwpspec avwroi++ avwnvols++ avw2ascii++ \
         avwsplit++ avwmerge++ avwcc++ avwinterleave++ \
         avwhd++ avwcpgeom++ avwcreatehd++ avwmaths++

SCRIPTS = avwval avwchpixdim avwanimate sliceanimate avwsize avwinfo avwedithd avwchfiletype
FSCRIPTS = fmriborient

RUNAVWS = avw2ascii avwmaths avwcc avwroi avwfill avwinterleave

all: ${XFILES}


avwcheck: avwcheck.cc
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o avwcheck avwcheck.cc ${LIBS}

avwconv: avwconv.cc
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o avwconv avwconv.cc ${LIBS}

avwstats++: avwstats++.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ avwstats++.o ${LIBS} 

avwfft: avwfft.cc
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o avwfft avwfft.cc ${LIBS}

avwmeants: avwmeants.cc
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o avwmeants avwmeants.cc ${LIBS}

avwfixfloat: avwfixfloat.cc
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o avwfixfloat avwfixfloat.cc ${LIBS}

avwcomplex: avwcomplex.cc
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o avwcomplex avwcomplex.cc ${LIBS}

avwnvols: avwnvols.c
	${CC} ${CFLAGS} ${LDFLAGS} -o avwnvols avwnvols.c ${LIBS}

avwsplit: avwsplit.c
	${CC} ${CFLAGS} ${LDFLAGS} -o avwsplit avwsplit.c ${LIBS}

avwmerge: avwmerge.c
	${CC} ${CFLAGS} ${LDFLAGS} -o avwmerge avwmerge.c ${LIBS}

avwhd: avwhd.c
	$(CC) $(CFLAGS) -o avwhd avwhd.c $(LDFLAGS) ${LIBS}

avwcreatehd: avwcreatehd.c
	$(CC) $(CFLAGS) -o avwcreatehd avwcreatehd.c $(LDFLAGS) ${SSLIBS}

avwcpgeom: avwcpgeom.c
	$(CC) $(CFLAGS) -o avwcpgeom avwcpgeom.c $(LDFLAGS) ${SSLIBS}

avwstats: avwstats.c
	$(CC) $(CFLAGS) -DFDT="float" -o avwstats avwstats.c $(LDFLAGS) ${SSLIBS}

avwroi_8UI: avwroi.c
	for FDT in 8UI 16SI 32SI 32R 64R ; do $(DATATYPES) ; \
	$(CC) $(CFLAGS) -DFDT="$$STR" -o avwroi_$$FDT avwroi.c $(LDFLAGS) -lss_$$FDT $(IOLIBS); \
	done

avwmaths_8UI: avwmaths.c
	for FDT in 8UI 16SI 32SI 32R 64R ; do $(DATATYPES) ; \
	$(CC) $(CFLAGS) -DFDT="$$STR" -o avwmaths_$$FDT avwmaths.c $(LDFLAGS) -lss_$$FDT $(IOLIBS); \
	done

avwcc_8UI: avwcc.c
	for FDT in 8UI 16SI 32SI 32R 64R ; do $(DATATYPES) ; \
	$(CC) $(CFLAGS) -DFDT="$$STR" -o avwcc_$$FDT avwcc.c $(LDFLAGS) -lss_$$FDT $(IOLIBS); \
	done

avwinterleave_8UI: avwinterleave.c
	for FDT in 8UI 16SI 32SI 32R 64R ; do $(DATATYPES) ; \
	$(CC) $(CFLAGS) -DFDT="$$STR" -o avwinterleave_$$FDT avwinterleave.c $(LDFLAGS) -lss_$$FDT $(IOLIBS); \
	done

avw2ascii_8UI: avw2ascii.c
	for FDT in 8UI 16SI 32SI 32R 64R ; do $(DATATYPES) ; \
	$(CC) $(CFLAGS) -DFDT="$$STR" -o avw2ascii_$$FDT avw2ascii.c $(LDFLAGS) -lss_$$FDT $(IOLIBS); \
	done

avwswapdim: avwswapdim.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ avwswapdim.o ${LIBS} 

avworient: avworient.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ avworient.o ${LIBS} 

avwslice: avwslice.cc
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o avwslice avwslice.cc ${LIBS}

avwpspec: avwpspec.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ avwpspec.o ${LIBS} 


#start of matts c++ rewrites
avwroi++: avwroi++.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ avwroi++.o ${LIBS}

avwnvols++: avwnvols++.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ avwnvols++.o ${LIBS}

avw2ascii++: avw2ascii++.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ avw2ascii++.o ${LIBS}

avwsplit++: avwsplit++.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ avwsplit++.o ${LIBS}

avwmerge++: avwmerge++.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ avwmerge++.o ${LIBS}

avwcc++: avwcc++.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ avwcc++.o ${LIBS}

avwinterleave++: avwinterleave++.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ avwinterleave++.o ${LIBS}

avwhd++: avwhd++.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ avwhd++.o ${LIBS}

avwcpgeom++: avwcpgeom++.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ avwcpgeom++.o ${LIBS}

avwcreatehd++: avwcreatehd++.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ avwcreatehd++.o ${LIBS}

avwmaths++: avwmaths++.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ avwmaths++.o ${LIBS}

