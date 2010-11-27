# $Id: Makefile,v 1.8 2005/11/03 15:17:37 duncan Exp $
include ${FSLCONFDIR}/default.mk

PROJNAME = bet2

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_ZLIB}

LIBS=-lmeshclass -lnewimage -lprob -lmiscmaths -lfslio -lniftiio -lznz -lnewmat -lutils -lz

BETOBJS=bet2.o 
BETSURFOBJS=betsurf.o 

XFILES=bet2 betsurf

SCRIPTS = betall

all: ${XFILES}

bet2: ${BETOBJS} 
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${BETOBJS} ${LIBS} 

betsurf: ${BETSURFOBJS} 
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${BETSURFOBJS} ${LIBS} 


