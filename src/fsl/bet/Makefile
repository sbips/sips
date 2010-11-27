include $(FSLCONFDIR)/default.mk

USRINCFLAGS = -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_ZLIB}

PROJNAME = bet

LIBS = -lss_$$FDT -lfslio -lniftiio -lznz -lm -lz

XFILES = bet_8UI bet_16SI bet_32SI bet_32R

RUNAVWS = bet

RUNTCLS = Bet

SCRIPTS = betfunc betpremask betsmallz

all:	bet_8UI

bet_8UI: bet.c
	for FDT in 8UI 16SI 32SI 32R ; do $(DATATYPES) ; \
        $(CC) $(CFLAGS) -DFDT="$$STR" -o bet_$$FDT bet.c $(LDFLAGS) $(LIBS) ; \
	done
