# Don't get CXXFLAGS from the environment variable, as the build doesn't work with -O2.
CXXFLAGS?=		-g -gdwarf-3 -fpermissive -Wall -O0
CXX?=		g++
CC = gcc
VPATH = .:FASTK:ClassPro
CFLAGS = -O2 -g -Wall -Wextra -Wno-unused-function
INCLUDES:=	-I. -I./FASTK -I./ClassPro
OBJS=		kthread.o bbf.o htab.o bseq.o misc.o sys.o \
		    kalloc.o paf.o hic_mapping.o hic_mapping_haplo.o hic_completeness.o hic_qv.o count.o hic_switch_error.o centro_asm.o
PROG=		pstools
LIBS=		-lm -lz -lpthread ./libminimap2.a ./libz.a

CLASSPRO_DIR = ClassPro
CLASSPRO_DEPS = ClassPro.c ClassPro.h kseq.h bessel.h libfastk.h gene_core.h \
				DB.h QV.h gsl/include/gsl/gsl_multifit.h
shared_OBJS = libfastk.o
CLASSPRO_OBJS = DB.o QV.o
CLASSPRO_LIBS = -L$(CLASSPRO_DIR)/gsl/lib -lgsl -lgslcblas -lm -lz -lpthread

FASTK_DIR = FASTK
FASTK_DEPS = FastK.c FastK.h io.c split.c FK_count.c table.c merge.c io.c gene_core.c gene_core.h MSDsort.c LSDsort.c \
			libfastk.c libfastk.h Profex.c
FASTK_OBJS = io.o split.o FK_count.o table.o merge.o MSDsort.o LSDsort.o Profex.o $(shared_OBJS)
FASTK_LIBS = $(FASTK_DIR)/LIBDEFLATE/libdeflate.a $(FASTK_DIR)/HTSLIB/libhts.a -lpthread -lz -lm -lbz2 -llzma -lcurl

INCLUDES+= -I$(CLASSPRO_DIR)/gsl/include -I$(FASTK_DIR)/HTSLIB

# $(FASTK_DIR)/deflate.lib: LIBDEFLATE
# 	cd $(FASTK_DIR)/LIBDEFLATE; make; cd ../..

# $(FASTK_DIR)/libhts.a: HTSLIB
# 	cd $(FASTK_DIR)/HTSLIB; make libhts.a; cd ../..

# $(FASTK_DIR)/HTSLIB/htslib_static.mk:
# 	cd $(FASTK_DIR)/HTSLIB; make htslib_static.mk; cd ../..

ifneq ($(asan),)
	CXXFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=addresss
endif

.SUFFIXES:.c .cpp .o
.PHONY:all clean test depend

all:$(PROG)

%.o: %.c
		$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

%.o: %.cpp
		$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

pstools: $(OBJS) main.o ClassPro.o $(CLASSPRO_OBJS) FASTK.o $(FASTK_OBJS)
	$(CXX) $(CXXFLAGS) main.cpp misc.o bbf.o bseq.o htab.o hic_mapping.o hic_mapping_haplo.o \
	hic_completeness.o hic_switch_error.o hic_qv.o count.o kalloc.o paf.o seqio.o seqhash.o centro_asm.o \
	ClassPro.o $(CLASSPRO_OBJS) FASTK.o $(FASTK_OBJS) $(UTILS_OBJS) $(LDFLAGS) $(LIBS) $(CLASSPRO_LIBS) \
	$(FASTK_LIBS) -o $@
clean:
		rm -fr gmon.out *.o ext/*.o a.out $(PROG) *~  *.dSYM session*
test:
		./test/run.sh

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CXXFLAGS) $(DFLAGS) -- *.c)



UTILS_OBJS=hash.o dict.o array.o utils.o
UTILS_HEADERS=utils.h array.h dict.h hash.h
$(UTILS_OBJS): utils.h $(UTILS_HEADERS)

seqhash.o: seqhash.h

seqio.o: seqio.h

kalloc.o: kalloc.h
paf.o: paf.h
# sketch.o: sketch.cpp
# htab.o: htab.cpp
bbf.o: yak.h
bseq.o: bseq.h kseq.h
htab.o: htab.cpp kthread.h yak-priv.h yak.h khashl.h
kthread.o: kthread.h
misc.o: yak-priv.h yak.h
sys.o: sys.cpp yak-priv.h yak.h
hic_mapping.o: htab.o hic_mapping.h kthread.h ketopt.h bseq.h yak-priv.h yak.h
hic_mapping_haplo.o: htab.o hic_mapping.h kthread.h ketopt.h bseq.h yak-priv.h yak.h
hic_completeness.o: htab.o hic_mapping.h kthread.h ketopt.h bseq.h yak-priv.h yak.h
hic_switch_error.o: htab.o hic_mapping.h kthread.h ketopt.h bseq.h yak-priv.h yak.h
hic_qv.o: kthread.h yak-priv.h yak.h bseq.h
centro_asm.o: centro_asm.cpp centro_asm.hpp FASTK.o ClassPro.o
count.o: kthread.h yak-priv.h yak.h kseq.h
# modmap.o: modmap.c seqio.o seqhash.o $(UTILS_OBJS)

ClassPro.o : $(wildcard $(CLASSPRO_DIR)/$(CLASSPRO_DEPS)) $(CLASSPRO_OBJS) $(shared_OBJS)
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

$(CLASSPRO_OBJS): %.o : $(CLASSPRO_DIR)/%.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

FASTK.o: $(wildcard $(FASTK_DIR)/$(FASTK_DEPS)) $(FASTK_OBJS)
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

$(FASTK_OBJS): %.o : $(FASTK_DIR)/%.c 
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

$(FASTK_DIR)/libfastk.c : $(FASTK_DIR)/gene_core.c
$(FASTK_DIR)/libfastk.h : $(FASTK_DIR)/gene_core.h

main.o: main.cpp ketopt.h bubble_chain.h paf_intersect.h seqio.o seqhash.o resolve_repeat_haplotype.h $(UTILS_OBJS)
