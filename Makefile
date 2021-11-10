# Don't get CXXFLAGS from the environment variable, as the build doesn't work with -O2.
CXXFLAGS?=		-g -gdwarf-3 -fpermissive -Wall -O0
CXX?=		g++
INCLUDES=	-I.
OBJS=		kthread.o bbf.o htab.o bseq.o misc.o sys.o \
		    kalloc.o paf.o hic_mapping.o hic_mapping_haplo.o hic_completeness.o hic_qv.o count.o hic_switch_error.o
PROG=		pstools
LIBS=		-lm -lz -lpthread ./libminimap2.a ./libz.a

ifneq ($(asan),)
	CXXFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.SUFFIXES:.c .o
.PHONY:all clean test depend

.c.o:
		$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

pstools: $(OBJS) main.o
	$(CXX) $(CXXFLAGS) main.cpp misc.o bbf.o bseq.o htab.o hic_mapping.o hic_mapping_haplo.o hic_completeness.o hic_switch_error.o hic_qv.o count.o kalloc.o paf.o seqio.o seqhash.o  $(UTILS_OBJS) $(LDFLAGS) $(LIBS) -o $@
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
count.o: kthread.h yak-priv.h yak.h kseq.h
# modmap.o: modmap.c seqio.o seqhash.o $(UTILS_OBJS)

main.o: main.cpp ketopt.h bubble_chain.h paf_intersect.h seqio.o seqhash.o resolve_repeat_haplotype.h $(UTILS_OBJS)
