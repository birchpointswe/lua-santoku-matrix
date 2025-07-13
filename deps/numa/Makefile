HAVE_NUMA := $(shell echo '#include <numa.h>' | $(CC) -x c - -o /dev/null 2>/dev/null && echo 1 || echo 0)

results.mk:
	touch results.mk
	[ "$(HAVE_NUMA)" -eq 1 ] && \
		echo 'LIB_LDFLAGS += -lnuma' >> results.mk || true

