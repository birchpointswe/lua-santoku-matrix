results.mk:
	@{ \
	  if [ -n "$$TK_MATHLIBS_CFLAGS" ]; then \
	    cflags="$$TK_MATHLIBS_CFLAGS"; \
	  elif pkg-config --cflags openblas 2>/dev/null 1>/dev/null; then \
	    cflags="$$(pkg-config --cflags openblas)"; \
	  elif command -v brew >/dev/null 2>&1 && brew --prefix openblas >/dev/null 2>&1; then \
	    cflags="-I$$(brew --prefix openblas)/include"; \
	  fi; \
	  if [ -n "$$TK_MATHLIBS_LDFLAGS" ]; then \
	    ldflags="$$TK_MATHLIBS_LDFLAGS"; \
	  elif pkg-config --libs openblas 2>/dev/null 1>/dev/null; then \
	    ldflags="$$(pkg-config --libs openblas)"; \
	  elif command -v brew >/dev/null 2>&1 && brew --prefix openblas >/dev/null 2>&1; then \
	    ldflags="-L$$(brew --prefix openblas)/lib -lopenblas"; \
	  else \
	    ldflags="-lopenblas"; \
	  fi; \
	  case "$$(uname -s)" in \
	    Linux*) ldflags="$$ldflags -Wl,-z,nodelete" ;; \
	  esac; \
	  echo "MATHLIBS_CFLAGS = $$cflags"; \
	  echo "MATHLIBS_LDFLAGS = $$ldflags"; \
	} > $@
