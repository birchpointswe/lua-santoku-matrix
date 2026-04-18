local env = {
  name = "santoku-matrix",
  version = "0.0.307-1",
  variable_prefix = "TK_MATRIX",
  license = "MIT",
  public = true,
  cflags = {
    "-std=gnu11", "-D_GNU_SOURCE", "-Wall", "-Wextra",
    "-Wstrict-overflow", "-Wsign-conversion", "-Wsign-compare",
    "-I$(shell luarocks show santoku --rock-dir)/include/",
  },
  ldflags = {
    "-lm",
  },
  native = {
    cflags = {
      "-fopenmp",
      "$(MATHLIBS_CFLAGS)",
    },
    ldflags = {
      "-fopenmp",
      "$(MATHLIBS_LDFLAGS)",
    },
  },
  build = {
    wasm = {
      ldflags = {
        "-sWASM_BIGINT",
      },
    },
  },
  test = {
    wasm = {
      ldflags = {
        "-sWASM_BIGINT",
      },
    },
  },
  dependencies = {
    "lua == 5.1",
    "santoku >= 0.0.328-1",
  },
}

env.homepage = "https://github.com/birchpointswe/lua-" .. env.name
env.tarball = env.name .. "-" .. env.version .. ".tar.gz"
env.download = env.homepage .. "/releases/download/" .. env.version .. "/" .. env.tarball

return { env = env }
