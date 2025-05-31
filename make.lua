local env = {

  name = "santoku-matrix",
  version = "0.0.54-1",
  variable_prefix = "TK_MATRIX",
  license = "MIT",
  public = true,

  cflags = {
    "-march=native", "-std=gnu11", "-O3", "-Wall", "-Wextra", "-Wstrict-overflow",
    "-I$(shell luarocks show santoku --rock-dir)/include/",
    "-I$(shell luarocks show santoku-threads --rock-dir)/include/",
  },

  ldflags = {
    "-march=native", "-O3", "-lm", "-lpthread", "-lnuma"
  },

  dependencies = {
    "lua == 5.1",
    "santoku >= 0.0.268-1",
    "santoku-threads >= 0.0.5-1",
  },

  test = {
    dependencies = {
      "luacov == 0.15.0-1",
    }
  },

}

env.homepage = "https://github.com/treadwelllane/lua-" .. env.name
env.tarball = env.name .. "-" .. env.version .. ".tar.gz"
env.download = env.homepage .. "/releases/download/" .. env.version .. "/" .. env.tarball

return {
  type = "lib",
  env = env,
}
