local env = {

  name = "santoku-matrix",
  version = "0.0.39-1",
  variable_prefix = "TK_MATRIX",
  license = "MIT",
  public = true,

  cflags = { "-I$(shell luarocks show santoku --rock-dir)/include" },
  ldflags = { },

  dependencies = {
    "lua == 5.1",
    "santoku >= 0.0.263-1",
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
