from conan import ConanFile
from conan.tools.cmake import cmake_layout

class OrbitalSimConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "CMakeDeps", "CMakeToolchain"

    def requirements(self):
        self.requires("sdl/3.2.20")
        self.requires("cjson/1.7.19")
        self.requires("opengl/system")

        if self.settings.os != "Emscripten":
            self.requires("glew/2.2.0")

    def configure(self):
        # Build SDL3 as shared library on Windows to avoid naming issues
        # (static builds produce SDL3-static.lib which CMakeDeps can't find)
        if self.settings.os == "Windows":
            self.options["sdl"].shared = True

    def layout(self):
        cmake_layout(self)
