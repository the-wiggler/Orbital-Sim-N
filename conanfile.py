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
    
    def layout(self):
        cmake_layout(self)
