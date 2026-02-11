
conan install . --output-folder=build/release --build=missing -s build_type=Release

mkdir -p build

cmake -B build

cd build

make

cd ..

