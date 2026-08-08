project('OpToCpp', 'cpp',
  version : '1.0',
  default_options : ['cpp_std=c++17', 'buildtype=debug']
)

# --- 1. OpenMP: The Smart Detection ---
# Meson 1.5+ handles the Apple Clang /usr/local vs /opt/homebrew search
omp_dep = dependency('openmp', required : false)

if host_machine.system() == 'darwin'
    # This dependency ONLY contains the compiler flags (-Xpreprocessor -fopenmp)
    # It will NOT be passed to the linker
    omp_compile_only = omp_dep.partial_dependency(compile_args: true, includes: true)
    
    # We manually find the library for the linker to avoid the trailing flag error
    cpp = meson.get_compiler('cpp')
    # Use find_library to get a clean link to 'omp'
    omp_lib_only = cpp.find_library('omp', dirs: '/opt/homebrew/opt/libomp/lib')
    
    omp_deps = [omp_compile_only, omp_lib_only]
else
    # Linux handles the unified dependency perfectly
    omp_deps = [omp_dep]
endif

# --- 2. Pybind11 ---
# Check system first; if missing, trigger the .wrap download
pybind11_dep = dependency('pybind11', required : false)


# --- 3. Python & Environment ---
py = import('python').find_installation()
inc = include_directories('include')

# --- 4. Building the Targets ---

# The C++ Backend (Static Library)

# Find all .cpp files in lib and its subfolders
run_command = run_command('find', 'lib', '-name', '*.cpp', check: true)
lib_sources = run_command.stdout().strip().split('\n')

cpp_backend = static_library('cpp_backend',
  lib_sources,
  include_directories : inc,
  dependencies : omp_deps
)

# The Main C++ Executable
executable('main.out',
  'src/main.cpp',
  link_with : cpp_backend,
  include_directories : inc,
  dependencies : omp_deps,
  install : true,
  install_dir: '.'
)

# The Python Module (.so)
py.extension_module('match_to_py',
  'src/pyBindings.cpp',
  include_directories : inc,
  link_with : cpp_backend,
  dependencies : [pybind11_dep] + omp_deps,
  install : true,
  install_dir: 'py'
)