import glob
import os
import shutil

pytest_cache_directories = glob.glob('.pytest_cache')
benchmark_directories = glob.glob('.benchmarks')
prof_directories = glob.glob('prof')
build_directories = glob.glob('build')
pycache_directories = glob.glob('__pycache__') + glob.glob('src/**/__pycache__', recursive = True) + glob.glob('tests/**/__pycache__', recursive = True)
c_files = glob.glob('src/**/*.c', recursive = True)
cpp_files = glob.glob('src/**/*.cpp', recursive = True)
so_files = glob.glob('src/**/*.so', recursive = True)
pyd_files = glob.glob('src/**/*.pyd', recursive = True)

print('Cleaning:')
print('  Would remove:')
for each in pytest_cache_directories + benchmark_directories + prof_directories + build_directories + pycache_directories:
    print('    ' + each)
for each in c_files + cpp_files + so_files + pyd_files:
    print('    ' + each)
print('Proceed (Y/n)? ', end = '')
answer = input()

if answer == 'Y':
    for each in pytest_cache_directories + benchmark_directories + prof_directories + build_directories + pycache_directories:
        shutil.rmtree(each)
    for each in c_files + cpp_files + so_files + pyd_files:
        os.remove(each)
