project('bessel-benchmarks', 'cpp',
        default_options : ['cpp_std=c++17'])
benchmark = dependency('benchmark')
executable('bessel_benchmarks.x', 'cpp/main.cpp',
           dependencies : benchmark,
           cpp_args : ['-O3'])
