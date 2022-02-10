# Readme

`xtess.f90` and `xqtess.f90` are the revised Fortran codes in double and quadruple precision based on the original Fortran codes `xtess.txt` and `xqtess.txt` (Fukushima 2018). 

They can accurately compute the gravitational field (i.e. gravitational potential, gravitational acceleration vector, and gravity gradient tensor) of a tesseroid both in double and quadruple precision no matter the computation point is located outside, near the surface of, on the surface of, or inside the tesseroid.

In the abstract of Fukushima (2018):
>  The method numerically integrates a surface integral representation of the gravitational potential of the tesseroid by conditionally splitting its line integration intervals and by using the double exponential quadrature rule. Then, it evaluates the gravitational acceleration vector and the gravity gradient tensor by numerically differentiating the numerically integrated potential. 
> 
> The achievable precision is 14–15 digits for the potential, 9–11 digits for the acceleration vector, and 6–8 digits for the gradient tensor in the double precision environment. The correct digits are roughly doubled if employing the quadruple precision computation.

The original sites of `xtess.txt` and `xqtess.txt` are:
- [xtess.txt (Fortran program to compute the gravitational field of a tesseroid accurately in double precision)](https://www.researchgate.net/publication/319442456_xtesstxt_Fortran_program_to_compute_the_gravitational_field_of_a_tesseroid_accurately_in_double_precision)
- [xqtess.txt (Fortran program to compute the gravitational field of a tesseroid accurately in quadruple precision)](https://www.researchgate.net/publication/319442459_xqtesstxt_Fortran_program_to_compute_the_gravitational_field_of_a_tesseroid_accurately_in_quadruple_precision)

# Fortran compilation suggestions

If the Fortran code `xtess.f90` (or `xqtess.f90`) runs on Linux system (e.g. Ubuntu) with the `ifort`, please use this command:

```
ifort xtess.f90 -o xtess && 
./xtess
```
or

```
ifort xqtess.f90 -o xqtess && 
./xqtess
```

If the Fortran code `xtess.f90` (or `xqtess.f90`) runs on the macOS system with the `ifort`, please use this command:

```
ifort -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib xtess.f90 -o xtess && 
./xtess
```

or

```
ifort -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib xqtess.f90 -o xqtess && 
./xqtess
```

Regarding the installation of the `ifort`, please refer to the Intel® oneAPI Base & HPC Toolkit:
- [Intel® oneAPI Base Toolkit: Essential oneAPI Tools & Libraries](https://www.intel.cn/content/www/cn/zh/developer/tools/oneapi/base-toolkit.html)
- [Intel® oneAPI HPC Toolkit: Cluster & HPC Development Tools](https://www.intel.cn/content/www/cn/zh/developer/tools/oneapi/hpc-toolkit.html)

# Reference

- Fukushima T (2018) Accurate computation of gravitational field of a tesseroid. Journal of Geodesy 92(12):1371–1386, https://doi.org/10.1007/s00190-018-1126-2