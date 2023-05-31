# catdap 1.3.7

* Removed C wrapper functions and registered entry points for the routines accessed by the `.Fortran` interface to call Fortran subroutines.

* Added legend to mosaic plot for single explanatory model.

* Fixed bugs in subroutine `catdap2mf` (catdap2f.f).
 
* Added a ‘`NEWS.md`’ file.


#  catdap 1.3.6

* Colorized plot output and added a new argument `gray.shade` to all functions.


#  catdap 1.3.5

* Fixed bugs in subroutines `catdap2mf`, `ac1p` and `mdap0` (catdap2f.f).

* Fixed bugs reproduced by checking with `--as-cran` using r-devel (r77865 or later), which causes `data.frame()` and `read.table()` to use a default `stringsAsFactors = FALSE`. (Reported by CRAN.)

* Changed to use the `warning()` function instead of the `cat()` function to print warning messages.

* In `catdap2()`, the values related to subsets of explanatory variables were collected in a list `subset`.


# catdap 1.3.4

* Fixed bugs in subroutines `catdap2mf` and `mdap0` (catdap2f.f). (reported in pre-test)

* Changed the argument name `explanatory.names` to `additional.output` and made it a list form in `catdap2()`.

* Added `missing` output value to `catdap2()` as a number indicating the type of missing value.

* Corrected the `Details' description and the comments on the examples in the for `catdap2()` help pages.

* Fixed a bug in the `ac1p' subroutine (catdap2f.f). (cause of Valgrind check error)


# catdap 1.2.4

* Fixed Fortran code according to the warning message when using gfortran with `-Wall -pedantic`.
