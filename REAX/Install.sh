# Install/unInstall package files in LAMMPS
# edit Makefile.package to include/exclude REAX library

if (test $1 = 1) then

  sed -i -e 's/[^ \t]*reax //' ../Makefile.package
  sed -i -e 's/[^ \t]*reax_[^ \t]*) //' ../Makefile.package
  sed -i -e 's|^PKG_INC =[ \t]*|&-I../../lib/reax |' ../Makefile.package
  sed -i -e 's|^PKG_PATH =[ \t]*|&-L../../lib/reax |' ../Makefile.package
  sed -i -e 's|^PKG_LIB =[ \t]*|&-lreax |' ../Makefile.package
  sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(reax_SYSPATH) |' ../Makefile.package
  sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(reax_SYSLIB) |' ../Makefile.package

  cp pair_reax.cpp ..
  cp pair_reax.h ..
  cp pair_reax_fortran.h ..

  cp fix_reax_bonds.h ..
  cp fix_reax_bonds.cpp ..

elif (test $1 = 0) then

  sed -i -e 's/[^ \t]*reax //' ../Makefile.package
  sed -i -e 's/[^ \t]*reax_[^ \t]*) //' ../Makefile.package

  rm ../pair_reax.cpp
  rm ../pair_reax.h
  rm ../pair_reax_fortran.h

  rm ../fix_reax_bonds.h
  rm ../fix_reax_bonds.cpp

fi
