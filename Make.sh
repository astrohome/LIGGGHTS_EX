# Make.sh = update Makefile.lib or Makefile.list or style_*.h files
# Syntax: sh Make.sh style
#         sh Make.sh Makefile.lib
#         sh Make.sh Makefile.list

# function to create one style_*.h file

style () {
  list=`grep -l $1 $2*.h`
  if (test -e style_$3.tmp) then
    rm -f style_$3.tmp
  fi
  for file in $list; do
    qfile="\"$file\""
    echo "#include $qfile" >> style_$3.tmp
  done
  if (test ! -e style_$3.tmp) then
    rm -f style_$3.h
    touch style_$3.h
  elif (test ! -e style_$3.h) then
    mv style_$3.tmp style_$3.h
    rm -f Obj_*/$4.d
  elif (test "`diff --brief style_$3.h style_$3.tmp`" != "") then
    mv style_$3.tmp style_$3.h
    rm -f Obj_*/$4.d
  else
    rm -f style_$3.tmp
  fi
}

# create individual style files
# called by "make machine"

if (test $1 = "style") then

  style ANGLE_CLASS     angle_      angle      force
  style ATOM_CLASS      atom_vec_   atom       atom
  style BOND_CLASS      bond_       bond       force
  style COMMAND_CLASS   ""          command    input
  style COMPUTE_CLASS   compute_    compute    modify
  style DIHEDRAL_CLASS  dihedral_   dihedral   force
  style DUMP_CLASS      dump_       dump       output
  style FIX_CLASS       fix_        fix        modify
  style IMPROPER_CLASS  improper_   improper   force
  style INTEGRATE_CLASS ""          integrate  update
  style KSPACE_CLASS    ""          kspace     force
  style MINIMIZE_CLASS  min_        minimize   update
  style PAIR_CLASS      pair_       pair       force
  style REGION_CLASS    region_     region     domain
  style CFD_DATACOUPLING_CLASS      cfd_datacoupling_  cfd_datacoupling  fix_cfd_coupling
  style CFD_REGIONMODEL_CLASS       cfd_regionmodel_  cfd_regionmodel  fix_cfd_coupling
  style LB_CLASS        ""          lb  
  style SPH_KERNEL_CLASS  sph_kernel_  sph_kernel  pair_sph-fix_sph

# edit Makefile.lib
# called by "make makelib"
# use current list of *.cpp and *.h files in src dir w/out main.cpp

elif (test $1 = "Makefile.lib") then

  list=`ls -1 *.cpp | sed s/^main\.cpp// | tr "[:cntrl:]" " "`
  sed -i -e "s/SRC =	.*/SRC =	$list/" Makefile.lib
  list=`ls -1 *.h | tr "[:cntrl:]" " "`
  sed -i -e "s/INC =	.*/INC =	$list/" Makefile.lib

# edit Makefile.list
# called by "make makelist"
# use current list of *.cpp and *.h files in src dir

elif (test $1 = "Makefile.list") then

  list=`ls -1 *.cpp | tr "[:cntrl:]" " "`
  sed -i -e "s/SRC =	.*/SRC =	$list/" Makefile.list
  list=`ls -1 *.h | tr "[:cntrl:]" " "`
  sed -i -e "s/INC =	.*/INC =	$list/" Makefile.list

fi
