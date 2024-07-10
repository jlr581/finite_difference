# System variables
f90 = gfortran 
FFLAGS = -O3 -march=native -fno-range-check 
# add the "-fopenmp" flag to FFLAGS for OPENMP parallelization with gfortran 
LINKFLAGS = -llapack -L/usr/lib64

# shouldn't need to edit below this point

rm_lst = data_types.o data_types_mod.mod \
         dq_coeff.o var_prec_extra.o var_prec_extra_mod.mod \
 	 ldl_mod.mod ldl.o var_prec.o dq_coeff_mod.mod \
	 solve.o solve_mod.mod hyper_dual.o hyper_dual_mod.mod \
         var_prec_mod.mod vector.o vector_mod.mod \
	 examples/poisson_64bit_native/poisson examples/poisson_64bit_native/*.o  examples/poisson_64bit_native/*.mod\
	 examples/poisson_128bit_native/poisson examples/poisson_128bit_native/*.o examples/poisson_128bit_native/*.mod\
	 examples/poisson_64bit_vp/poisson examples/poisson_64bit_vp/poisson.o \
	 examples/cylinder_64bit_vp/cylinder_64bit examples/cylinder_64bit_vp/cylinder_64bit.o \
	 examples/cylinder_32bit_vp/cylinder_32bit examples/cylinder_32bit_vp/cylinder_32bit.o \
	 examples/triangle_64bit_vp/triangle examples/triangle_64bit_vp/triangle.o \
	 examples/triangle_finite_difference/finite_difference examples/triangle_finite_difference/finite_difference.o 

default : examples/poisson_64bit_native/poisson examples/poisson_128bit_native/poisson examples/poisson_64bit_vp/poisson examples/cylinder_32bit_vp/cylinder_32bit examples/cylinder_64bit_vp/cylinder_64bit examples/triangle_64bit_vp/triangle examples/triangle_finite_difference/finite_difference

examples/poisson_64bit_native/poisson : examples/poisson_64bit_native/poisson.o examples/poisson_64bit_native/data_types.o examples/poisson_64bit_native/ldl.o examples/poisson_64bit_native/solve.o examples/poisson_64bit_native/hyper_dual.o examples/poisson_64bit_native/dq_coeff.o examples/poisson_64bit_native/vector.o
	$(f90) -o $@ $^ $(FFLAGS) 

examples/poisson_64bit_native/poisson.o : examples/poisson_64bit_native/poisson.f90 examples/poisson_64bit_native/data_types.o examples/poisson_64bit_native/ldl.o examples/poisson_64bit_native/solve.o examples/poisson_64bit_native/hyper_dual.o examples/poisson_64bit_native/dq_coeff.o examples/poisson_64bit_native/vector.o
	$(f90) -c $(FFLAGS) -o $@ $<

examples/poisson_64bit_native/data_types.o : examples/poisson_64bit_native/data_types.f90 
	$(f90) -c $(FFLAGS) -o $@ $<

examples/poisson_64bit_native/ldl.o : examples/poisson_64bit_native/ldl.f90 examples/poisson_64bit_native/data_types.o examples/poisson_64bit_native/hyper_dual.o 
	$(f90) -c $(FFLAGS) -o $@ $<

examples/poisson_64bit_native/solve.o : examples/poisson_64bit_native/solve.f90 examples/poisson_64bit_native/ldl.o examples/poisson_64bit_native/data_types.o examples/poisson_64bit_native/hyper_dual.o 
	$(f90) -c $(FFLAGS) -o $@ $<

examples/poisson_64bit_native/dq_coeff.o : examples/poisson_64bit_native/dq_coeff.f90 examples/poisson_64bit_native/ldl.o examples/poisson_64bit_native/data_types.o examples/poisson_64bit_native/hyper_dual.o 
	$(f90) -c $(FFLAGS) -o $@ $<

examples/poisson_64bit_native/hyper_dual.o : examples/poisson_64bit_native/hyper_dual.f90 examples/poisson_64bit_native/data_types.o 
	$(f90) -c $(FFLAGS) -o $@ $<

examples/poisson_64bit_native/vector.o : examples/poisson_64bit_native/vector.f90 examples/poisson_64bit_native/data_types.o 
	$(f90) -c $(FFLAGS) -o $@ $<

examples/poisson_128bit_native/poisson : examples/poisson_128bit_native/poisson.o examples/poisson_128bit_native/data_types.o examples/poisson_128bit_native/ldl.o examples/poisson_128bit_native/solve.o examples/poisson_128bit_native/hyper_dual.o examples/poisson_128bit_native/dq_coeff.o examples/poisson_128bit_native/vector.o
	$(f90) -o $@ $^ $(FFLAGS) 

examples/poisson_128bit_native/poisson.o : examples/poisson_128bit_native/poisson.f90 examples/poisson_128bit_native/data_types.o examples/poisson_128bit_native/ldl.o examples/poisson_128bit_native/solve.o examples/poisson_128bit_native/hyper_dual.o examples/poisson_128bit_native/dq_coeff.o examples/poisson_128bit_native/vector.o
	$(f90) -c $(FFLAGS) -o $@ $<

examples/poisson_128bit_native/data_types.o : examples/poisson_128bit_native/data_types.f90 
	$(f90) -c $(FFLAGS) -o $@ $<

examples/poisson_128bit_native/ldl.o : examples/poisson_128bit_native/ldl.f90 examples/poisson_128bit_native/data_types.o examples/poisson_128bit_native/hyper_dual.o 
	$(f90) -c $(FFLAGS) -o $@ $<

examples/poisson_128bit_native/solve.o : examples/poisson_128bit_native/solve.f90 examples/poisson_128bit_native/ldl.o examples/poisson_128bit_native/data_types.o examples/poisson_128bit_native/hyper_dual.o 
	$(f90) -c $(FFLAGS) -o $@ $<

examples/poisson_128bit_native/dq_coeff.o : examples/poisson_128bit_native/dq_coeff.f90 examples/poisson_128bit_native/ldl.o examples/poisson_128bit_native/data_types.o examples/poisson_128bit_native/hyper_dual.o 
	$(f90) -c $(FFLAGS) -o $@ $<

examples/poisson_128bit_native/hyper_dual.o : examples/poisson_128bit_native/hyper_dual.f90 examples/poisson_128bit_native/data_types.o 
	$(f90) -c $(FFLAGS) -o $@ $<

examples/poisson_128bit_native/vector.o : examples/poisson_128bit_native/vector.f90 examples/poisson_128bit_native/data_types.o 
	$(f90) -c $(FFLAGS) -o $@ $<

examples/poisson_64bit_vp/poisson : examples/poisson_64bit_vp/poisson.o data_types.o ldl.o solve.o hyper_dual.o var_prec_extra.o var_prec.o dq_coeff.o vector.o
	$(f90) -o $@ $^ $(FFLAGS)

examples/poisson_64bit_vp/poisson.o : examples/poisson_64bit_vp/poisson.f90 data_types.o ldl.o solve.o hyper_dual.o var_prec_extra.o var_prec.o dq_coeff.o vector.o
	$(f90) -c $(FFLAGS) -o $@ $<

examples/cylinder_32bit_vp/cylinder_32bit : examples/cylinder_32bit_vp/cylinder_32bit.o data_types.o ldl.o solve.o hyper_dual.o var_prec_extra.o var_prec.o dq_coeff.o vector.o
	$(f90) -o $@ $^ $(FFLAGS) $(LINKFLAGS) 

examples/cylinder_32bit_vp/cylinder_32bit.o : examples/cylinder_32bit_vp/cylinder_32bit.f90 data_types.o ldl.o solve.o hyper_dual.o var_prec_extra.o var_prec.o dq_coeff.o vector.o
	$(f90) -c $(FFLAGS) $(LINKFLAGS) -o $@  $<

examples/cylinder_64bit_vp/cylinder_64bit : examples/cylinder_64bit_vp/cylinder_64bit.o data_types.o ldl.o solve.o hyper_dual.o var_prec_extra.o var_prec.o dq_coeff.o vector.o
	$(f90) -o $@ $^ $(FFLAGS) $(LINKFLAGS) 

examples/cylinder_64bit_vp/cylinder_64bit.o : examples/cylinder_64bit_vp/cylinder_64bit.f90 data_types.o ldl.o solve.o hyper_dual.o var_prec_extra.o var_prec.o dq_coeff.o vector.o
	$(f90) -c $(FFLAGS) $(LINKFLAGS) -o $@  $<

examples/triangle_64bit_vp/triangle : examples/triangle_64bit_vp/triangle.o data_types.o ldl.o solve.o hyper_dual.o var_prec_extra.o var_prec.o dq_coeff.o vector.o
	$(f90) -o $@ $^ $(FFLAGS) $(LINKFLAGS) 

examples/triangle_64bit_vp/triangle.o : examples/triangle_64bit_vp/triangle.f90 data_types.o ldl.o solve.o hyper_dual.o var_prec_extra.o var_prec.o dq_coeff.o vector.o
	$(f90) -c $(FFLAGS) $(LINKFLAGS) -o $@  $<

examples/triangle_finite_difference/finite_difference : examples/triangle_finite_difference/finite_difference.o 
	$(f90) -o $@ $^ $(FFLAGS) $(LINKFLAGS) 

examples/triangle_finite_difference/finite_difference.o : examples/triangle_finite_difference/finite_difference.f90 
	$(f90) -c $(FFLAGS) $(LINKFLAGS) -o $@  $<

data_types.o : data_types.f90 
	$(f90) -c $(FFLAGS) $<

ldl.o : ldl.f90 data_types.o hyper_dual.o var_prec_extra.o var_prec.o
	$(f90) -c $(FFLAGS) $<

solve.o : solve.f90 ldl.o data_types.o hyper_dual.o var_prec_extra.o var_prec.o
	$(f90) -c $(FFLAGS) $<

dq_coeff.o : dq_coeff.f90 ldl.o data_types.o hyper_dual.o var_prec_extra.o var_prec.o
	$(f90) -c $(FFLAGS) $<

hyper_dual.o : hyper_dual.f90 data_types.o var_prec_extra.o var_prec.o
	$(f90) -c $(FFLAGS) $<

var_prec.o : var_prec.f90 var_prec_extra.o data_types.o 
	$(f90) -c $(FFLAGS) $<

var_prec_extra.o : var_prec_extra.f90 data_types.o 
	$(f90) -c $(FFLAGS) $<

vector.o : vector.f90 data_types.o 
	$(f90) -c $(FFLAGS) $<

.PHONY: clean
clean: 
	-@rm -f $(rm_lst)
	
