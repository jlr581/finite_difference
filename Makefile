# System variables
f90 = gfortran 
FFLAGS = -O3 -march=native -fno-range-check 
LINKFLAGS = -llapack -L/usr/lib64

# shouldn't need to edit below this point

rm_lst = data_types.o data_types_mod.mod \
         dq_coeff.o var_prec_extra.o var_prec_extra_mod.mod \
 	 ldl_mod.mod ldl.o var_prec.o dq_coeff_mod.mod \
	 solve.o solve_mod.mod hyper_dual.o hyper_dual_mod.mod \
         var_prec_mod.mod vector.o vector_mod.mod \
	 examples/poisson examples/poisson.o \
	 examples/poisson_regular examples/poisson_regular.o \
	 examples/cylinder examples/cylinder.o

default : examples/poisson examples/poisson_regular examples/cylinder

examples/poisson : examples/poisson.o data_types.o ldl.o solve.o hyper_dual.o var_prec_extra.o var_prec.o dq_coeff.o vector.o
	$(f90) -o $@ $^ $(FFLAGS) 

examples/poisson.o : examples/poisson.f90 data_types.o ldl.o solve.o hyper_dual.o var_prec_extra.o var_prec.o dq_coeff.o vector.o
	$(f90) -c $(FFLAGS) -o $@ $<

examples/poisson_regular : examples/poisson_regular.o data_types.o ldl.o solve.o hyper_dual.o var_prec_extra.o var_prec.o dq_coeff.o vector.o
	$(f90) -o $@ $^ $(FFLAGS) 

examples/poisson_regular.o : examples/poisson_regular.f90 data_types.o ldl.o solve.o hyper_dual.o var_prec_extra.o var_prec.o dq_coeff.o vector.o
	$(f90) -c $(FFLAGS) -o $@ $<

examples/cylinder : examples/cylinder.o data_types.o ldl.o solve.o hyper_dual.o var_prec_extra.o var_prec.o dq_coeff.o vector.o
	$(f90) -o $@ $^ $(FFLAGS) $(LINKFLAGS) 

examples/cylinder.o : examples/cylinder.f90 data_types.o ldl.o solve.o hyper_dual.o var_prec_extra.o var_prec.o dq_coeff.o vector.o
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
	
