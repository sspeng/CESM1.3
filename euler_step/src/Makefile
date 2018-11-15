OBJS = kinds.o hybrid_mod.o physical_constants.o control_mod.o dimensions_mod.o\
    element_mod.o hybvcoord_mod.o  derivative_mod.o edgetype_mod.o counter.o   \
		limiter_optim_iter_full.o 	slave_euler_step.o slave_euler_divergence.o    \
	  slave_euler_v.o  slave_euler_limit.o slave_biharmonic_wk_scalar.o          \
	  slave_biharmonic_wk_scalar_2.o slave_qdp_time_avg.o  slave_rhs_multiplier.o\
		slave_compute_overlap.o slave_update_qdp.o slave_edgeSpack_es.o            \
		slave_edgeSunpack_es.o \
		euler_step.o
all: euler_step
euler_step : $(OBJS)
	sw5f90 -hybrid kinds.o hybrid_mod.o physical_constants.o control_mod.o       \
	    dimensions_mod.o element_mod.o hybvcoord_mod.o  derivative_mod.o         \
      edgetype_mod.o counter.o limiter_optim_iter_full.o slave_euler_step.o    \
	    slave_euler_divergence.o slave_euler_v.o slave_euler_limit.o             \
	    slave_biharmonic_wk_scalar.o slave_biharmonic_wk_scalar_2.o              \
			euler_step.o slave_qdp_time_avg.o slave_rhs_multiplier.o                 \
			slave_compute_overlap.o slave_update_qdp.o slave_edgeSpack_es.o          \
			slave_edgeSunpack_es.o \
	    -o ../bin/euler_step
element_mod.o: element_mod.F90
	sw5f90 -host -c element_mod.F90
kinds.o: kinds.F90
	sw5f90 -host -c kinds.F90
dimensions_mod.o: dimensions_mod.F90
	sw5f90 -host -c dimensions_mod.F90
control_mod.o: control_mod.F90
	sw5f90 -host -c control_mod.F90
physical_constants.o: physical_constants.F90
	sw5f90 -host -c physical_constants.F90
derivative_mod.o : derivative_mod.F90
	sw5f90 -host -c derivative_mod.F90
edgetype_mod.o : edgetype_mod.F90
	sw5f90 -host -c edgetype_mod.F90
hybrid_mod.o: hybrid_mod.F90
	sw5f90 -host -c hybrid_mod.F90
hybvcoord_mod.o: hybvcoord_mod.F90
	sw5f90 -host -c hybvcoord_mod.F90
counter.o: counter.c
	sw5cc -host -c counter.c
limiter_optim_iter_full.o: limiter_optim_iter_full.F90
	sw5f90 -slave -c limiter_optim_iter_full.F90 -OPT:IEEE_arithmetic=2
slave_euler_step.o: slave_euler_step.c dma_macros.h ldm_alloc.h
	sw5cc -slave -c slave_euler_step.c -msimd
slave_euler_v.o: slave_euler_v.c dma_macros.h ldm_alloc.h
	sw5cc -slave -c slave_euler_v.c -msimd -OPT:IEEE_arithmetic=1
slave_euler_divergence.o: slave_euler_divergence.c dma_macros.h ldm_alloc.h
	sw5cc -slave -c slave_euler_divergence.c -msimd
slave_euler_limit.o: slave_euler_limit.c dma_macros.h ldm_alloc.h
	sw5cc -slave -c slave_euler_limit.c -msimd
slave_biharmonic_wk_scalar.o : slave_biharmonic_wk_scalar.c
	sw5cc -slave -c slave_biharmonic_wk_scalar.c -msimd
slave_biharmonic_wk_scalar_2.o : slave_biharmonic_wk_scalar_2.c
	sw5cc -slave -c slave_biharmonic_wk_scalar_2.c -msimd
slave_qdp_time_avg.o : slave_qdp_time_avg.c
	sw5cc -slave -c slave_qdp_time_avg.c -msimd
slave_rhs_multiplier.o : slave_rhs_multiplier.c
	sw5cc -slave -c slave_rhs_multiplier.c -msimd
slave_compute_overlap.o : slave_compute_overlap.c
	sw5cc -slave -c slave_compute_overlap.c -msimd
slave_update_qdp.o : slave_update_qdp.c
	sw5cc -slave -c slave_update_qdp.c -msimd
slave_edgeSpack_es.o : slave_edgeSpack_es.c
	sw5cc -slave -c slave_edgeSpack_es.c -msimd
slave_edgeSunpack_es.o : slave_edgeSunpack_es.c
	sw5cc -slave -c slave_edgeSunpack_es.c -msimd
euler_step.o: euler_step.F90
	sw5f90 -host -c euler_step.F90
clean:
	rm *.o  *.mod ../bin/euler_step
