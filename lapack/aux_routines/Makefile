F90 = gfortran
#F90 = nagfor
#F90 = ifort
FFLAGS = 
AR = ar
RM = rm -f

TARGET = libaux.a

OBJ = blas_interfaces.o \
lapack_example_aux.o \
lapack_interfaces.o \
lapack_precision.o \
nagf_blas_damax_val.o \
nagf_blas_ddscl.o \
nagf_blas_dpyth.o \
nagf_blas_zddscl.o \
nagf_blas_zmload.o \
nagf_file_print_matrix_complex_band.o \
nagf_file_print_matrix_complex_band_comp.o \
nagf_file_print_matrix_complex_gen.o \
nagf_file_print_matrix_complex_gen_comp.o \
nagf_file_print_matrix_complex_packed_comp.o \
nagf_file_print_matrix_integer_comp.o \
nagf_file_print_matrix_real_band.o \
nagf_file_print_matrix_real_gen.o \
nagf_file_print_matrix_real_gen_comp.o \
nagf_file_print_matrix_real_packed.o \
nagf_sort_cmplxvec_rank_rearrange.o \
nagf_sort_realmat_rank_rows.o \
nagf_sort_realvec_rank.o \
nagf_sort_realvec_rank_rearrange.o \
p01abf.o \
x04aafn.o \
x04abfn.o \
x04bafn.o \
x04cafn.o \
x04cbfn.o \
x04cbzn.o \
x04ccfn.o \
x04cdfn.o \
x04cefn.o \
x04cffn.o \
x04dafn.o \
x04dbfn.o \
x04ddfn.o \
x04defn.o \
x04dffn.o \
x04ebfn.o

MOD = lapack_precision.mod \
blas_interfaces.mod \
lapack_example_aux.mod \
lapack_interfaces.mod

.SUFFIXES:
.SUFFIXES: .f90 .o 
.SUFFIXES: .f90 .mod 

.f90.o:
	$(F90) $(FFLAGS) -c $< -o $@
.f90.mod:
	$(F90) $(FFLAGS) -c $< -o $*.o

$(TARGET): $(OBJ)
	$(AR) r $(TARGET) $(OBJ)


#DEPENDENCY
blas_interfaces.o: blas_interfaces.f90
blas_interfaces.o: lapack_precision.mod
blas_interfaces.mod: blas_interfaces.f90
blas_interfaces.mod: lapack_precision.mod
lapack_example_aux.o: lapack_example_aux.f90
lapack_example_aux.o: lapack_precision.mod
lapack_example_aux.mod: lapack_example_aux.f90
lapack_example_aux.mod: lapack_precision.mod
lapack_interfaces.o: lapack_interfaces.f90
lapack_interfaces.o: lapack_precision.mod
lapack_interfaces.mod: lapack_interfaces.f90
lapack_interfaces.mod: lapack_precision.mod
lapack_precision.o: lapack_precision.f90
lapack_precision.mod: lapack_precision.f90
nagf_blas_damax_val.o: nagf_blas_damax_val.f90
nagf_blas_damax_val.o: lapack_precision.mod
nagf_blas_ddscl.o: nagf_blas_ddscl.f90
nagf_blas_ddscl.o: blas_interfaces.mod
nagf_blas_ddscl.o: lapack_precision.mod
nagf_blas_dpyth.o: nagf_blas_dpyth.f90
nagf_blas_dpyth.o: lapack_precision.mod
nagf_blas_zddscl.o: nagf_blas_zddscl.f90
nagf_blas_zddscl.o: blas_interfaces.mod
nagf_blas_zddscl.o: lapack_precision.mod
nagf_blas_zmload.o: nagf_blas_zmload.f90
nagf_blas_zmload.o: lapack_precision.mod
nagf_file_print_matrix_complex_band.o: nagf_file_print_matrix_complex_band.f90
nagf_file_print_matrix_complex_band.o: lapack_precision.mod
nagf_file_print_matrix_complex_band_comp.o: nagf_file_print_matrix_complex_band_comp.f90
nagf_file_print_matrix_complex_band_comp.o: lapack_precision.mod
nagf_file_print_matrix_complex_gen.o: nagf_file_print_matrix_complex_gen.f90
nagf_file_print_matrix_complex_gen.o: lapack_precision.mod
nagf_file_print_matrix_complex_gen_comp.o: nagf_file_print_matrix_complex_gen_comp.f90
nagf_file_print_matrix_complex_gen_comp.o: lapack_precision.mod
nagf_file_print_matrix_complex_packed_comp.o: nagf_file_print_matrix_complex_packed_comp.f90
nagf_file_print_matrix_complex_packed_comp.o: lapack_precision.mod
nagf_file_print_matrix_integer_comp.o: nagf_file_print_matrix_integer_comp.f90
nagf_file_print_matrix_real_band.o: nagf_file_print_matrix_real_band.f90
nagf_file_print_matrix_real_band.o: lapack_precision.mod
nagf_file_print_matrix_real_gen.o: nagf_file_print_matrix_real_gen.f90
nagf_file_print_matrix_real_gen.o: lapack_precision.mod
nagf_file_print_matrix_real_gen_comp.o: nagf_file_print_matrix_real_gen_comp.f90
nagf_file_print_matrix_real_gen_comp.o: lapack_precision.mod
nagf_file_print_matrix_real_packed.o: nagf_file_print_matrix_real_packed.f90
nagf_file_print_matrix_real_packed.o: lapack_precision.mod
nagf_sort_cmplxvec_rank_rearrange.o: nagf_sort_cmplxvec_rank_rearrange.f90
nagf_sort_cmplxvec_rank_rearrange.o: lapack_precision.mod
nagf_sort_realmat_rank_rows.o: nagf_sort_realmat_rank_rows.f90
nagf_sort_realmat_rank_rows.o: lapack_precision.mod
nagf_sort_realvec_rank.o: nagf_sort_realvec_rank.f90
nagf_sort_realvec_rank.o: lapack_precision.mod
nagf_sort_realvec_rank_rearrange.o: nagf_sort_realvec_rank_rearrange.f90
nagf_sort_realvec_rank_rearrange.o: lapack_precision.mod
p01abf.o: p01abf.f90
x04aafn.o: x04aafn.f90
x04abfn.o: x04abfn.f90
x04bafn.o: x04bafn.f90
x04cafn.o: x04cafn.f90
x04cafn.o: lapack_example_aux.mod
x04cafn.o: lapack_precision.mod
x04cbfn.o: x04cbfn.f90
x04cbfn.o: lapack_example_aux.mod
x04cbfn.o: lapack_precision.mod
x04cbzn.o: x04cbzn.f90
x04ccfn.o: x04ccfn.f90
x04ccfn.o: lapack_example_aux.mod
x04ccfn.o: lapack_precision.mod
x04cdfn.o: x04cdfn.f90
x04cdfn.o: lapack_example_aux.mod
x04cdfn.o: lapack_precision.mod
x04cefn.o: x04cefn.f90
x04cefn.o: lapack_example_aux.mod
x04cefn.o: lapack_precision.mod
x04cffn.o: x04cffn.f90
x04cffn.o: lapack_example_aux.mod
x04cffn.o: lapack_precision.mod
x04dafn.o: x04dafn.f90
x04dafn.o: lapack_example_aux.mod
x04dafn.o: lapack_precision.mod
x04dbfn.o: x04dbfn.f90
x04dbfn.o: lapack_example_aux.mod
x04dbfn.o: lapack_precision.mod
x04ddfn.o: x04ddfn.f90
x04ddfn.o: lapack_example_aux.mod
x04ddfn.o: lapack_precision.mod
x04defn.o: x04defn.f90
x04defn.o: lapack_example_aux.mod
x04defn.o: lapack_precision.mod
x04dffn.o: x04dffn.f90
x04dffn.o: lapack_example_aux.mod
x04dffn.o: lapack_precision.mod
x04ebfn.o: x04ebfn.f90
x04ebfn.o: lapack_example_aux.mod



clean:
	$(RM) $(OBJ) $(MOD) $(TARGET)
