FC = ifx
FCOPTS = -O3 -qopenmp
LD = ifx
LDOPTS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread 
LDLIBS = -L/storage/knapp/wigxjpf/lib -lwigxjpf
WIGX = -I/storage/knapp/wigxjpf/mod

EXENAME = ../run/HF_BCS
OBJS =HF_BCS.o type_defs.o declarations.o read_input.o basis.o HFBCS_solver.o v2b_bin.o v3b_no2b_bin.o hamiltonian.o Tcm_2b.o coupling_coefficients.o 

$(EXENAME): $(OBJS)
	$(LD) $(LDOPTS) -o $(EXENAME) $(OBJS) $(LDLIBS)

HF_BCS.o: declarations.o read_input.o basis.o HFBCS_solver.o hamiltonian.o Tcm_2b.o coupling_coefficients.o HF_BCS.f90 
	$(FC) $(FCOPTS) $(WIGX) -c HF_BCS.f90

type_defs.o: type_defs.f90 
	$(FC) $(FCOPTS) -c type_defs.f90

declarations.o: type_defs.o declarations.f90 
	$(FC) $(FCOPTS) -c declarations.f90

read_input.o: declarations.o read_input.f90 
	$(FC) $(FCOPTS) -c read_input.f90

basis.o: declarations.o basis.f90 
	$(FC) $(FCOPTS) -c basis.f90

HFBCS_solver.o: declarations.o basis.o v2b_bin.o v3b_no2b_bin.o HFBCS_solver.f90 
	$(FC) $(FCOPTS) -c HFBCS_solver.f90

hamiltonian.o: declarations.o v2b_bin.o v3b_no2b_bin.o Tcm_2b.o coupling_coefficients.o hamiltonian.f90 
	$(FC) $(FCOPTS) $(WIGX) -c hamiltonian.f90

v2b_bin.o: declarations.o v2b_bin.f90 
	$(FC) $(FCOPTS) -c v2b_bin.f90

v3b_no2b_bin.o: declarations.o v3b_no2b_bin.f90 
	$(FC) $(FCOPTS) -c v3b_no2b_bin.f90

Tcm_2b.o: declarations.o coupling_coefficients.o Tcm_2b.f90 
	$(FC) $(FCOPTS) $(WIGX) -c Tcm_2b.f90

coupling_coefficients.o: coupling_coefficients.f90
	$(FC) $(FCOPTS) $(WIGX) -c coupling_coefficients.f90 

clean:
	rm -f *.o *.mod $(EXENAME)

#clobber:
#	rm -f $(EXENAME)
