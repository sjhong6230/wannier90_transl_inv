include ../../make.inc
.cc.o:
	$(CXX) -g -I ../../src -c $<  -fpermissive
wannier_c.x : wannier.o ../../libwan2.a
	$(CXX) -g  wannier.o ../../libwan2.a $(LIBS) -lgfortran -o $@ -fpermissive
clean:
	@rm wannier.o wannier_c.x
