CFLAGS=  -m64 -g -w
CXX=g++


ILOG= /opt/ibm/ILOG/CPLEX_Studio1271/
CPPFLAGS= -DIL_STD -I$(ILOG)/cplex/include -I$(ILOG)/concert/include
CPLEXLIB=-L$(ILOG)/cplex/lib/x86-64_linux/static_pic -lilocplex -lcplex -L$(ILOG)/concert/lib/x86-64_linux/static_pic -lconcert -lm -lpthread

comp-mono:  
	$(CXX) $(CFLAGS) $(CPPFLAGS) -o mono  01FA2-mono.cpp   $(CPLEXLIB) 

clean:
	rm -f  *.out *.aux *.log *.nav *.snm *.out *.toc 
