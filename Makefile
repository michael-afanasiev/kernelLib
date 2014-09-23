NETCDF_LIB = -L/usr/local/lib -lnetcdf_c++4
NETCDF_INC = -I/usr/local/include

CXX      = g++ 
CC       = gcc
CXXFLAGS = -O3 -std=c++11 $(NETCDF_INC)
CFLAGS   = -O3

LDFLAGS = $(NETCDF_LIB)

OBJS_GENERIC     = \
	./obj/main.o
								 
./obj/%.o: ./src/%.c
	@$(CC) 						$(CFLAGS) -c -o $@ $<
	@echo "CC  $<"

./obj/%.o: ./src/%.cpp
	@$(CXX) $(CXXFLAGS) -c -o $@ $<
	@echo "CXX $<"

###################
all: main

main: $(OBJS_GENERIC)
	@$(CXX) $(LDFLAGS) -o ./bin/test.exe $(OBJS_GENERIC) $(LDFLAGS)
	@echo "Linking ... $<"
	
###################	
clean:
	$(RM) ./obj/*.o
