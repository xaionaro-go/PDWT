NVCC=nvcc
#CFLAGS="-arch=sm_30"
LDFLAGS=-lcublas

PDWTCORE=src/wt.cu src/common.cu src/utils.cu src/separable.cu src/nonseparable.cu src/haar.cu src/filters.cpp
PDWTOBJ=build/wt.o build/common.o build/utils.o build/separable.o build/nonseparable.o build/haar.o build/filters.o

PROTOC = protoc
GRPC_CPP_PLUGIN = grpc_cpp_plugin
GRPC_CPP_PLUGIN_PATH ?= `which $(GRPC_CPP_PLUGIN)`
PROTOS_PATH = src/protobuf

vpath %.proto $(PROTOS_PATH)

#
# Using constant memory accross several files requires to use separate compilation (relocatable device code),
# Otherwise a new constant memory buffer is created for each file (even if the symbol is defined in a common file),
# since __constant__ variables have a file scope linkage. This was fine until the introduction of Wavelets::set_filters().
# As constant memory is managed through the use of "symbols" rather than buffers, another strategy would be to
# get the pointer address with cudaGetSymbolAddress(), which is not recommended.
#
# Separate compilation might be the way to go for better modularity, easier refactoring and compilation speed.
# However, cython does not offer flexibility to make two linkage steps (one "nvcc --dlink" to link the cuda ".o" together,
# the other to link the pyx ".o" with the linked cuda ".o").
#
# If you still want to use separate compilation :
#   - replace "-c $^" with "-dc $^" in the Makefile targets rules
#   - uncomment the definition of SEPARATE_COMPILATION in filters.h
#
demo: $(PDWTCORE) src/demo.cpp src/io.cpp
	mkdir -p build
	$(NVCC) -g $(CFLAGS) -odir build -c $^
	$(NVCC) $(CFLAGS) -o demo $(PDWTOBJ) build/demo.o build/io.o $(LDFLAGS)

%.pb.o: %.pb.cc
	mkdir -p build
	g++ -I build/ -c -o build/$@ build/$<

%.grpc.pb.cc: %.proto
	mkdir -p build
	$(PROTOC) -I $(PROTOS_PATH) --grpc_out=build/ --plugin=protoc-gen-grpc=$(GRPC_CPP_PLUGIN_PATH) $<

%.pb.cc: %.proto
	mkdir -p build
	$(PROTOC) -I $(PROTOS_PATH) --cpp_out=build/ $<

service: $(PDWTCORE) src/service.cpp service.pb.o service.grpc.pb.o
	mkdir -p build
	protoc --grpc_out=build/ --plugin=protoc-gen-grpc=`which grpc_cpp_plugin` src/protobuf/service.proto
	protoc --cpp_out=build/ src/protobuf/service.proto
	$(NVCC) -g $(CFLAGS) -I build -odir build -c $^
	$(NVCC) $(CFLAGS) $(shell pkg-config --libs protobuf) $(shell pkg-config --libs grpc++) build/service.pb.o build/service.grpc.pb.o -o service $(PDWTOBJ) build/service.o $(LDFLAGS)

libpdwt.so: $(PDWTCORE)
	mkdir -p build
	$(NVCC) $(CFLAGS) --ptxas-options=-v --compiler-options '-fPIC' -odir build -c $^
	$(NVCC) $(CFLAGS)  -o $@ --shared $(PDWTOBJ) $(LDFLAGS)


# Double precision library
libpdwtd.so: $(PDWTCORE)
	mkdir -p build
	$(NVCC) --ptxas-options=-v --compiler-options '-fPIC' -DDOUBLEPRECISION -odir build -c $^
	$(NVCC) $(CFLAGS)  -o $@ --shared $(PDWTOBJ) $(LDFLAGS)



clean:
	rm -rf build demo libpdwt*.so service
