
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wt.h"

#include <grpc/grpc.h>
#include <grpcpp/security/server_credentials.h>
#include <grpcpp/server.h>
#include <grpcpp/server_builder.h>
#include <grpcpp/server_context.h>
#include "service.grpc.pb.h"

using namespace std;
using grpc::Server;
using grpc::ServerBuilder;
using grpc::ServerContext;
using grpc::ServerReader;
using grpc::ServerReaderWriter;
using grpc::ServerWriter;
using grpc::Status;

const char *waveletName(Wavelet wavelet) {
    switch(wavelet) {
    case Wavelet::DB2:
        return "db2";
    case Wavelet::DB3:
        return "db3";
    case Wavelet::DB4:
        return "db4";
    case Wavelet::DB5:
        return "db5";
    case Wavelet::DB6:
        return "db6";
    case Wavelet::DB7:
        return "db7";
    case Wavelet::DB8:
        return "db8";
    case Wavelet::DB9:
        return "db9";
    case Wavelet::DB10:
        return "db10";
    case Wavelet::DB11:
        return "db11";
    case Wavelet::DB12:
        return "db12";
    case Wavelet::DB13:
        return "db13";
    case Wavelet::DB14:
        return "db14";
    case Wavelet::DB15:
        return "db15";
    case Wavelet::DB16:
        return "db16";
    case Wavelet::DB17:
        return "db17";
    case Wavelet::DB18:
        return "db18";
    case Wavelet::DB19:
        return "db19";
    case Wavelet::DB20:
        return "db20";
    case Wavelet::SYM2:
        return "sym2";
    case Wavelet::SYM3:
        return "sym3";
    case Wavelet::SYM4:
        return "sym4";
    case Wavelet::SYM5:
        return "sym5";
    case Wavelet::SYM6:
        return "sym6";
    case Wavelet::SYM7:
        return "sym7";
    case Wavelet::SYM8:
        return "sym8";
    case Wavelet::SYM9:
        return "sym9";
    case Wavelet::SYM10:
        return "sym10";
    case Wavelet::SYM11:
        return "sym11";
    case Wavelet::SYM12:
        return "sym12";
    case Wavelet::SYM13:
        return "sym13";
    case Wavelet::SYM14:
        return "sym14";
    case Wavelet::SYM15:
        return "sym15";
    case Wavelet::SYM16:
        return "sym16";
    case Wavelet::SYM17:
        return "sym17";
    case Wavelet::SYM18:
        return "sym18";
    case Wavelet::SYM19:
        return "sym19";
    case Wavelet::SYM20:
        return "sym20";
    case Wavelet::COIF1:
        return "coif1";
    case Wavelet::COIF2:
        return "coif2";
    case Wavelet::COIF3:
        return "coif3";
    case Wavelet::COIF4:
        return "coif4";
    case Wavelet::COIF5:
        return "coif5";
    case Wavelet::BIOR1_3:
        return "bior1.3";
    case Wavelet::BIOR1_5:
        return "bior1.5";
    case Wavelet::BIOR2_2:
        return "bior2.2";
    case Wavelet::BIOR2_4:
        return "bior2.4";
    case Wavelet::BIOR2_6:
        return "bior2.6";
    case Wavelet::BIOR2_8:
        return "bior2.8";
    case Wavelet::BIOR3_1:
        return "bior3.1";
    case Wavelet::BIOR3_3:
        return "bior3.3";
    case Wavelet::BIOR3_5:
        return "bior3.5";
    case Wavelet::BIOR3_7:
        return "bior3.7";
    case Wavelet::BIOR3_9:
        return "bior3.9";
    case Wavelet::BIOR4_4:
        return "bior4.4";
    case Wavelet::BIOR5_5:
        return "bior5.5";
    case Wavelet::BIOR6_8:
        return "bior6.8";
    case Wavelet::RBIO1_3:
        return "rbio1.3";
    case Wavelet::RBIO1_5:
        return "rbio1.5";
    case Wavelet::RBIO2_2:
        return "rbio2.2";
    case Wavelet::RBIO2_4:
        return "rbio2.4";
    case Wavelet::RBIO2_6:
        return "rbio2.6";
    case Wavelet::RBIO2_8:
        return "rbio2.8";
    case Wavelet::RBIO3_1:
        return "rbio3.1";
    case Wavelet::RBIO3_3:
        return "rbio3.3";
    case Wavelet::RBIO3_5:
        return "rbio3.5";
    case Wavelet::RBIO3_7:
        return "rbio3.7";
    case Wavelet::RBIO3_9:
        return "rbio3.9";
    case Wavelet::RBIO4_4:
        return "rbio4.4";
    case Wavelet::RBIO5_5:
        return "rbio5.5";
    case Wavelet::RBIO6_8:
        return "rbio6.8";
    case Wavelet::HAAR:
        return "haar";
    }
    return "";
}

#if DTYPE != float
#error "only float supported as DTYPE"
// see memcpy-s below
#endif
class ServiceImpl final : public WTService::Service {
public:
    Status Forward(
        ServerContext *context,
        const WTForwardRequest *request,
        WTForwardResponse *response
    ) override {
        float *data = (float*)calloc(request->values_size(), sizeof(float));
        memcpy(data, request->values().data(), sizeof(float)*request->values_size());
        Wavelets W(data, request->values_size(), 1, waveletName(request->wavelet()), request->levels(), 1, 1, 0, 0, 1);
        W.forward();
        Coefficients *coefficients = new Coefficients();
        float *values = (float *)calloc(request->values_size(), sizeof(float));
        for (int levelNum = 0; levelNum<=request->levels(); levelNum++) {
            int n = W.get_coeff(values, levelNum);
            Level *level = coefficients->add_levels();
            google::protobuf::RepeatedField< float >*level_value;
            level_value->Resize(n, 0);
            for (int idx = 0; idx<n; idx++) {
                level_value->Set(idx, values[idx]);
            }

        }
        free(values);
        response->set_allocated_coefficients(coefficients);
        return Status::OK;
    }
};


int main(int argc, char **argv) {
    std::string server_address("0.0.0.0:11210");
    ServiceImpl service;

    ServerBuilder builder;
    builder.AddListeningPort(server_address, grpc::InsecureServerCredentials());
    builder.RegisterService(&service);
    std::unique_ptr<Server> server(builder.BuildAndStart());
    std::cout << "Server listening on " << server_address << std::endl;
    server->Wait();
    return 0;
}
