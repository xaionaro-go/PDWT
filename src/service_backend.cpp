
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

using namespace std;

float *read_file(string path, unsigned int *amount_of_values) {
    FILE *fp;

    if((fp=fopen(path.c_str(), "rb"))==NULL) {
        fprintf(stderr, "unable to open the input file '%s'\n", path.c_str());
        exit(1);
    }

    fseek(fp, 0, SEEK_END);
    long size = ftell(fp);
    *amount_of_values = size / sizeof(float);
    fseek(fp, 0, SEEK_SET);

    float *f = (float *)malloc(size);
    if(f==NULL) {
        fprintf(stderr, "unable to allocate memory\n");
        exit(1);
    }

    if(fread(f, sizeof(char), size, fp)!=size) {
        fprintf(stderr, "unable to read the input file '%s'\n", path.c_str());
        exit(1);
    }

    fclose(fp);
    return f;
}

void write_file(string path, unsigned int amount_of_values, float *f) {
    FILE *fp;

    if((fp=fopen(path.c_str(), "wb"))==NULL) {
        fprintf(stderr, "unable to open the output file '%s'\n", path.c_str());
        exit(1);
    }

    if(fwrite(f, sizeof(float), amount_of_values, fp)!=amount_of_values) {
        fprintf(stderr, "unable to write to the output file '%s'\n", path.c_str());
        exit(1);
    }

    fclose(fp);
}

bool str2bool(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    std::istringstream is(str);
    bool b;
    is >> std::boolalpha >> b;
    return b;
}

int main(int argc, char **argv) {
    for (std::string line; std::getline(std::cin, line, '\n');) {
        string wavelet_func = "";
        unsigned int levels = 0;
        bool is_inverse = 0;
        string input_file_path = "";
        string output_file_path = "";

        string word;
        stringstream words_stream(line);
        int word_num = 0;
        while (std::getline(words_stream, word, '\000')) {
            switch (word_num) {
            case 0:
                wavelet_func = word;
                break;
            case 1:
                levels = atoi(word.c_str());
                break;
            case 2:
                is_inverse = str2bool(word);
                break;
            case 3:
                input_file_path = word;
                break;
            case 4:
                output_file_path = word;
                break;
            case 5:
                fprintf(stderr, "too many words\n");
                exit(1);
                break;
            }
            word_num++;
        }

        unsigned int amount_of_values;
        float *data = read_file(input_file_path, &amount_of_values);

        if (is_inverse) {
            Wavelets W(data, 1, amount_of_values, wavelet_func.c_str(), levels, 1, 1, 0, 0, 1);
            W.inverse();
        } else {
            Wavelets W(data, 1, amount_of_values, wavelet_func.c_str(), levels, 1, 1, 0, 0, 1);
            W.forward();
            float *coefficients = (float*)calloc(amount_of_values*levels, sizeof(float));
            float *coefficients_ptr = coefficients;
            int amount_of_coefficients = 0;
            for(int level = 0; level < levels; level++) {
                int n = W.get_coeff(coefficients_ptr, level);
                coefficients_ptr = &coefficients[n];
                amount_of_coefficients += n;
            }
            write_file(output_file_path, amount_of_coefficients, coefficients);
            printf("\n");
        }
    }

    return 0;
}
