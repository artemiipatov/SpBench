////////////////////////////////////////////////////////////////////////////////////
// MIT License                                                                    //
//                                                                                //
// Copyright (c) 2021 Egor Orachyov                                               //
//                                                                                //
// Permission is hereby granted, free of charge, to any person obtaining a copy   //
// of this software and associated documentation files (the "Software"), to deal  //
// in the Software without restriction, including without limitation the rights   //
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell      //
// copies of the Software, and to permit persons to whom the Software is          //
// furnished to do so, subject to the following conditions:                       //
//                                                                                //
// The above copyright notice and this permission notice shall be included in all //
// copies or substantial portions of the Software.                                //
//                                                                                //
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     //
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       //
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE    //
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         //
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  //
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE  //
// SOFTWARE.                                                                      //
////////////////////////////////////////////////////////////////////////////////////

#include "benchmark_base.hpp"
#include "matrix_loader.hpp"
#include "args_processor.hpp"
#include "profile_mem.hpp"

extern "C"
{
#include "suitesparse/LAGraph.h"
};

#define BENCH_DEBUG
#define GrB_CHECK(func) do { auto s = func; assert(s == GrB_SUCCESS); } while(0);

namespace benchmark {
    class Multiply : public BenchmarkBase {
    public:

        Multiply(int argc, const char** argv) {
            argsProcessor.parse(argc, argv);
            assert(argsProcessor.isParsed());

            benchmarkName = "SuiteSparse-Multiply";
            experimentsCount = argsProcessor.getExperimentsCount();
        }

        ~Multiply() {
            output_mem_profile(benchmarkName + "-Mem.txt", argsProcessor.getInputString());
        }

    protected:

        void setupBenchmark() override {
            GrB_CHECK(GrB_init(GrB_BLOCKING));
        }

        void tearDownBenchmark() override {
            GrB_CHECK(GrB_finalize());
        }

        void setupExperiment(size_t experimentIdx, size_t& iterationsCount, std::string& name) override {
//            std::cout << "setup" << std::endl;

            src_count = 50;

            src = (GrB_Index*)std::malloc(sizeof(GrB_Index) * src_count);

            for (auto i = 0; i < src_count; i++) {
                src[i] = (GrB_Index) i;
            }

            auto& entry = argsProcessor.getEntries()[experimentIdx];

            iterationsCount = entry.iterations;
            name = entry.name;

            const auto& file = entry.name;
            const auto& type = entry.isUndirected;

            MatrixLoader loader(file, type);
            loader.loadData();
            input = std::move(loader.getMatrix());

#ifdef BENCH_DEBUG
            log << ">   Load matrix: \"" << file << "\" isUndirected: " << type << std::endl
                << "                 size: " << input.nrows << " x " << input.ncols << " nvals: " << input.nvals << std::endl;
#endif // BENCH_DEBUG

            size_t n = input.nrows;
            assert(input.nrows == input.ncols);

            GrB_CHECK(GrB_Matrix_new(&A, GrB_INT32, n, n));

            std::vector<GrB_Index> I(input.nvals);
            std::vector<GrB_Index> J(input.nvals);

            int* X = (int*)std::malloc(sizeof(int) * input.nvals);

            for (auto i = 0; i < input.nvals; i++) {
                I[i] = input.rows[i];
                J[i] = input.cols[i];
                X[i] = 1;
            }

            GrB_CHECK(GrB_Matrix_build_INT32(A, I.data(), J.data(), X, input.nvals, GrB_FIRST_INT32));

            std::free(X);

            std::string msgString = "graph has not been created";
            char msg[27];
            strcpy(msg, msgString.c_str());

            GrB_CHECK(LAGraph_New(&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg))
        }

        void tearDownExperiment(size_t experimentIdx) override {
            input = Matrix{};

            GrB_CHECK(GrB_Matrix_free(&A));
            A = nullptr;
        }

        void setupIteration(size_t experimentIdx, size_t iterationIdx) override {
            GrB_CHECK(GrB_Matrix_new(&L, GrB_INT32, src_count, input.ncols));
            GrB_CHECK(GrB_Matrix_new(&P, GrB_INT32, src_count, input.ncols));
        }

        void execIteration(size_t experimentIdx, size_t iterationIdx) override {
            std::string msgString = "msbfs failure";
            char msg[14];
            strcpy(msg, msgString.c_str());

            GrB_CHECK(LAGr_BreadthFirstSearch(&L, &P, G, src, src_count, msg));
        }

        void tearDownIteration(size_t experimentIdx, size_t iterationIdx) override {
            GrB_Index nrows;
            GrB_Index ncols;
            GrB_Index nvals;

            GrB_CHECK(GrB_Matrix_nrows(&nrows, L));
            GrB_CHECK(GrB_Matrix_ncols(&ncols, L));
            GrB_CHECK(GrB_Matrix_nvals(&nvals, L));

#ifdef BENCH_DEBUG
            log << "   Result matrix: size " << nrows << " x " << ncols
                << " nvals " << nvals << std::endl;
#endif

            GrB_CHECK(GrB_Matrix_free(&L));
            GrB_CHECK(GrB_Matrix_free(&P));
            L = nullptr;
            P = nullptr;
        }

    protected:

        GrB_Matrix A;
        LAGraph_Graph G;
        GrB_Matrix L;
        GrB_Matrix P;
        GrB_Index *src;
        int src_count;

        ArgsProcessor argsProcessor;
        Matrix input;
    };

}

int main(int argc, const char** argv) {
    benchmark::Multiply multiply(argc, argv);
    multiply.runBenchmark();
    return 0;
}