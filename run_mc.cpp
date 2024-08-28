#include "./mc_subroutine/mc_read_load_compute.hpp"
#include "./potentialFunction/potentialFunctionPrototype.hpp"

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cout << "wrong arguments" << std::endl;
        std::exit(2);
    }
//    unsigned int numThreads = std::thread::hardware_concurrency();
//    std::cout << "numThreads=" << numThreads << std::endl;
    auto mcObj=mc_computation(std::string(argv[1]));
    int N=10;
    std::shared_ptr<double> xVec(new double[N], std::default_delete<double[]>());
    for (int j=0;j<N;j++)
    {
        xVec.get()[j]=static_cast<double>(j);
    }
    mcObj.print_shared_ptr(xVec,N);

    mcObj.init_and_run();




}