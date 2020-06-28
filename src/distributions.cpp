//
// Created by roryh on 22/06/2020.
//
#include <random>

std::mt19937 &global_engine() {
    static std::random_device rd{};
    static std::mt19937 gen(rd());
    //static std::default_random_engine e{};
    return gen;
}
