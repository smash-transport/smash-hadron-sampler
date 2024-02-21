#include "virtest/vir/test.h"

#include "../src/include/gen.h"

using namespace gen;

TEST(test) {
    int variable = test(7.9);
    VERIFY(variable == 7);
}