#include "titanlibcustom.h"
#include <iostream>

using namespace titanlibcustom;

float* titanlibcustom::test_array(float* v, int n) {
    int count = 0;
    for(int i = 0; i < n; i++)
        count++;
    return v;
 }
void titanlibcustom::test_not_implemented_exception() {
    throw titanlibcustom::not_implemented_exception();
}
