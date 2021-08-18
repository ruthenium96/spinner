//
// Created by radium on 17.08.2021.
//

#include "Task.h"
#include "tz_sorter.h"
#include "c2_symmetrizer.h"


int main() {

    std::vector<int> mults = {2, 2, 2, 2};

    Task T(mults);

    T.print();

    T = tz_sorter(T);

    T.print();

    T = c2_symmetrizer(T, 2);

    T.print();

}