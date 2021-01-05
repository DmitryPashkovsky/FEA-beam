#include <iostream>
#include "structure.h"


int main( void )
{
    init f("sys1.txt");
  //  f.info();
    bsystem S(f);
    
    S.create_right_vector();
    S.create_connection_matrix();
    S.create_global_matrix();
    S.total_condence();
    S.solve();
    S.dsp.print();
    S.all_diagram_M_N();

    S.save_result();

    system("python plot.py");
    return 0;
}

