#include <iostream>
#include <stdlib.h>
using namespace std;
int main(/*int argc, char *argv[]*/)
{
    bool run = true;
    char c;
    unsigned int n,N;
    float sup_float, sdown_float;
    double sup_double, sdown_double;
    while(run) {
        cout << "N = ";
        cin >> N;
        sup_float = 0;
        sup_double = 0;
        for(n = 1; n <= N; n++) {
            sup_float += 1.0f/n;
            sup_double += 1.0/n;
        }
        sdown_float = 0;
        sdown_double = 0;
        for(n = N; n; n--) {
            sdown_float += 1.0f/n;
            sdown_double += 1.0/n;
        }
        cout << "sup float = " << sup_float << endl;
        cout << "sup double = " << sup_double << endl;
        cout << "sdown float = " << sdown_float << endl;
        cout << "sdown double = " << sdown_double << endl;
        cout << "Continue (y/n):";
        cin >> c;
        if(c != 'y')
            run = false;
    }
    return 0;
}

