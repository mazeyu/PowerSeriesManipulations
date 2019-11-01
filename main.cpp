#include "powerSeries.cpp"



int main() {
    a[0] = get(1, 0);
    a[1] = get(1, 0);
    b[0] = get(1, 0);
    b[1] = get(1, 0);
    multiply(a, b, 4);
    for (int i = 0; i < 16; i++) printf("%lf ", res[i].r);
    printf("\n");
    inverse(a, 4);
    for (int i = 0; i < 16; i++) printf("%lf ", res[i].r);
    printf("\n");
    derivative(a, 4);
    for (int i = 0; i < 16; i++) printf("%lf ", res[i].r);
    printf("\n");
    integ(a, 4);
    for (int i = 0; i < 16; i++) printf("%lf ", res[i].r);
    printf("\n");
    ln(a, 4);
    for (int i = 0; i < 16; i++) printf("%lf ", res[i].r);
    printf("\n");

    a[0] = get(0, 0);
    exp(a, 4);
    for (int i = 0; i < 16; i++) printf("%lf ", res[i].r);
    printf("\n");
}