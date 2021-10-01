#include <numa.h>
#include <stdio.h>
int main() {
    if (numa_available() <0) {
        printf("Your system does not support NUMA API\n");
    } else {
        printf("Your system supports NUMA API\n");
    }
    return 0;
}
