#include <stdio.h>
#include <stdint.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>

int main(int argc, char *argv[])
{
    struct stat sb;
    
    int fd = open(argv[1], O_RDONLY);

    // get the size in bytes of the file
    fstat (fd, &sb);

    // map the file in a memory area
    void *pointer = mmap (0, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);

    char* p = (char*) pointer;
    // print 3 char of the file to demostrate it is loaded ;)
    printf("first 3 chars of the file: %c %c %c\n", p[0], p[1], p[2]);

    close(fd);
    while (true)
    {
        sleep(1000);
    }

    // detach
    //munmap(p, sb.st_size);
}
