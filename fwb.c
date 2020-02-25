#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>

int main() {
   // read in offsets
   int64_t* offsets = malloc(43928650*8);
   FILE* o = fopen("/mnt/j/proj/mut/LPWGS/20190612_readlevel/offsets.bin", "r");
   fread(offsets, 8, 43928650, o);
   fclose(o);

   // memory map file
   // TODO: consider using huge pages?
   struct stat sb;
   int fwb_fd = open("/mnt/j/db/hg19/map/align75.fwb", O_RDONLY);
   fstat(fwb_fd, &sb);
   uint8_t* map = (uint8_t*) mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fwb_fd, 0);

   // read file
   uint8_t width = 2;
   
   uint8_t* buf = calloc(43928650, width);
   if(buf == NULL) {
      fprintf(stderr, "Couldn't allocate output buffer!\n");
      exit(1);
   }
   for(int64_t i = 0; i < 43928650; i++) {
      // FIXME: handle user-specified null value? or can we just do this in Python?
      if(offsets[i] < 0) {
	 buf += width;
	 continue;
      }
      for(int64_t j = 0; j < width; j++) *(buf++) = map[offsets[i] + j];
   }
}
