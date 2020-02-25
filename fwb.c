#include <Python.h>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>

#include <numpy/arrayobject.h>

static PyObject* query_mmap(PyObject* NPY_UNUSED(self), PyObject* args) {
   // parse arguments
   char* fwb_path;
   uint8_t* width;
   PyArray_Dims* offsets_obj;
   if(!PyArg_ParseTuple(
     args, "sbO&",
     &fwb_path, &width, PyArray_IntpConverter, &offsets_obj
   )) return NULL;

   // memory map file
   // TODO: consider using huge pages?
   struct stat sb;
   int fwb_fd = open(fwb_path, O_RDONLY);
   fstat(fwb_fd, &sb);
   uint8_t* map = (uint8_t*) mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fwb_fd, 0);

   // copy into buffer
   npy_intp n_offsets = offsets_obj->len;
   npy_intp* offsets = offsets_obj->ptr;
   uint8_t* buf = calloc(n_offsets, *width);
   if(buf == NULL) {
      fprintf(stderr, "Couldn't allocate output buffer!\n");
      exit(1);
   }
   for(int64_t i = 0; i < n_offsets; i++) {
      // FIXME: handle user-specified null value? or can we just do this in Python?
      if(offsets[i] < 0) {
	 buf += *width;
	 continue;
      }
      for(int64_t j = 0; j < *width; j++) *(buf++) = map[offsets[i] + j];
   }

   // create numpy array
   npy_intp dims[1] = {1};
   int output_type;
   switch(*width) {
      case 2: { output_type = NPY_UINT16; break; }
      case 4: { output_type = NPY_UINT32; break; }
      case 8: { output_type = NPY_UINT64; break; }
   }

   PyObject* buf_np = PyArray_SimpleNewFromData(n_offsets, dims, output_type, buf);

   // clean up
   PyDimMem_FREE(offsets);
   munmap(map, sb.st_size);
   close(fwb_fd);

   return buf_np;
}
