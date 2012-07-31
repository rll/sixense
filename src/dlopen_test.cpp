#include <iostream>
#include <dlfcn.h>
#include <sixense.h>
using namespace std;

int main() {

  typedef void (*VoidFunc)();
  sixenseAllControllerData acd;
  void* handle = dlopen("/home/joschu/Dropbox/myros/sixense/lib/libsixense_x64.so",RTLD_LAZY);
  VoidFunc _sixenseGetAllNewestData = (VoidFunc)dlsym(handle, "sixenseGetAllNewestData");
  //_sixenseGetAllNewestData(&acd);

}
