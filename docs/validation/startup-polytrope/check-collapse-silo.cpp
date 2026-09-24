#include <silo.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
int main(int argc, char** argv) {
 if(argc != 2) return 1;
 auto* f=DBOpen(argv[1],DB_HDF5,DB_READ);
 if(!f) return 2;
 auto* blocks=DBGetMultivar(f,"density");
 if(!blocks) return 3;
 double peak=0, error=0;
 int cells=0;
 for(int b=0;b<blocks->nvars;++b) {
  auto* q=DBGetQuadvar(f,blocks->varnames[b]);
  auto* mesh=DBGetQuadmesh(f,("block"+std::to_string(b)+"/mesh").c_str());
  if(!q || !mesh || q->datatype!=DB_DOUBLE || mesh->datatype!=DB_DOUBLE) return 4;
  auto* rho=static_cast<double*>(q->vals[0]);
  auto* x=static_cast<double*>(mesh->coords[0]);
  auto* y=static_cast<double*>(mesh->coords[1]);
  auto* z=static_cast<double*>(mesh->coords[2]);
  int i=0;
  for(int iz=0;iz<mesh->dims[2]-1;++iz)
   for(int iy=0;iy<mesh->dims[1]-1;++iy)
    for(int ix=0;ix<mesh->dims[0]-1;++ix,++i) {
     double cx=(x[ix]+x[ix+1])/2,cy=(y[iy]+y[iy+1])/2,cz=(z[iz]+z[iz+1])/2;
     double expected=1e4*(0.01+std::exp(-(cx*cx+cy*cy+cz*cz)/(2*3e8*3e8)));
     peak=std::max(peak,rho[i]);
     error=std::max(error,std::abs(rho[i]-expected));
     ++cells;
    }
  DBFreeQuadvar(q); DBFreeQuadmesh(mesh);
 }
 std::cout<<std::setprecision(15)<<"blocks="<<blocks->nvars<<" cells="<<cells<<" peakDensity="<<peak<<" maxAbsDensityError="<<error<<'\n';
 DBFreeMultivar(blocks); DBClose(f);
 return error<1e-9 && peak>9900 ? 0 : 5;
}
