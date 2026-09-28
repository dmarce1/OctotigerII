#include "octotigerII/problems/radiatingStarStructure.hpp"
#include <iostream>
#include <iomanip>
#include <numbers>
#include <cmath>
using namespace octotigerII;
int main(){
 problems::RadiatingStarStructure s({}); const double G=units::value(constants::G), c=units::value(constants::c), kap=.34;
 const double Re=units::value(s.equatorialRadius()), Rp=units::value(s.polarRadius()), M=units::value(s.mass()), om=units::value(s.angularVelocity()), dyn=std::sqrt(Re*Re*Re/(G*M));
 auto center=s.eos().atDensity(s.eos().parameters().centralDensity);
 std::cout<<std::setprecision(12)<<"M "<<M<<" Re "<<Re<<" Rp "<<Rp<<" Tc "<<units::value(center.temperature)<<" Pc "<<units::value(center.pressure)<<" Er_c "<<units::value(center.radiationEnergy)<<" Omega "<<om<<" tdyn "<<dyn<<" period "<<2*std::numbers::pi/om<<" beta_v "<<om*Re/c<<"\n";
 for(double mu:{0.,1.}){
  double R=mu==0?Re:Rp;int nr=20000;double dr=R/nr,tau=0,I=0;
  for(int i=0;i<=nr;i++){double r=i*dr;auto q=s.sample(units::Length::from_value(r*std::sqrt(1-mu*mu)),units::Length::from_value(r*mu));double weight=(i==0||i==nr?.5:1);tau+=weight*kap*units::value(q.thermodynamics.density)*dr;I+=weight*kap*units::value(q.thermodynamics.density)*r*dr;}
  std::cout<<"mu "<<mu<<" central_tau "<<tau<<" diff_simple "<<3*tau*R/c<<" diff_integral "<<3*I/c<<"\n";
  for(double f:{0.,.5,.9,.95,.98,.99,.995,.998,.999,.9995}){auto q=s.sample(units::Length::from_value(f*R*std::sqrt(1-mu*mu)),units::Length::from_value(f*R*mu));auto d=s.diffusion(q,kap);double u=1.5*units::value(q.thermodynamics.gasPressure)+units::value(q.thermodynamics.radiationEnergy); std::cout<<"f "<<f<<" rho "<<units::value(q.thermodynamics.density)<<" f0 "<<d.fluxFactor<<" tth_dyn "<<u/d.heating/dyn<<" tau_cell128 "<<kap*units::value(q.thermodynamics.density)*(2.4*Re/128)<<"\n";}
 }
 double U=0,Ug=0,Ur=0,m90=0,m95=0,m98=0;int nr=2048,nmu=64;double dr=Re/nr;
 for(int j=0;j<=nmu;j++){double mu=double(j)/nmu,wa=(j==0||j==nmu?.5:1)/nmu;
 for(int i=0;i<=nr;i++){double r=i*dr,wr=(i==0||i==nr?.5:1)*dr;auto q=s.sample(units::Length::from_value(r*std::sqrt(1-mu*mu)),units::Length::from_value(r*mu)).thermodynamics;double dv=4*std::numbers::pi*r*r*wr*wa;if(r>.9*Re)m90+=units::value(q.density)*dv;if(r>.95*Re)m95+=units::value(q.density)*dv;if(r>.98*Re)m98+=units::value(q.density)*dv;Ug+=1.5*units::value(q.gasPressure)*dv;Ur+=units::value(q.radiationEnergy)*dv;}}
 std::cout<<"mass outside .9/.95/.98 Re "<<m90/M<<" "<<m95/M<<" "<<m98/M<<"\n";U=Ug+Ur;double L=4*std::numbers::pi*G*M*c/kap;std::cout<<"Ug "<<Ug<<" Ur "<<Ur<<" U "<<U<<" Edd_L "<<L<<" U/L_dyn "<<U/L/dyn<<"\n";
}
