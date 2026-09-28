#include "octotigerII/problems/radiatingStarAtmosphere.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <sstream>

namespace octotigerII::problems {
namespace {
constexpr Real pi = std::numbers::pi_v<Real>;
Real val(auto q) { return units::value(q); }
Real eddington(Real f) { return Real(1)/3+2*f*f/(2+std::sqrt(4-3*f*f)); }
Real eddingtonDerivative(Real f) { return 2*f/std::sqrt(4-3*f*f); }
Real hermite(Real a, Real da, Real b, Real db, Real h, Real t) {
	return (1-t)*(1-t)*((1+2*t)*a+h*t*da)+t*t*((3-2*t)*b+h*(t-1)*db);
}
bool positive(Real x) { return std::isfinite(x) && x>0; }
}

RadiatingStarAtmosphere::Thermodynamics RadiatingStarAtmosphere::thermodynamics(Real density) const {
	auto const s=eos_.atDensity(units::Density::from_value(std::max(density,Real(1e-30))*val(parameters_.eos.centralDensity)));
	Real const pc=val(eos_.centralPressure()), gamma=1+1/parameters_.eos.index;
	Real const derivative=gamma*std::pow(std::max(density,Real(1e-30)),gamma-1);
	Real const second=derivative*(gamma-1)/std::max(density,Real(1e-30));
	Real const db=s.radiationPressureSecondDerivative*val(eos_.centralIntegralH())/parameters_.eos.index
		*std::pow(std::max(density,Real(1e-30)),1/parameters_.eos.index-1);
	return {val(s.temperature),val(s.radiationEnergy)/pc,3*derivative*s.radiationPressureDerivative,
		val(s.gasPressure)/pc,derivative*(1-s.radiationPressureDerivative),
		3*(second*s.radiationPressureDerivative+derivative*db),second*(1-s.radiationPressureDerivative)-derivative*db};
}

Real RadiatingStarAtmosphere::opacity(Real density) const {
	if (density<=parameters_.cutoffDensityFraction) return 0;
	Real const taper=1-parameters_.cutoffDensityFraction/density;
	return parameters_.opacityLengthScale*taper*taper;
}

std::array<Real,3> RadiatingStarAtmosphere::coreRhs(Real r,std::array<Real,3> const& y) const {
	Real const rho=y[0],m=y[1],f=y[2];
	if (!(rho>0) || !(f>0 && f<1.05)) throw std::runtime_error("Radiating-star core left its physical integration domain");
	auto const t=thermodynamics(rho);
	Real const chi=opacity(rho)*rho, flux=f*t.energy, d=eddington(f), dp=eddingtonDerivative(f);
	Real const rhoPrime=(-(parameters_.eos.index+1)*rho*m/(r*r)+chi*flux)/t.pressureDerivative;
	Real const fPrime=(-chi*flux-d*t.energyDerivative*rhoPrime-6*f*f/(2+std::sqrt(4-3*f*f))*t.energy/r)/(t.energy*dp);
	return {rhoPrime,r*r*rho,fPrime};
}

std::array<Real,3> RadiatingStarAtmosphere::envelopeRhs(Real r,std::array<Real,3> const& y) const {
	Real const h=std::max(y[0],Real(0));
	Real const rho=std::pow(h*(envelopeGamma_-1)/(envelopeGamma_*envelopeK_),1/(envelopeGamma_-1));
	Real const f=vacuumFluxFactor(r),d=eddington(f),dp=eddingtonDerivative(f);
	return {-(parameters_.eos.index+1)*y[1]/(r*r),r*r*rho,f/r*(d-1)/(d-f*dp)};
}

RadiatingStarAtmosphere::RadiatingStarAtmosphere(Parameters p) : parameters_(p),eos_(p.eos) {
	if (!positive(p.opacityLengthScale) || !(p.cutoffDensityFraction>0 && p.cutoffDensityFraction<1)
		|| !positive(p.tolerance) || !positive(p.maximumRadialStep)
		|| !(p.centerSeriesRadius>0 && p.centerSeriesRadius<.05) || p.maximumSteps<10)
		throw std::invalid_argument("Invalid radiating-star atmosphere integration parameters");
	auto const tc=thermodynamics(1);
	Real const pgE=tc.pressureDerivative/tc.energyDerivative;
	Real const qa=Real(2.5)/tc.energy,qb=opacity(1)*(1+1/(3*pgE));
	Real const qd=(p.eos.index+1)/(9*pgE);
	centerFluxDerivative_=2*qd/(qb+std::sqrt(qb*qb+4*qa*qd));
	centerDensitySecondDerivative_=(-(p.eos.index+1)/3+opacity(1)*centerFluxDerivative_)/tc.pressureDerivative;
	Real const rstart=p.centerSeriesRadius;
	core_.push_back({0,0,{1,0,0},{0,0,centerFluxDerivative_/tc.energy}});
	using Vector=boost::numeric::ublas::vector<Real>;
	using Matrix=boost::numeric::ublas::matrix<Real>;
	auto integrate=[&](std::vector<Node>& nodes, std::array<Real,3> initial, Real initialRadius, Real target, bool core) {
		auto rhs=[&](Vector const& x,Vector& dx,Real r) {
			auto const result=core?coreRhs(r,{x[0],x[1],x[2]}):envelopeRhs(r,{x[0],x[1],x[2]});
			for(int k=0;k<3;++k)dx[k]=result[k];
		};
		auto jacobian=[&](Vector const& x,Matrix& jac,Real r,Vector& dr) {
			for(int i=0;i<3;++i)for(int j=0;j<3;++j)jac(i,j)=0;
			Real const n1=p.eos.index+1;
			if(core){
				Real const rho=x[0],m=x[1],f=x[2],chi=opacity(rho)*rho;
				Real const chiPrime=rho>p.cutoffDensityFraction?p.opacityLengthScale*(1-std::pow(p.cutoffDensityFraction/rho,2)):0;
				auto const t=thermodynamics(rho);auto const derivative=coreRhs(r,{rho,m,f});
				Real const d=eddington(f),dp=eddingtonDerivative(f),dpp=8/std::pow(4-3*f*f,Real(1.5));
				Real const curvature=6*f*f/(2+std::sqrt(4-3*f*f));
				jac(0,0)=(-n1*m/(r*r)+chiPrime*f*t.energy+chi*f*t.energyDerivative-derivative[0]*t.pressureSecondDerivative)/t.pressureDerivative;
				jac(0,1)=-n1*rho/(r*r*t.pressureDerivative);jac(0,2)=chi*t.energy/t.pressureDerivative;
				jac(1,0)=r*r;
				Real const denominator=t.energy*dp;
				jac(2,0)=(-chiPrime*f*t.energy-chi*f*t.energyDerivative-d*(t.energySecondDerivative*derivative[0]+t.energyDerivative*jac(0,0))
					-curvature*t.energyDerivative/r-derivative[2]*t.energyDerivative*dp)/denominator;
				jac(2,1)=-d*t.energyDerivative*jac(0,1)/denominator;
				jac(2,2)=(-chi*t.energy-dp*t.energyDerivative*derivative[0]-d*t.energyDerivative*jac(0,2)
					-3*dp*t.energy/r-derivative[2]*t.energy*dpp)/denominator;
				dr[0]=2*n1*rho*m/(r*r*r*t.pressureDerivative);dr[1]=2*r*rho;
				dr[2]=(-d*t.energyDerivative*dr[0]+curvature*t.energy/(r*r))/denominator;
			}else{
				Real const h=std::max(x[0],Real(0));
				Real const rho=std::pow(h*(envelopeGamma_-1)/(envelopeGamma_*envelopeK_),1/(envelopeGamma_-1));
				jac(0,1)=-n1/(r*r);jac(1,0)=h>0?r*r*rho/((envelopeGamma_-1)*h):0;
				Real const f=vacuumFluxFactor(r),d=eddington(f),dp=eddingtonDerivative(f),dpp=8/std::pow(4-3*f*f,Real(1.5));
				Real const denominator=d-f*dp,a=(d-1)/denominator,ap=(dp*denominator+(d-1)*f*dpp)/(denominator*denominator);
				dr[0]=2*n1*x[1]/(r*r*r);dr[1]=2*r*rho;dr[2]=f/(r*r)*(a*a+f*ap*a-a);
			}
		};
		// The radiation-flux equation is a stiff spatial relaxation equation.
		// SDIRK2 is L-stable; step doubling controls its actual local error.
		auto stage=[&](Vector const& base,Vector const& guess,Real r,Real h){
			Vector z=guess,f(3),dr(3);Matrix j(3,3);
			for(int iteration=0;iteration<24;++iteration){
				rhs(z,f,r);jacobian(z,j,r,dr);
				Real a[3][4]{};
				for(int row=0;row<3;++row){
					for(int col=0;col<3;++col)a[row][col]=(row==col?1:0)-h*j(row,col);
					a[row][3]=base[row]+h*f[row]-z[row];
				}
				for(int k=0;k<3;++k){
					int pivot=k;for(int row=k+1;row<3;++row)if(std::abs(a[row][k])>std::abs(a[pivot][k]))pivot=row;
					for(int col=k;col<4;++col)std::swap(a[k][col],a[pivot][col]);
					if(!(std::abs(a[k][k])>0))throw std::runtime_error("Singular radiating-star radial Newton stage");
					for(int row=k+1;row<3;++row){Real const ratio=a[row][k]/a[k][k];for(int col=k;col<4;++col)a[row][col]-=ratio*a[k][col];}
				}
				Real change[3]{};for(int k=2;k>=0;--k){Real q=a[k][3];for(int col=k+1;col<3;++col)q-=a[k][col]*change[col];change[k]=q/a[k][k];}
				Real norm=0,factor=1;
				for(int k=0;k<3;++k)norm=std::max(norm,std::abs(change[k])/(p.tolerance*std::max(std::abs(z[k]),Real(1e-8))));
				if(core)for(int tries=0;tries<32;++tries){if(z[0]+factor*change[0]>0 && z[2]+factor*change[2]>0 && z[2]+factor*change[2]<1.05)break;factor*=.5;}
				for(int k=0;k<3;++k)z[k]+=factor*change[k];
				if(norm<.02)return z;
			}
			throw std::runtime_error("Radiating-star radial Newton stage did not converge");
		};
		auto step=[&](Vector const& x,Real r,Real h){
			Real const gamma=1-1/std::sqrt(Real(2));
			Vector first=stage(x,x,r+gamma*h,gamma*h),f(3),base=x;
			rhs(first,f,r+gamma*h);for(int k=0;k<3;++k)base[k]+=h*(1-gamma)*f[k];
			return stage(base,first,r+h,gamma*h);
		};
		auto fineStep=[&](Vector const& x,Real r,Real h){return step(step(x,r,h/2),r+h/2,h/2);};
		Vector x(3);for(int k=0;k<3;++k)x[k]=initial[k];
		Real r=initialRadius,h=std::min(rstart,p.maximumRadialStep);
		nodes.push_back({r,0,initial,core?coreRhs(r,initial):envelopeRhs(r,initial)});
		for(int iteration=0;iteration<p.maximumSteps;++iteration){
			Vector next(3),coarse(3);Real error=0;
			try{coarse=step(x,r,h);next=fineStep(x,r,h);}
			catch(std::runtime_error const&){h*=.5;if(h<1e-14*std::max(Real(1),r))throw;continue;}
			for(int k=0;k<3;++k)error=std::max(error,std::abs(next[k]-coarse[k])/(3*p.tolerance*std::max({std::abs(next[k]),std::abs(x[k]),Real(1e-8)})));
			if(!(error<=1)){h*=std::max(Real(.2),Real(.85)*std::pow(error,Real(-1)/3));continue;}
			bool const finished=next[0]<=target;
			if(finished){
				Real lo=0,hi=h;
				for(int n=0;n<48;++n){Real const mid=(lo+hi)/2;next=fineStep(x,r,mid);if(next[0]>target)lo=mid;else hi=mid;}
				h=(lo+hi)/2;next=fineStep(x,r,h);next[0]=target;
			}
			r+=h;x=next;std::array<Real,3> const state{x[0],x[1],x[2]};
			auto derivative=core?coreRhs(r,state):envelopeRhs(r,state);
			nodes.push_back({r,0,state,derivative});
			if(core){
				auto const t=thermodynamics(x[0]);
				Real const heat=t.energy*derivative[2]+x[2]*t.energyDerivative*derivative[0]+2*x[2]*t.energy/r;
				if(!(heat>=0) || !(x[2]<1))throw std::runtime_error("Radiating-star parameters require negative heating or nonrealizable radiation");
			}
			if(finished)return;
			h=std::min(p.maximumRadialStep,h*std::clamp(Real(.85)*std::pow(std::max(error,Real(1e-12)),Real(-1)/3),Real(.5),Real(2)));
			if(!(r<1000) || !std::isfinite(x[0]))break;
		}
		std::ostringstream message;message<<"Radiating-star "<<(core?"core":"envelope")
			<<" exceeded its radial integration budget at r="<<r<<", state="<<x[0]<<", dt="<<h;
		throw std::runtime_error(message.str());
	};
	// rho=1+A r^2+B r^4 and F=F1 r+F3 r^3 satisfy both momentum
	// equations through r^3; the mass integral includes its matching r^7 term.
	Real const a=centerDensitySecondDerivative_/2,f1=centerFluxDerivative_;
	Real const chi=opacity(1),chiPrime=p.opacityLengthScale*(1-p.cutoffDensityFraction*p.cutoffDensityFraction);
	Real const gasC=-8*(p.eos.index+1)*a/15+chiPrime*a*f1-2*tc.pressureSecondDerivative*a*a;
	Real const radC=-chiPrime*a*f1-Real(2)/3*tc.energySecondDerivative*a*a
		+Real(3.5)*f1*f1*tc.energyDerivative*a/(tc.energy*tc.energy)-Real(21)/32*std::pow(f1,4)/std::pow(tc.energy,3);
	Real const f3=(radC-tc.energyDerivative*gasC/(3*tc.pressureDerivative))
		/(chi*(1+tc.energyDerivative/(3*tc.pressureDerivative))+7*f1/tc.energy);
	Real const b=(gasC+chi*f3)/(4*tc.pressureDerivative),r2=rstart*rstart;
	Real const rhoStart=1+a*r2+b*r2*r2;
	integrate(core_,{rhoStart,rstart*r2*(Real(1)/3+a*r2/5+b*r2*r2/7),(f1+f3*r2)*rstart/thermodynamics(rhoStart).energy},rstart,p.cutoffDensityFraction,true);
	auto const& end=core_.back();transitionRadius_=end.radius;transitionFluxFactor_=end.value[2];
	if(!(transitionFluxFactor_>2*std::sqrt(Real(3))/5 && transitionFluxFactor_<1))
		throw std::runtime_error("Radiating-star opacity cutoff must reach the outward M1 vacuum branch");
	auto const te=thermodynamics(p.cutoffDensityFraction);
	luminosity_=transitionRadius_*transitionRadius_*transitionFluxFactor_*te.energy;
	envelopeGamma_=te.pressureDerivative*p.cutoffDensityFraction/te.pressure;
	envelopeK_=te.pressure/std::pow(p.cutoffDensityFraction,envelopeGamma_);
	Real const h=envelopeGamma_/(envelopeGamma_-1)*te.pressure/p.cutoffDensityFraction;
	integrate(envelope_,{h,end.value[1],transitionFluxFactor_},transitionRadius_,0,false);
	surfaceRadius_=envelope_.back().radius;mass_=envelope_.back().value[1];
	// Integrate Phi' = m/r^2 inward. Cubic Hermite mass interpolation makes the
	// Simpson interval rule independent of the hydrostatic pressure equations.
	auto fillPotential=[&](std::vector<Node>& nodes,Real outerPotential){
		nodes.back().potential=outerPotential;
		for(std::size_t i=nodes.size()-1;i>0;--i){
			auto& a=nodes[i-1];auto const& b=nodes[i];Real const h=b.radius-a.radius,rm=(a.radius+b.radius)/2;
			Real const mm=hermite(a.value[1],a.derivative[1],b.value[1],b.derivative[1],h,Real(0.5));
			Real const ga=a.radius>0?a.value[1]/(a.radius*a.radius):0,gb=b.value[1]/(b.radius*b.radius);
			a.potential=b.potential-h*(ga+4*mm/(rm*rm)+gb)/6;
		}
	};
	fillPotential(envelope_,-mass_/surfaceRadius_);fillPotential(core_,envelope_.front().potential);
	diagnostics_.radialSteps=static_cast<int>(core_.size()+envelope_.size()-2);
	diagnostics_.transitionFluxFactor=transitionFluxFactor_;
	diagnostics_.envelopeGamma=envelopeGamma_;diagnostics_.coreMassFraction=end.value[1]/mass_;
	diagnostics_.minimumHeatingCgs=std::numeric_limits<Real>::max();
	for(std::size_t i=0;i<core_.size();++i){
		auto const& a=core_[i];auto const state=sample(units::Length::from_value(a.radius*val(eos_.scaleLength())));
		if(i+1<core_.size())diagnostics_.minimumHeatingCgs=std::min(diagnostics_.minimumHeatingCgs,state.heatingCgs);
		if(i){auto const& b=core_[i-1];diagnostics_.centralOpticalDepth+=(a.radius-b.radius)*(opacity(a.value[0])*a.value[0]+opacity(b.value[0])*b.value[0])/2;}
	}
}

Real RadiatingStarAtmosphere::vacuumFluxFactor(Real r) const {
	if(r<=transitionRadius_)return transitionFluxFactor_;
	Real const d0=eddington(transitionFluxFactor_),target=std::log(r/transitionRadius_);
	Real lo=d0,hi=1;
	for(int k=0;k<60;++k){
		Real const d=(lo+hi)/2;
		if(d==1)break;
		Real const logRatio=-Real(0.5)*std::log((1-d)/(1-d0))-Real(0.25)*std::log((3*d-1)/(3*d0-1))
			+Real(0.75)*std::log((3-d)/(3-d0));
		if(logRatio<target)lo=d;else hi=d;
	}
	Real const d=(lo+hi)/2;
	return std::sqrt((3*d-1)*(3-d)/4);
}

std::array<Real,3> RadiatingStarAtmosphere::interpolate(std::vector<Node> const& nodes,Real r) const {
	auto b=std::upper_bound(nodes.begin(),nodes.end(),r,[](Real x,Node const& n){return x<n.radius;});
	if(b==nodes.begin())return b->value;
	if(b==nodes.end())return nodes.back().value;
	auto const& a=*(b-1);Real const h=b->radius-a.radius,t=(r-a.radius)/h;
	std::array<Real,3> y{};for(int k=0;k<3;++k)y[k]=hermite(a.value[k],a.derivative[k],b->value[k],b->derivative[k],h,t);
	return y;
}

std::array<Real,3> RadiatingStarAtmosphere::interpolateDerivative(std::vector<Node> const& nodes,Real r) const {
	auto b=std::upper_bound(nodes.begin(),nodes.end(),r,[](Real x,Node const& n){return x<n.radius;});
	if(b==nodes.begin())return b->derivative;
	if(b==nodes.end())return nodes.back().derivative;
	auto const& a=*(b-1);Real const h=b->radius-a.radius,t=(r-a.radius)/h;
	std::array<Real,3> derivative{};
	for(int k=0;k<3;++k)derivative[k]=6*t*(1-t)*(b->value[k]-a.value[k])/h
		+(1-4*t+3*t*t)*a.derivative[k]+(3*t*t-2*t)*b->derivative[k];
	return derivative;
}

Real RadiatingStarAtmosphere::potential(std::vector<Node> const& nodes,Real r) const {
	auto b=std::upper_bound(nodes.begin(),nodes.end(),r,[](Real x,Node const& n){return x<n.radius;});
	if(b==nodes.begin())return b->potential;
	if(b==nodes.end())return nodes.back().potential;
	auto const& a=*(b-1);Real const h=b->radius-a.radius,t=(r-a.radius)/h;
	return hermite(a.potential,a.radius>0?a.value[1]/(a.radius*a.radius):0,b->potential,b->value[1]/(b->radius*b->radius),h,t);
}

RadiatingStarAtmosphere::State RadiatingStarAtmosphere::sample(units::Length radius) const {
	Real const length=val(eos_.scaleLength()),r=val(radius)/length,rhoc=val(parameters_.eos.centralDensity),pc=val(eos_.centralPressure()),c=val(constants::c);
	if(!std::isfinite(r) || r<0)throw std::invalid_argument("Radiating-star sample radius must be finite and nonnegative");
	Real rho=0,m=mass_,f=0,pg=0,temp=0,energy=0,flux=0,rhop=0,pgp=0,ep=0,fp=0,heat=0,kappa=0,phi=0;
	if(r<transitionRadius_){
		auto y=interpolate(core_,r);rho=y[0];m=y[1];f=y[2];auto const t=thermodynamics(rho);
		temp=t.temperature;energy=t.energy;pg=t.pressure;flux=f*energy;kappa=opacity(rho);
		// Derive the fixed heater from the same continuous interpolant used for
		// initialization. Re-evaluating the stiff radial RHS here would amplify
		// its small integration error and spoil the interpolated energy balance.
		if(r>0){auto const derivative=interpolateDerivative(core_,r);rhop=derivative[0];ep=t.energyDerivative*rhop;pgp=t.pressureDerivative*rhop;fp=energy*derivative[2]+f*ep;heat=fp+2*flux/r;}
		else {fp=centerFluxDerivative_;heat=3*centerFluxDerivative_;}
		phi=potential(core_,r);
	}else{
		f=vacuumFluxFactor(r);flux=luminosity_/(r*r);energy=flux/f;fp=-2*flux/r;
		Real const d=eddington(f),dp=eddingtonDerivative(f),fprime=f/r*(d-1)/(d-f*dp);ep=-2*energy/r-energy*fprime/f;
		if(r<surfaceRadius_){
			auto const y=interpolate(envelope_,r);m=y[1];Real const h=std::max(y[0],Real(0));
			rho=std::pow(h*(envelopeGamma_-1)/(envelopeGamma_*envelopeK_),1/(envelopeGamma_-1));pg=envelopeK_*std::pow(rho,envelopeGamma_);
			Real const gasConstant=val(constants::boltzmann/constants::atomicMassUnit)/parameters_.eos.meanMolecularWeight;
			temp=pc*pg/(rhoc*rho*gasConstant);pgp=-(parameters_.eos.index+1)*rho*m/(r*r);rhop=pgp/(envelopeGamma_*pg/rho);
			phi=potential(envelope_,r);
		}else phi=-mass_/r;
	}
	Real const d=eddington(f),massScale=4*pi*rhoc*length*length*length;
	State result;
	result.density=units::Density::from_value(rho*rhoc);result.temperature=units::Temperature::from_value(temp);
	result.gasPressure=units::Pressure::from_value(pg*pc);result.radiationEnergy=units::EnergyDensity::from_value(energy*pc);
	result.radiationFlux=units::EnergyFlux::from_value(flux*c*pc);result.enclosedMass=units::Mass::from_value(m*massScale);
	result.radiationPressureRadial=units::Pressure::from_value(d*energy*pc);result.radiationPressureTangential=units::Pressure::from_value((1-d)*energy*pc/2);
	result.potential=units::VelocitySquared::from_value(phi*val(eos_.centralIntegralH()));
	result.gravityMagnitude=units::Acceleration::from_value(r>0?val(eos_.centralIntegralH())/length*m/(r*r):0);
	result.densityDerivative=rhop*rhoc/length;result.gasPressureDerivative=pgp*pc/length;result.radiationEnergyDerivative=ep*pc/length;result.radiationFluxDerivative=fp*c*pc/length;
	result.opacityCgs=kappa/(rhoc*length);result.heatingCgs=heat*c*pc/length;result.fluxFactor=f;
	return result;
}

units::Length RadiatingStarAtmosphere::transitionRadius() const{return units::Length::from_value(transitionRadius_*val(eos_.scaleLength()));}
units::Length RadiatingStarAtmosphere::surfaceRadius() const{return units::Length::from_value(surfaceRadius_*val(eos_.scaleLength()));}
units::Mass RadiatingStarAtmosphere::mass() const{Real const a=val(eos_.scaleLength());return units::Mass::from_value(4*pi*val(parameters_.eos.centralDensity)*a*a*a*mass_);}
Real RadiatingStarAtmosphere::luminosityCgs() const{Real const a=val(eos_.scaleLength());return 4*pi*a*a*val(constants::c)*val(eos_.centralPressure())*luminosity_;}

} // namespace octotigerII::problems
