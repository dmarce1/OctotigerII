"""Spherical full-M1 reference used by the rotating-star table generator.

Units are 4*pi*G = R_gas*T_c = rho_c = c = 1.  Fixed model parameters
are beta_c=0.8 and structural polytropic index n=3.5.  The c=1 convention
here scales radiation flux by the physical light speed; this is a stationary
reference calculation and does not make a reduced-speed approximation.

The optically coupled core has p_total=1.25*rho**(9/7).  Its opacity tapers
to exactly zero at rho=0.015.  A continuously matched transparent gas
polytrope and the outgoing vacuum M1 branch complete the finite-radius star.
"""
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
N=3.5; Gam=1+1/N; K=1.25; A=.25
def eos(rho):
    rho=max(rho,1.e-30); p=K*rho**Gam
    T=brentq(lambda x: rho*x+A*x**4-p,0,(p/A)**.25,xtol=1e-25)
    Tr=(K*Gam*rho**(Gam-1)-T)/(rho+4*A*T**3)
    E=3*A*T**4; Er=12*A*T**3*Tr; pg=rho*T; pgr=T+rho*Tr
    return T,E,Er,pg,pgr

def run(k0,rcut,rtol=2e-8):
    def kap(rho): return k0*max(0,1-rcut/max(rho,1e-30))**2
    T,E,Er,pg,pgr=eos(1.); pgE=pgr/Er; chi=kap(1.)
    a=2.5/E; b=chi*(1+1/(3*pgE)); d=1/(9*pgE)
    f1=2*d/(b+np.sqrt(b*b+4*a*d))
    rho2=(-1/3+chi*f1)/pgr
    r0=1e-5; y0=[1+.5*rho2*r0*r0,r0**3/3,f1*r0/E]
    def deriv(r,y):
        rho,m,f=y; T,E,Er,pg,pgr=eos(rho); chi=kap(rho)*rho; F=f*E
        D=1/3+2*f*f/(2+np.sqrt(max(1e-20,4-3*f*f)))
        Dp=2*f/np.sqrt(max(1e-20,4-3*f*f))
        rhop=(-rho*m/(r*r)+chi*F)/pgr
        fp=(-chi*F-D*Er*rhop-(3*D-1)*E/r)/(E*Dp)
        return [rhop,r*r*rho,fp]
    def end(r,y):return y[0]-rcut
    end.terminal=True;end.direction=-1
    def cone(r,y):return y[2]-.999
    cone.terminal=True;cone.direction=1
    s=solve_ivp(deriv,[r0,100],y0,method='Radau',rtol=rtol,dense_output=True,atol=[1e-12,1e-12,1e-12],events=[end,cone],max_step=.02)
    js=[]
    for r,y in zip(s.t,s.y.T):
        rho,m,f=y; T,E,Er,pg,pgr=eos(rho); rr,mm,ff=deriv(r,y)
        js.append(E*ff+f*Er*rr+2*f*E/r)
    print('k,cut',k0,rcut,'success',s.success,'R,rho,f',s.t[-1],s.y[0,-1],s.y[2,-1],'jmin/max',min(js),max(js),'steps',len(js),'why',s.message)
    return s,np.array(js),deriv


class PiecewiseRadial:
    """PCHIP branches with a shared value and independent derivatives at joins."""
    def __init__(self, radii, values, joins):
        edges=[radii[0],*joins,radii[-1]]
        self.joins=joins
        self.branches=[PchipInterpolator(radii[(radii>=left)&(radii<=right)], values[(radii>=left)&(radii<=right)]) for left,right in zip(edges[:-1],edges[1:])]
    def __call__(self,r):
        r=np.asarray(r)
        answer=self.branches[-1](r)
        for i in range(len(self.joins)-1,-1,-1): answer=np.where(r<=self.joins[i],self.branches[i](r),answer)
        return answer
    def derivative(self,nu=1):
        answer=object.__new__(PiecewiseRadial);answer.joins=self.joins;answer.branches=[b.derivative(nu) for b in self.branches]
        return answer

def spherical(kappa=100., cut=.015):
    core,_,rhs=run(kappa,cut,rtol=1e-10)
    rc=core.t[-1]; rhoc,mc,fc=core.y[:,-1]
    tc,ec,_,pc,pcr=eos(rhoc)
    gam=pcr*rhoc/pc; kg=pc/rhoc**gam
    hc=gam/(gam-1)*tc; lum=rc*rc*fc*ec
    def outer(r,y):
        h,m,f=y
        rho=(max(h,0)*(gam-1)/(gam*kg))**(1/(gam-1))
        root=np.sqrt(max(1e-14,4-3*f*f)); d=1/3+2*f*f/(2+root); dp=2*f/root
        return [-m/r**2,r*r*rho,f/r*(d-1)/(d-f*dp)]
    def surface(r,y):return y[0]
    surface.terminal=True;surface.direction=-1
    envelope=solve_ivp(outer,[rc,100],[hc,mc,fc],events=surface,rtol=1e-11,atol=1e-12,dense_output=True,max_step=.04)
    rs=envelope.t[-1]; ms=envelope.y[1,-1]
    vacuum=solve_ivp(outer,[rs,4*rs],envelope.y[:,-1],rtol=1e-11,atol=1e-12,dense_output=True,max_step=.1)
    rg=np.unique(np.r_[np.linspace(0,rc,6001),np.linspace(rc,rs,4001),np.linspace(rs,4*rs,4001)])
    inner=rg<rc; env=(rg>=rc)&(rg<rs); vac=rg>=rs
    rho=np.zeros_like(rg); temp=rho.copy(); e=rho.copy(); f=rho.copy(); m=rho.copy()
    rhoi,mi,fi=core.sol(np.maximum(rg[inner],core.t[0])); rhoi[0]=1;mi[0]=0;fi[0]=0
    rho[inner]=rhoi;m[inner]=mi;f[inner]=fi
    vals=np.array([eos(r) for r in rhoi]);temp[inner]=vals[:,0];e[inner]=vals[:,1]
    he,me,fe=envelope.sol(rg[env]);rho[env]=(he*(gam-1)/(gam*kg))**(1/(gam-1));temp[env]=kg*rho[env]**(gam-1);m[env]=me;f[env]=fe
    _,mv,fv=vacuum.sol(rg[vac]);m[vac]=ms;f[vac]=fv
    e[~inner]=lum/(rg[~inner]**2*f[~inner])
    force=np.divide(m,rg*rg,out=np.zeros_like(m),where=rg>0)
    integ=cumulative_trapezoid(force,rg,initial=0)
    phi=-ms/rg[-1]-(integ[-1]-integ)
    # Barotropic gas enthalpy, continuously joined to the transparent envelope.
    densities=np.geomspace(cut,2,6000); es=np.array([eos(x) for x in densities])
    integrals=hc+cumulative_trapezoid(es[:,4]/densities,densities,initial=0)
    h_of_rho=PchipInterpolator(densities,integrals)
    h=np.zeros_like(rg);h[inner]=h_of_rho(rho[inner]);h[env]=he;h[vac]=-ms/rs+ms/rg[vac]
    # EOS on the whole h interval, including vacuum.
    denlow=np.geomspace(1e-18,cut,2500)[:-1]
    den=np.r_[0,denlow,densities]
    htab=np.r_[0,gam/(gam-1)*kg*denlow**(gam-1),integrals]
    ttab=np.r_[0,kg*denlow**(gam-1),es[:,0]]
    info=dict(R=rs,M=ms,rc=rc,cut=cut,kappa=kappa,gammaEnvelope=gam,lum=lum,hc=h[0])
    interp={key:PiecewiseRadial(rg,value,[rc,rs]) for key,value in dict(rho=rho,h=h,E=e,f=f,phi=phi,mass=m).items()}
    interp['density']=PchipInterpolator(htab,den,extrapolate=True)
    interp['temperature']=PchipInterpolator(htab,ttab,extrapolate=True)
    return info,interp
