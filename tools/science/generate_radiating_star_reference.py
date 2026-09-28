"""Offline rotating gas/gravity/M1 stellar-reference generator.

The reference solve uses F_h(U)-F_h(U_spherical).  This is an analytic-
background discretization, not a force or source added to the evolution.
Unknowns are h/h_c, log(E/E_spherical), f_r, q, Phi/h_c, with
f_theta=-sqrt(1-mu**2)*q.  The gas enthalpy h belongs to the prescribed
gas barotrope.  Equatorial rotation is imposed; the rotation elsewhere is
recovered from radial force balance and checked for nonnegative Omega².

The radial mesh is graded on both sides of the fixed opacity cutoff.
Outward radial pressure derivatives suppress the decoupled collocation
modes.  The exact spherical solution, core/exterior interpolation split,
and explicit continuum-refinement checks are part of the construction.
"""
import argparse,time
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.sparse import lil_matrix,diags,kron,eye,bmat,csr_matrix
from scipy.sparse.linalg import splu
from scipy.special import eval_legendre
from scipy.optimize._numdiff import approx_derivative,group_columns
from radiating_star_spherical import spherical

class SphericalModel:
 def __init__(self,nr=49,nm=17,spin=0,extent=2.):
  self.info,self.ref=spherical();self.nr=nr;self.nm=nm;self.N=nr*nm
  nc=(nr-1)//2; t=np.linspace(0,1,nc+1); self.x=np.r_[self.info['rc']*(2*t-t*t),self.info['rc']+(extent-1)*self.info['rc']*t[1:]**2];self.mu=np.linspace(0,1,nm)
  self.dr=self.x[1];self.dm=self.mu[1];self.r,self.MU=np.meshgrid(self.x,self.mu,indexing='ij');self.S=1-self.MU**2
  self.R=self.r*np.sqrt(self.S);self.Z=self.r*self.MU
  self.omega2=spin**2*self.info['M']/self.info['R']**3
  self.hscale=self.info['hc'];self.Eref=self.ref['E'](self.r);self.forceScale=np.maximum(.02,self.info['M']/(1+self.r**2));self.radScale=np.maximum(1e-7,self.Eref/np.maximum(self.r,1))
  self.initial=np.array([self.ref['h'](self.r)/self.hscale,np.zeros_like(self.r),self.ref['f'](self.r),np.zeros_like(self.r),self.ref['phi'](self.r)/self.hscale])
  self.mask=self.r<=self.info['rc']+1e-10
  self.kappa=self.info['kappa']*np.maximum(0,1-self.info['cut']/np.maximum(self.ref['rho'](self.r),1e-30))**2
  self.kappa[~self.mask]=0
  self.rho0=self.decode(self.initial.ravel())[2];self.phib=self.ref['phi'](self.r[-1]);self.phib0=self.phib.copy()
  self.ops();self.baseline=self.rawresidual(self.initial.ravel())
 def decode(self,y):
  h,l,f,q,p=y.reshape(5,self.nr,self.nm);h=h*self.hscale;p=p*self.hscale
  rho=np.maximum(self.ref['density'](np.maximum(h,0)),0);T=np.maximum(self.ref['temperature'](np.maximum(h,0)),0)
  E=self.Eref*np.exp(l)
  return h,p,rho,T,E,E*f,E*q
 def ops(self):
  nr,nm=self.nr,self.nm;N=self.N
  def central(n,dx,odd=False):
   a=lil_matrix((n,n))
   for i in range(1,n-1):a[i,i-1]=-.5/dx;a[i,i+1]=.5/dx
   a[-1,-3]=.5/dx;a[-1,-2]=-2/dx;a[-1,-1]=1.5/dx
   if odd:a[0,1]=1/dx
   return a.tocsr()
  def backward(n,dx):
   a=lil_matrix((n,n));a[1,0]=-2/dx;a[1,1]=2/dx
   for i in range(2,n):a[i,i]=1.5/dx;a[i,i-1]=-2/dx;a[i,i-2]=.5/dx
   return a.tocsr()
  self.re=kron(central(nr,self.dr),eye(nm),format='csr');self.ro=kron(central(nr,self.dr,True),eye(nm),format='csr')
  self.me=kron(eye(nr),central(nm,self.dm),format='csr');self.mo=kron(eye(nr),central(nm,self.dm,True),format='csr')
  self.rb=kron(backward(nr,self.dr),eye(nm),format='csr');self.mb=kron(eye(nr),backward(nm,self.dm),format='csr')
  ir=np.divide(1,self.r.ravel(),out=np.zeros(N),where=self.r.ravel()>0);self.ir=diags(ir);self.sm=diags(self.S.ravel());self.mm=diags(self.MU.ravel())
  def fdweights(nodes, at, degree=1):
   offsets=nodes-at
   a=np.array([offsets**k for k in range(len(nodes))]);rhs=np.zeros(len(nodes));rhs[degree]=1 if degree==1 else 2
   return np.linalg.solve(a,rhs)
  self.weights=fdweights
  def radial_derivative(back=False,odd=False):
   a=lil_matrix((nr,nr))
   for i in range(1,nr):
    if back and i==1:a[i,0]=-2/self.x[1];a[i,1]=2/self.x[1];continue
    sites=np.arange(i-2,i+1) if back or i==nr-1 else np.arange(i-1,i+2)
    a[i,sites]=fdweights(self.x[sites],self.x[i])
   if odd:a[0,1]=1/self.x[1]
   return kron(a.tocsr(),eye(nm),format='csr')
  self.re=radial_derivative();self.ro=radial_derivative(odd=True);self.rb=radial_derivative(back=True)
  self.gas=self.mm@self.rb+self.ir@self.sm@self.mb
  integ=lil_matrix((nr,nr))
  for i in range(1,nr):
   integ[i,0]=(self.x[1]-self.x[0])/2;integ[i,i]=(self.x[i]-self.x[i-1])/2
   for j in range(1,i):integ[i,j]=(self.x[j+1]-self.x[j-1])/2
  select=lil_matrix((nr,N));select[np.arange(nr),np.arange(nr)*nm]=1
  self.eqint=integ.tocsr()@select.tocsr()
  lap=lil_matrix((N,N));self.center=np.arange(nm);self.equator=np.arange(nr)*nm;self.outer=np.arange((nr-1)*nm,N)
  weights=np.full(nm,self.dm);weights[[0,-1]]*=.5
  for i in range(nr-1):
   for j in range(nm):
    k=i*nm+j;r=self.x[i];mu=self.mu[j]
    if i==0:
     if j==0:
      lap[k,k]=-6/self.dr**2
      for z in range(nm):lap[k,nm+z]=6*weights[z]/self.dr**2
     else:lap[k,k]=1;lap[k,0]=-1
     continue
    sites=np.arange(i-1,i+2);wr=fdweights(self.x[sites],r,2)+2/r*fdweights(self.x[sites],r);lap[k,sites*nm+j]=wr
    if j==0:lap[k,k]-=2/(r*r*self.dm**2);lap[k,k+1]+=2/(r*r*self.dm**2)
    elif j==nm-1:
     lap[k,k]+=-3/(r*r*self.dm);lap[k,k-1]+=4/(r*r*self.dm);lap[k,k-2]+=-1/(r*r*self.dm)
    else:
     lap[k,k]-=2*(1-mu*mu)/(r*r*self.dm**2);lap[k,k-1]+=(1-mu*mu)/(r*r*self.dm**2)+mu/(r*r*self.dm);lap[k,k+1]+=(1-mu*mu)/(r*r*self.dm**2)-mu/(r*r*self.dm)
  self.lap=lap.tocsr()
 def rawresidual(self,y):
  h,p,rho,T,E,F,Q=[a.ravel() for a in self.decode(y)];kap=self.kappa.ravel();r=self.r.ravel();mu=self.MU.ravel();ss=self.S.ravel();ir=self.ir.diagonal();fs=self.forceScale.ravel();rs=self.radScale.ravel();S=h+p
  f=F/E;q=Q/E;f2=f*f+ss*q*q;b=3/(2+np.sqrt(4-3*f2));A=E*(1/3-b*f2/3);P=A+E*b*f*f;C=E*b*f*q;D=E*b*q*q
  gas=(self.gas@S-kap*(mu*F+ss*Q))/fs
  eq=self.equator;gas[eq]=(S[eq]-self.hscale-p[0]-.5*self.omega2*self.x**2-self.eqint@(kap*F))/self.hscale
  gas[self.center]=h[self.center]/self.hscale-1
  rr=(self.rb@P+ir*(2*P-2*A-ss*D+ss*(self.mo@C)-2*mu*C)+kap*rho*F)/rs
  rz=(self.rb@C+ir*(3*C+self.me@A+ss*(self.me@D)-3*mu*D)+kap*rho*Q)/rs
  rr[self.center]=f[self.center];rz[eq]=q[eq];rz[self.center]=q[self.center]
  divf=self.rb@F+ir*(2*F+ss*(self.mo@Q)-2*mu*Q)
  thermal=np.where(self.mask.ravel(),(E-.75*T**4)/np.maximum(E,1e-8),divf/rs)
  grav=self.lap@p-rho;grav[self.center[1:]]=(p[self.center[1:]]-p[0])/self.hscale;grav[self.outer]=(p[self.outer]-self.phib)/self.hscale
  return np.array([gas,thermal,rr,rz,grav]).ravel()
 def residual(self,y):return self.rawresidual(y)-self.baseline
 def update_boundary(self,y):
  rho=self.decode(y)[2]-self.rho0;self.phib=self.phib0.copy()
  for ell in range(0,min(16,self.nm),2):
   pl=eval_legendre(ell,self.mu);proj=np.trapz(rho*pl,self.mu,axis=1);moment=np.trapz(proj*self.x**(ell+2),self.x)
   self.phib-=moment/self.x[-1]**(ell+1)*pl
 def pattern(self):
  N=self.N;nm=self.nm;nr=self.nr
  stencil=(abs(self.re)+abs(self.ro)+abs(self.me)+abs(self.mo)+abs(self.gas)+abs(self.lap)+eye(N)).astype(bool).tolil()
  for i in range(nr):stencil[i*nm,np.arange(i+1)*nm]=1;stencil[i*nm,0]=1
  stencil[self.center,:]=0
  for i in self.center:stencil[i,i]=1;stencil[i,0]=1
  stencil[0,nm:2*nm]=1
  return bmat([[stencil]*5]*5,format='csr')

def cone_step(m,y,d):
 f=y.reshape(5,m.nr,m.nm)[2:4];df=d.reshape(5,m.nr,m.nm)[2:4];weights=np.array([np.ones_like(m.S),m.S]);a=np.sum(weights*df*df,axis=0);b=2*np.sum(weights*f*df,axis=0);c=np.sum(weights*f*f,axis=0)-(1-1e-9)**2
 disc=np.sqrt(np.maximum(0,b*b-4*a*c));roots=np.full_like(a,np.inf);active=a>0;pos=active&(b>=0);neg=active&(b<0);roots[pos]=-2*c[pos]/(b[pos]+disc[pos]);roots[neg]=(-b[neg]+disc[neg])/(2*a[neg])
 return min(1.,.99*np.min(roots),2/max(2,np.max(np.abs(d.reshape(5,m.nr,m.nm)[1]))))

def solve(m,y,iterations=30):
 pattern=m.pattern();groups=group_columns(pattern);print('COLORS',max(groups)+1,flush=True)
 for it in range(iterations):
  t=time.monotonic();f=m.residual(y);norm=np.linalg.norm(f);print('NEWTON',it,norm,max(abs(f)),flush=True)
  if max(abs(f))<2e-9:return y,True
  jac=approx_derivative(m.residual,y,method='3-point',sparsity=(pattern,groups))
  d=splu(jac.tocsc()).solve(-f);alpha=cone_step(m,y,d)
  accepted=False
  for line in range(15):
   trial=y+alpha*d;fn=m.residual(trial);nn=np.linalg.norm(fn)
   if np.isfinite(nn) and nn<(1-1e-4*alpha)*norm:y=trial;accepted=True;break
   alpha*=.5
  print('STEP',alpha,'accepted',accepted,'sec',time.monotonic()-t,flush=True)
  if not accepted:return y,False
 return y,False


def export_include(model, state, toroidal, spin, c_dimensionless, path):
    """Write the agreed C++ reader format as one raw string literal.

    The cutoff row has identical values in both pieces, but independent
    radial derivatives. Values are corrections to the spherical background,
    except W, whose spherical value is zero. h and Phi use physical
    dimensionless units rather than the Newton solver's h_c scaling.
    """
    import io
    from pathlib import Path
    from scipy.interpolate import RectBivariateSpline

    delta = state.reshape(5, model.nr, model.nm) - model.initial
    delta[[0, 4]] *= model.hscale
    values = np.concatenate((delta, toroidal[None]), axis=0)
    cut = (model.nr - 1) // 2
    pieces = [slice(0, cut + 1), slice(cut, model.nr)]
    output = io.StringIO()
    output.write('OCTOII_RADIATING_STAR 1\n')
    output.write(' '.join(format(v, '.17e') for v in
                          (3.5, 0.8, 100., .015, spin,
                           model.x[-1]/model.info['rc'], c_dimensionless)) + '\n')
    output.write(f'{cut+1} {model.nr-cut} {model.nm} 6\n')
    for axis in (model.mu, model.x[:cut+1], model.x[cut:]):
        output.write(' '.join(format(v, '.17e') for v in axis) + '\n')
    for piece in pieces:
        radius = model.x[piece]
        jets = np.empty((len(radius), model.nm, 6, 4))
        for field in range(6):
            spline = RectBivariateSpline(radius, model.mu, values[field, piece],
                                         kx=3, ky=3, s=0)
            for slot, (dr, dmu) in enumerate(((0, 0), (1, 0), (0, 1), (1, 1))):
                jets[:, :, field, slot] = spline(radius, model.mu, dx=dr, dy=dmu)
            jets[:, :, field, 0] = values[field, piece]
            # Equatorial parity is exact, rather than a spline extrapolation.
            if field != 3:
                jets[:, 0, field, 2:] = 0
            else:
                jets[:, 0, field, :2] = 0
            if radius[0] == 0 and field in (0, 1, 4):
                jets[0, :, field, 1] = 0
                jets[0, :, field, 2:] = 0
        np.savetxt(output, jets.reshape(-1, 24), fmt='%.17e')
    Path(path).write_text('R"OCTOIIRADSTAR(\n' + output.getvalue()
                          + ')OCTOIIRADSTAR"\n')


def main():
    """Reproduce the nonlinear solve, moving correction, and table export."""
    from pathlib import Path
    from radiating_star_rotation import (
        angular_velocity_squared, light_speed, moving_residual, toroidal_flux,
    )

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nr', type=int, default=193)
    parser.add_argument('--nm', type=int, default=65)
    parser.add_argument('--extent', type=float, default=4.)
    parser.add_argument('--spin', type=float, default=.1)
    parser.add_argument('--outer', type=int, default=5)
    parser.add_argument('--iterations', type=int, default=25)
    parser.add_argument('--resume', type=Path)
    parser.add_argument('--output', type=Path, required=True,
                        help='NumPy state and metadata for independent audits')
    parser.add_argument('--include', type=Path,
                        help='Optional generated C++ raw-text include')
    parser.add_argument('--static-only', action='store_true',
                        help='Diagnostic meridional solve without moving corrections')
    args = parser.parse_args()
    if args.nr < 17 or args.nr % 2 == 0 or args.nm < 5 or args.extent <= 1:
        parser.error('Use an odd radial count >=17, angular count >=5, extent >1')
    if args.spin < 0 or args.spin > .2:
        parser.error('This validated continuation family currently covers 0 <= spin <= .2')
    model = SphericalModel(args.nr, args.nm, 0, args.extent)
    state = model.initial.ravel().copy()
    if args.resume:
        saved = np.load(args.resume)
        if saved['solution'].shape != model.initial.shape:
            raise ValueError('The resumed state has a different grid shape')
        if not np.allclose(saved['x'], model.x, rtol=0, atol=1e-10):
            raise ValueError('The resumed state has a different radial grid')
        state = saved['solution'].ravel().copy()
    model.omega2 = args.spin**2 * model.info['M'] / model.info['R']**3
    for iteration in range(args.outer):
        model.update_boundary(state)
        state, converged = solve(model, state, args.iterations)
        if not converged:
            raise RuntimeError('The stationary meridional reference did not converge')
    c_dimensionless = light_speed()
    toroidal = np.zeros((model.nr, model.nm))
    if not args.static_only:
        for iteration in range(4):
            omega2 = angular_velocity_squared(model, state, toroidal, c_dimensionless)
            toroidal = toroidal_flux(model, state, omega2, c_dimensionless, toroidal)
            beta = model.R * np.sqrt(omega2) / c_dimensionless
            model.residual = lambda trial: moving_residual(model, trial, toroidal, beta)
            state, converged = solve(model, state, args.iterations)
            if not converged:
                raise RuntimeError('The retained moving-gas reference did not converge')
            model.update_boundary(state)
    omega2 = angular_velocity_squared(model, state, toroidal, c_dimensionless)
    energy = model.decode(state)[4]
    f2 = state.reshape(5, model.nr, model.nm)[2]**2 + model.S * (
        state.reshape(5, model.nr, model.nm)[3]**2 + (toroidal/energy)**2)
    if np.any(f2 >= 1):
        raise RuntimeError('The converged reference violates radiation realizability')
    residual = model.residual(state).reshape(5, model.nr, model.nm)
    np.savez(args.output, solution=state.reshape(model.initial.shape), W=toroidal,
             omega2=omega2, x=model.x, mu=model.mu, phib=model.phib,
             residual=residual, spin=args.spin, cDimensionless=c_dimensionless,
             extent=args.extent, beta=.8, index=3.5, kappa=100., cut=.015,
             movingTerms=not args.static_only)
    if args.include:
        if args.static_only:
            raise ValueError('A shipped include requires the retained moving-gas terms')
        export_include(model, state, toroidal, args.spin, c_dimensionless, args.include)
    print('FINAL', float(np.max(np.abs(residual))), str(args.output), flush=True)


if __name__ == '__main__':
    main()
