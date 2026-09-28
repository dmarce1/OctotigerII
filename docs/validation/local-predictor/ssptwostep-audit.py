"""Audit of the user-supplied variable-step SSP two-step RK draft.

Standalone diagnostics only; this is not a production integrator.
Uses printed closed-form coefficients and its fitted c4 exactly as supplied.
The fit is not assumed valid outside the audited step-ratio interval.
"""
import math

def coefficients(a, c4=None):
    if c4 is None:
        d=a-.47732381
        c4=-3.41125065*d**.5+8.16487469*d**(1/3)-4.01095172*d**.2-.001
    r=math.sqrt(3*a**9*c4**3*(2+8*a+12*a*a+71*a**5*c4+16*a**6*c4*c4+2*a**4*(1+71*c4)+a**3*(8+71*c4)))
    s3=1+15*a*a+408*a**8*c4*c4+64*a**9*c4**3+4*a**3*(5+12*c4)+48*a**7*c4*(1+17*c4)+6*a**5*(1+48*c4)+3*a**4*(5+64*c4)+a**6*(1+192*c4+408*c4*c4)+a*(6-24*r)-24*r
    s=math.cbrt(s3)
    v=(1+32*a**5*c4+16*a**6*c4*c4+a**4*(1+64*c4)+a**3*(4-8*c4*(-4+s))+s+s*s+2*a*(2+s)+a*a*(6+s))/(24*a*a*c4*s)
    c3=((1+a)**2*(4*a*(-1+v)*v-2*v*v+a*a*(-1+2*v)))/(2*a**3*v*(a+2*v))
    v5=-.5*(a+v+3*a*a*v*(1-2*c4*v)-a**3*(1+2*(-2+2*c3+c4)*v))/(a*a*(-1+c3+c4)*(a+3*v))
    c1=(1-a**3*(-1+c3)-3*a*a*(c4*(v-v5)+v5-c3*v5))/(a*a*c4*(a-3*v))
    return c1,1-c1,c3,c4,1-c3-c4,v,v5


def order_residuals(a):
    c1,c2,c3,c4,c5,v,v5=coefficients(a)
    stage_time=c1*v-a*c2
    stage_second=a*a*c2/2
    old_third=-a**3/6
    shared=c5*(old_third+v5*a*a/2)
    return (
        c3+c4+c5-1,
        c4*(stage_time+v)+c5*(-a+v5)-1,
        c4*(stage_second+v*stage_time)+c5*(a*a/2-a*v5)-.5,
        c4*(c2*old_third+v*stage_time**2/2)+shared-1/6,
        c4*(c2*old_third+v*stage_second)+shared-1/6,
    )


def nonlinear_error(steps, variable):
    weights=[(.8 if i%2==0 else 1.2) if variable else 1.0 for i in range(steps)]
    widths=[.5*w/sum(weights) for w in weights]
    h_previous=widths[0]
    # Exact startup isolates the multistep method's error; y'=y^2, y(0)=1.
    old,current=1.,1/(1-h_previous)
    for h in widths[1:]:
        c1,c2,c3,c4,c5,v,v5=coefficients(h_previous/h)
        stage=c1*(current+v*h*current**2)+c2*old
        old,current=current,c3*current+c4*(stage+v*h*stage**2)+c5*(old+v5*h*old**2)
        h_previous=h
    return abs(current-2.)


def report(label, errors):
    print(label)
    print('  errors:', ' '.join(f'{e:.9g}' for e in errors))
    print('  orders:', ' '.join(f'{math.log2(a/b):.6f}' for a,b in zip(errors,errors[1:])))


if __name__=='__main__':
    import cmath
    for a,c4 in ((.5,.5),(.75,.5),(.5,.75),(.75,.75),(1,.75)):
        print('draft table alpha,c4=',a,c4,'(c1,c2,c3,c4,c5,v,v5)=',coefficients(a,c4))
    ratios=[i/1000 for i in range(700,2001)]
    print('sampled alpha interval [.7,2], spacing .001:')
    print('  minimum coefficient:',min(min(coefficients(a)) for a in ratios))
    print('  maximum third-order condition residual:',max(abs(r) for a in ratios for r in order_residuals(a)))
    print('fit at alpha=10:',coefficients(10))
    for variable in (False,True):
        report('nonlinear ODE, variable steps='+str(variable),[nonlinear_error(n,variable) for n in (20,40,80,160)])
    cells,mode,stop=32,4,.25
    dx=1/cells; theta=2*math.pi*mode/cells
    spatial=-(1-cmath.exp(-1j*theta))*(1+.5j*math.sin(theta))/dx
    c1,c2,c3,c4,c5,v,v5=coefficients(1)
    errors=[]
    for steps in (128,256,512,1024):
        h=stop/steps; z=h*spatial
        old,current=1+0j,1+z+z*z/2+z**3/6  # SSPRK3 startup.
        for _ in range(1,steps):
            stage=c1*(current+v*z*current)+c2*old
            old,current=current,c3*current+c4*(stage+v*z*stage)+c5*(old+v5*z*old)
        errors.append(abs(current-cmath.exp(stop*spatial)))
    report('fixed-mesh centered PLM/upwind advection, exact spatial exponential reference',errors)
    print('alpha=1 SSP coefficient:',1/max(v,v5))
    print('alpha=1 stage time fraction:',c1*v-c2)
    print('RHS-only speed ratio vs SSPRK3 at SSP limits:',3/(2*max(v,v5)))
