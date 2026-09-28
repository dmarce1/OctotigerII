#!/usr/problem/env python3
"""Compare with untouched Timmes routines from explicitly supplied source directories.
Requires gfortran. No downloads. Does not modify the reference sources.
"""
import argparse
import math
from pathlib import Path
import random
import re
import subprocess
import tempfile

p=argparse.ArgumentParser()
p.add_argument('--helmholtz-source',type=Path,required=True)
p.add_argument('--timmes-source',type=Path,required=True)
p.add_argument('--probe',type=Path,required=True)
a=p.parse_args()
header=Path(__file__).parents[1]/'include/octotigerII/helmholtz/helmholtz.hpp'
fields=re.findall(r'double (\w+)\{\}',header.read_text())

def compile_reference(work,source,routine,outputs):
    src=(source/('helmholtz.f90' if routine=='helmeos' else 'eosfxt.f90')).read_text()
    src=src[src.index('      subroutine '):]
    driver="""program compare
include 'implno.dek'
include 'vector_eos.dek'
integer ios
jlo_eos=1
jhi_eos=1
"""
    if routine=='helmeos': driver+='call read_helm_table\n'
    driver+='do\nread(*,*,iostat=ios) den_row(1),temp_row(1),abar_row(1),zbar_row(1)\nif(ios/=0) exit\n'
    driver+='call '+routine+'\nif(eosfail) stop 2\n'
    driver+="write(*,'(*(ES26.17E3,1X))') "+', &\n'.join(outputs)+'\nend do\nend program\n'
    f=work/(routine+'.f90'); f.write_text(driver+src)
    exe=work/routine
    subprocess.run(['gfortran','-O2','-ffree-line-length-none','-I',str(source),str(f),'-o',str(exe)],check=True)
    return exe

def run(command,states,cwd):
    text=''.join(' '.join(format(v,'.17g') for v in s)+'\n' for s in states)
    r=subprocess.run(list(map(str,command)),input=text,text=True,capture_output=True,cwd=cwd,check=True)
    return [list(map(float,line.split())) for line in r.stdout.splitlines()]

def compare(x,y,names,states,label):
    assert len(x)==len(y)==len(states)
    worst=(0,None)
    failures=[]
    matched_nan=0
    for k,(r,s) in enumerate(zip(x,y)):
        assert len(r)==len(s)==len(names)
        for name,u,v in zip(names,r,s):
            if math.isnan(u) and math.isnan(v):
                matched_nan+=1
                continue # inherited invalid extreme-state diagnostics
            if u==v: continue
            # Near-zero derivatives suffer cancellation. The following absolute
            # floor is only appropriate for dimensionless Maxwell diagnostics.
            atol=3e-12 if name in ('dse','dpe','dsp') else 0
            scale=max(abs(u),abs(v),1e-290)
            if label.startswith('Direct'):
                # Dimensional conditioning scale: higher derivatives and pair-
                # dominated charge differences can be zero to roundoff. Compare
                # their absolute error to Q/(rho^i T^j), not to a cancelled Q'.
                d,t,_,_=states[k]
                values=dict(zip(names,r))
                # Entropy itself subtracts large degenerate energy/chemical-
                # potential terms; e/T is its cancellation scale.
                values['s']=max(abs(values['s']),abs(values['e'])/t)
                if name=='s': scale=max(scale,values['s'])
                derivatives={'pd':('p',1,0),'pt':('p',0,1),'pdd':('p',2,0),
                    'pdt':('p',1,1),'ptt':('p',0,2),'st':('s',0,1),'sdd':('s',2,0),
                    'etad':('eta',1,0),'etat':('eta',0,1),'etadt':('eta',1,1),
                    'nd':('n',1,0),'nt':('n',0,1),'ndt':('n',1,1)}
                if name in derivatives:
                    base,nd,nt=derivatives[name]
                    scale=max(scale,abs(values[base])/d**nd/t**nt)
            error=abs(u-v)/scale
            if error>worst[0] and name not in ('dse','dpe','dsp'): worst=(error,(states[k],name,u,v))
            if not(math.isfinite(u) and math.isfinite(v)) or abs(u-v)>atol+3e-10*scale:
                failures.append((states[k],name,u,v,error))
    print(label,'states=',len(states),'fields=',len(names),'worst=',worst,'failures=',len(failures),'matched nonfinite diagnostics=',matched_nan)
    for f in failures[:12]: print('  ',f)
    return len(failures)

rng=random.Random(61423)
states=[(2e9, 1e4, 12, 6)]
for i in range(400):
    abar,zbar=rng.choice([(1,1),(4,2),(12,6),(16,8),(56,26),(13.714285714,6.857142857)])
    states.append((10**rng.uniform(-11.9,14.9)*abar/zbar,10**rng.uniform(3.01,12.99),abar,zbar))
for d in (1e-12,1e-4,1,1e6,1e15):
    for t in (1e3,1e7,1e10,1e13): states.append((2*d,t,12,6))
with tempfile.TemporaryDirectory(prefix='helmholtz-reference-') as tmp:
    work=Path(tmp)
    (work/'helm_table.dat').symlink_to(a.helmholtz_source.resolve()/'helm_table.dat')
    ref=compile_reference(work,a.helmholtz_source.resolve(),'helmeos',[n+'_row(1)' for n in fields])
    x=run([ref],states,work); y=run([a.probe.resolve(),work/'helm_table.dat'],states,work)
    failures=compare(x,y,fields,states,'Tabulated Fortran/C++')
    names='p e s pd pt pdd pdt ptt st sdd eta etad etat etadt n nd nt ndt'.split()
    outputs=['pele_row(1)+ppos_row(1)','eele_row(1)+epos_row(1)','sele_row(1)+spos_row(1)',
        'dpepd_row(1)','dpept_row(1)','dpepdd_row(1)','dpepdt_row(1)','dpeptt_row(1)',
        'dsept_row(1)','dsepdd_row(1)','etaele_row(1)','detad_row(1)','detat_row(1)','detadt_row(1)',
        'xne_row(1)+xnp_row(1)','dxned_row(1)','dxnet_row(1)','dxnedt_row(1)']
    ref=compile_reference(work,a.timmes_source.resolve(),'eosfxt',outputs)
    # Focus on classical, degenerate, relativistic, and pair-producing states.
    direct_states=[(d,t,1,1) for d in (1e-6,1,1e3,1e6,1e9) for t in (1e4,1e6,1e8,1e9,1e10)]
    x=run([ref],direct_states,work); y=run([a.probe.resolve(),'--direct'],direct_states,work)
    failures+=compare(x,y,names,direct_states,'Direct Fortran/C++')
    raise SystemExit(1 if failures else 0)
