from sympy import *  

# Creates the symbols 
dPC = Symbol('dPC',real=true)
dSC = Symbol('dSC',real=true)
dPF = Symbol('dPF',real=true)
dSF = Symbol('dSF',real=true)
AOA = Symbol('AOA',real=true)
rho = Symbol('rho',real=true)
AsurfF = Symbol('AsurfF',real=true)
AsurfC = Symbol('AsurfC',real=true)
V = Symbol('V',real=true)
B = Symbol('B',real=true)
Momx = Symbol('Momx',real=true)
Momy = Symbol('Momy',real=true)
Momz = Symbol('Momz',real=true)
Fx = Symbol('Fx',real=true)
Fy = Symbol('Fy',real=true)
Fz = Symbol('Fz',real=true)
PP = Symbol('PP',real=true)
# Defines the surface coefficients from flat plate thy [delPortCanard, delStarCanard, delPortFin, delStarFin]

# Cls = 2*pi*delSurf
# Cds = 1.28*sin(delSurf)

# (1 - 0.5*(dPC - pi/2)**2 + (1/24)*(dPC - pi/2)**4 - (1/720)*(dPC- pi/2)**6)
# (1 - 0.5*(dSC - pi/2)**2 + (1/24)*(dSC - pi/2)**4 - (1/720)*(dSC- pi/2)**6)
# (1 - 0.5*(dPF - pi/2)**2 + (1/24)*(dPF - pi/2)**4 - (1/720)*(dPF- pi/2)**6)
# (1 - 0.5*(dSF - pi/2)**2 + (1/24)*(dSF - pi/2)**4 - (1/720)*(dSF- pi/2)**6)

# Canard forces (portCanard (pc), starboardCanard (sc))
Fpcx = (-1.28*cos(AOA)*cos(B) - 2*pi*dPC*sin(AOA)*sin(B))*(0.5*rho*(V**2)*AsurfC)
Fpcy = (1.28*cos(AOA)*sin(B) - 2*pi*dPC*sin(AOA)*cos(B))*(0.5*rho*(V**2)*AsurfC)
Fpcz = (-1.28*sin(AOA))*(0.5*rho*(V**2)*AsurfC)

Fscx = (-1.28*cos(AOA)*cos(B) + 2*pi*dSC**sin(AOA)*sin(B))*(0.5*rho*(V**2)*AsurfC)
Fscy = (1.28*cos(AOA)*sin(B) + 2*pi*dSC**sin(AOA)*cos(B))*(0.5*rho*(V**2)*AsurfC)
Fscz = (-1.28*sin(AOA))*(0.5*rho*(V**2)*AsurfC)

# Fin forces (portFin (pf), starboardFin (sf))
Fpfx = (-1.28*cos(AOA)*cos(B) - 2*pi*dPF*sin(AOA)*sin(B))*(0.5*rho*(V**2)*AsurfF)
Fpfy = (1.28*cos(AOA)*sin(B) - 2*pi*dPF*sin(AOA)*cos(B))*(0.5*rho*(V**2)*AsurfF)
Fpfz = (-1.28*sin(AOA))*(0.5*rho*(V**2)*AsurfF)

Fsfx = (-1.28*cos(AOA)*cos(B) + 2*pi*dSF*sin(AOA)*sin(B))*(0.5*rho*(V**2)*AsurfF)
Fsfy = (1.28*cos(AOA)*sin(B) + 2*pi*dSF*sin(AOA)*cos(B))*(0.5*rho*(V**2)*AsurfF)
Fsfz = (-1.28*pi*dSF*sin(AOA))*(0.5*rho*(V**2)*AsurfF)


Fsx = Fpcx + Fscx + Fpfx + Fsfx
Fsy = Fpcy + Fscy + Fpfy + Fsfy
Fsz = Fpcz + Fscz + Fpfz + Fsfz

rpc = Matrix(([15.25, 6.535, 0]))
rsc = Matrix(([15.25, -6.535, 0])) 
rpf = Matrix(([-19.375, 6.75, 0]))
rsf = Matrix(([-19.375, -6.75, 0]))


# Sets the force vectors
Fpc = Matrix(([Fpcx,Fpcy,Fpcz]))
Fsc = Matrix(([Fscx,Fscy,Fscz]))
Fpf = Matrix(([Fpfx,Fpfy,Fpfz]))
Fsf = Matrix(([Fsfx,Fsfy,Fsfz]))


# Calculates the radius vector of the center of pressure
PortCanard = Fpc.cross(rpc.cross(Fpc))
StarCanard = Fpc.cross(rsc.cross(Fsc))
PortFin = Fpf.cross(rpf.cross(Fpf))
StarFin = Fsf.cross(rsf.cross(Fsf))

A = Matrix(([Fsx,Fsy,Fsz]))
normsqr = (A.norm())**2 
COPrad = (PortCanard + StarCanard + PortFin + StarFin)/normsqr
COPx = COPrad[0]
CopNorm = COPrad.norm()
Mpc = rpc.cross(Fpc)
Msc = rsc.cross(Fsc)
Mpf = rpf.cross(Fpf)
Msf = rsf.cross(Fsf)
M = Mpc + Msc + Mpf + Msf

Mx = M[0]
My = M[1]
Mz = M[2]


Inverted = solve([Mx-Momx, Mz - Momz, My - Momy, Fsx-Fx], [dPC, dSC, dPF, dSF], dict=True) 
# dSCInverted = solve((Fsx-Fx, Fsy - Fy, Fsz - Fz, CopNorm-PP), dSC)
# dPFInverted = solve((Fsx-Fx, Fsy - Fy, Fsz - Fz, CopNorm-PP), dPF)
# dSFInverted = solve((Fsx-Fx, Fsy - Fy, Fsz - Fz, CopNorm-PP), dSF)


print(Inverted)


