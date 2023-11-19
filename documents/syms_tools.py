from sympy import cos, sin, Matrix, latex, zeros
import os

cwd = os.path.join(os.path.dirname(__file__), "sympy_output")
def write(s, _file:str):
    global cwd
    if not isinstance(s, str): s = latex(s)

    _file = os.path.join(cwd, _file)
    with open(_file, "w") as file: file.write(s)

def Euler2Quaternion(phi, theta, psi):
    """
    Converts an euler angle attitude to a quaternian attitude
    :param euler: Euler angle attitude in a np.matrix(phi, theta, psi)
    :return: Quaternian attitude in np.array(e0, e1, e2, e3)
    """

    e0 = cos(psi/2.0) * cos(theta/2.0) * cos(phi/2.0) + sin(psi/2.0) * sin(theta/2.0) * sin(phi/2.0)
    e1 = cos(psi/2.0) * cos(theta/2.0) * sin(phi/2.0) - sin(psi/2.0) * sin(theta/2.0) * cos(phi/2.0)
    e2 = cos(psi/2.0) * sin(theta/2.0) * cos(phi/2.0) + sin(psi/2.0) * cos(theta/2.0) * sin(phi/2.0)
    e3 = sin(psi/2.0) * cos(theta/2.0) * cos(phi/2.0) - cos(psi/2.0) * sin(theta/2.0) * sin(phi/2.0)

    return Matrix([e0, e1, e2, e3])


def writeMathModule(name:str, destination="controller", **kwargs):
    file = os.path.join(os.path.dirname(os.path.dirname(__file__)), destination, f"{name}.h")
    print("writing math module in:", file)

    with open(file, "w") as f:
        for arg in kwargs:
            f.write(f"{arg} = {kwargs[arg]}\n")




"""
DEPRICATED:

def Diff(f):
    return diff(f, t)



"""