from sympy import cos, sin, Matrix, latex, zeros
import os

cwd = os.path.join(os.path.dirname(__file__), "sympy_output")
to_replace = {
    "starboard" : "s",
    "port" : "p",
    "canard" : "c",
    "fin" : "f"
}
def write(s, _file:str):
    global cwd
    if not isinstance(s, str):
        s = latex(s)
        for item in to_replace:
            s = s.replace(item, to_replace[item])


    _file = os.path.join(cwd, _file)
    with open(_file, "w") as file: file.write(s)

def Euler2Quaternion(phi, theta, psi):
    """
    Converts an euler angle attitude to a quaternian attitude
    :param euler: Euler angle attitude in a np.matrix(phi, theta, psi)
    :return: Quaternian attitude in np.array(e0, e1, e2, e3)
    """

    e0 = cos(psi/2) * cos(theta/2) * cos(phi/2) + sin(psi/2) * sin(theta/2) * sin(phi/2)
    e1 = cos(psi/2) * cos(theta/2) * sin(phi/2) - sin(psi/2) * sin(theta/2) * cos(phi/2)
    e2 = cos(psi/2) * sin(theta/2) * cos(phi/2) + sin(psi/2) * cos(theta/2) * sin(phi/2)
    e3 = sin(psi/2) * cos(theta/2) * cos(phi/2) - cos(psi/2) * sin(theta/2) * sin(phi/2)

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