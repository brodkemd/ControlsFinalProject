import numpy as np

def Matrix(arg):
    m = np.array(arg, dtype=float)
    if m.shape[1] == 1:
        return np.atleast_2d(m).T[0]
    return m

def LoadModule(file:str, obj) -> None:
        options = {
            "sqrt"   : np.sqrt,
            "Matrix" : Matrix,
            "self"   : obj
        }
        with open(file, 'r') as f:
            for line in f.readlines():
                exec(f"self.{line}", options , vars(obj))
