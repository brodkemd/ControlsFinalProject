import numpy as np

def arrToStr(arr:np.ndarray, joins=[], starts=[], ends=[]):
    def checks(obj, default):
        if not isinstance(obj, list): obj = [obj]
        if not len(obj):              obj = [default]
        return obj

    joins   = checks(joins, ",")
    starts  = checks(starts, "[")
    ends    = checks(ends, "]")

    arr_str = []
    for i in range(len(arr)):
        if isinstance(arr[i], np.ndarray):
            arr_str.append(arrToStr(arr[i], joins[1:], starts[1:], ends[1:]))
        else:
            if np.iscomplex(arr[i]):
                if arr[i].imag < 0.0:
                    arr_str.append("{:.16f}-{:.16f}j".format(arr[i].real, abs(arr[i].imag)))
                else:
                    arr_str.append("{:.16f}+{:.16f}j".format(arr[i].real, arr[i].imag))
            else:
                arr_str.append("{:.16f}".format(arr[i]))
    return starts[0] + joins[0].join(arr_str) + ends[0]

def numToStr(num, format=":.16f"):
    num_str = ""
    if np.iscomplex(num):
        if num.imag < 0.0:
            num_str = (f"{{{format}}}-{{{format}}}j").format(num.real, abs(num.imag))
        else:
            num_str = (f"{{{format}}}+{{{format}}}j").format(num.real, num.imag)
    else:
        num_str = (f"{{{format}}}").format(num)