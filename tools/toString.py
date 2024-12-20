import numpy as np

def arrToStr(arr:np.ndarray, joins=[], starts=[], ends=[], num_format=":.16f", pad="", shift_for_minus=False):
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
            arr_str.append(arrToStr(arr[i], joins[1:], starts[1:], ends[1:], num_format=num_format, pad=pad, shift_for_minus=shift_for_minus))
        else:
            arr_str.append(numToStr(arr[i], format=num_format, pad=pad, shift_for_minus=shift_for_minus))
    return starts[0] + joins[0].join(arr_str) + ends[0]

def numToStr(num, format=":.16f", pad="", shift_for_minus=False):
    num_str = ""
    front_pad = ""
    if np.iscomplex(num):
        if shift_for_minus:
            if num.real > 0.0:
                front_pad = " "
        if num.imag < 0.0:
            num_str = (f"{front_pad}{{{format}}}{pad}-{pad}{{{format}}}j").format(num.real, abs(num.imag))
        else:
            num_str = (f"{front_pad}{{{format}}}{pad}+{pad}{{{format}}}j").format(num.real, num.imag)
    else:
        if shift_for_minus:
            if num > 0.0:
                front_pad = " "
        num_str = (f"{front_pad}{{{format}}}").format(num.real)
    return num_str