
def info(level, *args, **kwargs):
    space = "  "
    given = ""
    if not isinstance(level, int):
        given = level
        level = 0
    if not len(given):
        print(level*space, *args, **kwargs)
    else:
        print(level*space, given, *args, **kwargs)

def error(*args, **kwargs):
    print("Error:", *args, **kwargs)
    print("Exiting ...")
    exit()