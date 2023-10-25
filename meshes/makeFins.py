from tools.makeSimpleMeshes import makeSimpleMeshes as Mesher
import os

cwd = os.path.dirname(__file__)
body_file = os.path.join(cwd, "Body.vtu")

mesher = Mesher()

body = mesher.loadVTU(body_file)
canard = mesher.rectangularPrism(5, 2, 0.5, 10)
canard = mesher.transform(canard, 10, 0, 0, 0, 0, -45)
whole = mesher.appendObjects(body, canard)
mesher.render(whole)
#mesher.exportToVtu(file, body)