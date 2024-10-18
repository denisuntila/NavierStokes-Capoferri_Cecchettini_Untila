import math

class NacaAirfoil:

  def __init__(self, file_name):
    with open(file_name, 'r') as file:
      lines = [line.strip() for line in file]
      self.name = lines[0]
      lines.pop(0)
      self.data = []
      self.angle = 0

      for line in lines:
        temp = []
        strings = line.split()
        temp.append(float(strings[0]) - 0.5)
        temp.append(float(strings[1]))

        self.data.append(temp)
      self.chord = 1.0

  def resize(self, chord):
    temp = chord / self.chord
    for i in range(len(self.data)):
      for j in range(2):
        self.data[i][j] *= temp
    
    self.chord = chord
  
  def rotate(self, angle):
    self.angle -= angle
    cosalpha = math.cos(-angle)
    sinalpha = math.sin(-angle)

    for i in range(len(self.data)):
      point = self.data[i]
      self.data[i][0] = cosalpha * point[0] - sinalpha * point[1]
      self.data[i][1] = sinalpha * point[0] + cosalpha * point[1]
    
      

class Mesh:
  def __init__(self):
    self.lx = 2.2
    self.ly = 0.41
    self.dx = 0.4
    self.dy = 0.2
    self.h = 0.01
    self.points = []
    self.lines = []
    self.loops = []


  def add_airfoil(self, airfoil_data):
    for point in airfoil_data:
      self.points.append(point)
    
    for i in range(len(airfoil_data) - 1):
      self.lines.append([i, i + 1])
    self.lines.append([len(airfoil_data) - 1, 0])

    temp = []
    for i in range(len(airfoil_data)):
      temp.append(i)
    self.loops.append(temp)

    



  
  
  def write_file(self, file_name):
    buffer = "// Domain size\n"
    buffer += ("Lx = " + "{:10.4f}".format(self.lx) + ";\n")
    buffer += ("Ly = " + "{:10.4f}".format(self.ly) + ";\n")

    buffer += ("\n// Coordinates of the center of the airfoil\n")
    buffer += ("Dx = " + "{:10.4f}".format(self.dx) + ";\n")
    buffer += ("Dy = " + "{:10.4f}".format(self.dy) + ";\n")

    buffer += ("\nh = " + "{:10.4f}".format(self.h) + ";\n")

    buffer += ("\n// Points\n")
    buffer += ("Point(0) = {0, 0, 0, h};\n")
    buffer += ("Point(1) = {Lx, 0, 0, h};\n")
    buffer += ("Point(2) = {Lx, Ly, 0, h};\n")
    buffer += ("Point(3) = {0, Ly, 0, h};\n")

    buffer += "\n"
    for i, point in enumerate(self.points):
      sign = ["+", "+"]
      if point[0] < 0:
        sign[0] = ""
      if point[1] < 0:
        sign[1] = ""

      buffer += ("Point(" + str(i + 4) + ") = {Dx " + sign[0] + "{:.5f}".format(point[0])
      + ", Dy " + sign[1] + "{:.5f}".format(point[1]) + ", 0, h};\n")


    buffer += ("\n\n// Lines\n")
    buffer += ("Line(0) = {0, 1};\n")
    buffer += ("Line(1) = {1, 2};\n")
    buffer += ("Line(2) = {2, 3};\n")
    buffer += ("Line(3) = {3, 0};\n")

    buffer += "\n"
    for i, line in enumerate(self.lines):
      buffer += ("Line(" + str(i + 4) + ") = {" + str(line[0] + 4) + ", "
      + str(line[1] + 4) + "};\n")


    buffer += ("\n\n// Loops\n")
    buffer += ("Line Loop(1) = {0, 1, 2, 3};\n")
    buffer += "\n"

    for i, loop in enumerate(self.loops):
      buffer += ("Line Loop(" + str(i + 2) + ") = {") 
      for line_index in loop[:-1]:
        buffer += (str(line_index + 4) + ", ") 
      buffer += (str(loop[-1] + 4) + "};\n")
    
    buffer += "\n\n// Surfaces\n"
    buffer += ("Plane Surface(0) = {1, 2};\n")


    buffer += ("\n\n// Physical entities\n")
    buffer += ("Physical Line(0) = {0};\n")
    buffer += ("Physical Line(1) = {1};\n")
    buffer += ("Physical Line(2) = {2};\n")
    buffer += ("Physical Line(3) = {3};\n")

    for i, loop in enumerate(self.loops):
      buffer += ("Physical Line(" + str(i + 4) + ") = {") 
      for line_index in loop[:-1]:
        buffer += (str(line_index + 4) + ", ") 
      buffer += (str(loop[-1] + 4) + "};\n")

    buffer += ("Physical Surface(10) = {0};\n")


    buffer += "\n\nMesh 2;\n"


    with open(file_name, 'w') as file:
      file.write(buffer)




if __name__ == "__main__":
    a = NacaAirfoil("naca.dat")
    a.resize(0.5)

    angle = 13 #degrees

    file_name = a.name.replace(" ", "_") + ".geo"
    a.rotate(angle * math.pi / 180.0)
    mesh = Mesh()
    mesh.add_airfoil(a.data)
    mesh.write_file(file_name)
    

