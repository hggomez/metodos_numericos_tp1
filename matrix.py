class SystemOfEquations:
  def __init__(self, matrix: [[float]], independent_terms: [float]):
    self.augmented_matrix = matrix
    for idx in range(len(matrix)):
      self.augmented_matrix[idx].append(independent_terms[idx])
    self.height = len(self.augmented_matrix)
    self.width = len(self.augmented_matrix[0])

  def plain_gaussian_elimination(self):
    for col_idx in range(self.height):
      pivot = self.augmented_matrix[col_idx][col_idx]
      for row_idx in range(col_idx+1, self.height):
        elem_to_zero = self.augmented_matrix[row_idx][col_idx]
        k = elem_to_zero / pivot
        for elem_idx in range(self.width):
          self.augmented_matrix[row_idx][elem_idx] -= k * self.augmented_matrix[col_idx][elem_idx]
  
  def backward_substitution(self):
    solution = [0] * self.height
    for row_idx in reversed(range(self.height)):
      # Evalúa las variables a la derecha de la diagonal y suma sus valores
      solved_terms = sum(map(lambda x: x[0] * x[1], zip(solution[row_idx+1:self.height], self.augmented_matrix[row_idx][row_idx+1:self.height])))
      # Pasa los términos sumados, despeja la variable de la diagonal y resuelve
      solution[row_idx] = (self.augmented_matrix[row_idx][self.width-1] - solved_terms) / self.augmented_matrix[row_idx][row_idx]
    return solution

  def partial_pivoting_gaussian_elimination(self):
    for col_idx in range(self.height):

      #looks for the max posible pivot value
      pivot_row = col_idx
      pivot = self.augmented_matrix[col_idx][col_idx]
      
      for potential_pivot_row_idx in range(col_idx, self.height):
        if abs(self.augmented_matrix[potential_pivot_row_idx][col_idx]) > abs(pivot):
          pivot = self.augmented_matrix[potential_pivot_row_idx][col_idx]
          pivot_row = potential_pivot_row_idx

      # Swapirili
      self.augmented_matrix[pivot_row], self.augmented_matrix[col_idx] = self.augmented_matrix[col_idx], self.augmented_matrix[pivot_row]

      for row_idx in range(col_idx+1, self.height):
        elem_to_zero = self.augmented_matrix[row_idx][col_idx]
        k = elem_to_zero / pivot
        for elem_idx in range(self.width):
          self.augmented_matrix[row_idx][elem_idx] -= k * self.augmented_matrix[col_idx][elem_idx]
  
  def __str__(self):
    return "\n".join([str(row) for row in self.augmented_matrix])

class tridiagonalSOE:
  def __init__(self, a: [float], b: [float], c: [float], d: [float]):
    self.a = a # we assume a[0] = 0
    self.b = b 
    self.c = c # we assume c[n-1] = 0
    self.d = d
    self.k = [1 * i for i in range(len(a))]
    self.n = len(a)

  def plain_gaussian_elimination(self):
    for i in range(1, self.n):
      k = self.a[i] / self.b[i-1]
      self.k[i] = k # save k for later
      self.a[i] = self.a[i] - (k * self.b[i-1])
      self.b[i] = self.b[i] - (k * self.c[i-1])
      self.d[i] = self.d[i] - (k * self.d[i-1])

  def backward_substitution(self):
    x = [0] * self.n
    x[self.n-1] = self.d[self.n-1] / self.b[self.n-1]
    for i in range(self.n-2, -1, -1):
      x[i] = (self.d[i] - (self.c[i] * x[i+1])) / self.b[i]
    return x
  
  def set_new_d(self, d: [float]):
    self.d = d
    for i in range(1, self.n):
      self.d[i] = self.d[i] - (self.k[i] * self.d[i-1])

  def __str__(self):
    return "a: " + str(self.a) + "\nb: " + str(self.b) + "\nc: " + str(self.c) + "\nd: " + str(self.d)

print("=== PLAIN GAUSSIAN ===")
soe = SystemOfEquations([[1,2,3],[4,2,2],[1,1,1]], [6,6,6])
soe.plain_gaussian_elimination()
print("SOE: ")
print(soe)
print("Solution: ")
print(soe.backward_substitution())

try:
  soe = SystemOfEquations([[0,2,3],[4,2,2],[1,1,1]], [6,6,6])
  soe.plain_gaussian_elimination()
  print("\nSOE with [0][0] zero: ")
  print(soe)
  print("Solution: ")
  print(soe.backward_substitution())
except ZeroDivisionError:
  print("ERROR: Division by zero")

soe = SystemOfEquations([[4,2,2],[0,2,3],[1,1,1]], [6,6,6])
soe.plain_gaussian_elimination()
print("\nSOE with [1][0] zero: ")
print(soe)
print("Solution: ")
print(soe.backward_substitution())

print("\n=== PARTIAL PIVOTING GAUSSIAN ===")
soe = SystemOfEquations([[1,2,3],[4,2,2],[1,1,1]], [6,6,6])
soe.partial_pivoting_gaussian_elimination()
print("SOE: ")
print(soe)
print("Solution: ")
print(soe.backward_substitution())

soe = SystemOfEquations([[0,2,3],[4,2,2],[1,1,1]], [6,6,6])
soe.partial_pivoting_gaussian_elimination()
print("\nSOE with [0][0] zero: ")
print(soe)
print("Solution: ")
print(soe.backward_substitution())

soe = SystemOfEquations([[4,2,2],[0,2,3],[1,1,1]], [6,6,6])
soe.partial_pivoting_gaussian_elimination()
print("\nSOE with [1][0] zero: ")
print(soe)
print("Solution: ")
print(soe.backward_substitution())

print("\n=== TRIDIAGONAL SOE ===")
print(f"A:\n[1,2,0]\n[3,4,5]\n[0,6,7]")
print("\nSolutions for [6, 6, 6]: ")
soe = SystemOfEquations([[1,2,0],[3,4,5],[0,6,7]], [6,6,6])
soe.plain_gaussian_elimination()
print("SystemOfEquations class: ", end="")
print(soe.backward_substitution())
print("Tridiagonal class: ", end="")
triag_soe = tridiagonalSOE([0,3,6], [1,4,7], [2,5,0], [6,6,6])
triag_soe.plain_gaussian_elimination()
print(triag_soe.backward_substitution())

print("\nSolutions for [10, 10, 10]: ")
soe_2 = SystemOfEquations([[1,2,0],[3,4,5],[0,6,7]], [10,10,10])
soe_2.plain_gaussian_elimination()
print("SystemOfEquations class: ", end="")
print(soe_2.backward_substitution())
print("Tridiagonal class, using set_new_d: ", end="")
triag_soe.set_new_d([10,10,10])
print(triag_soe.backward_substitution())
print()
