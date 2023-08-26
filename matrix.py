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
  

print("=== PLAIN GAUSSIAN ===")
soe = SystemOfEquations([[1,2,3],[4,2,2],[1,1,1]], [6,6,6])
soe.plain_gaussian_elimination()
print("SOE: ")
print(soe)
print("SOLUTION: ")
print(soe.backward_substitution())

try:
  soe = SystemOfEquations([[0,2,3],[4,2,2],[1,1,1]], [6,6,6])
  soe.plain_gaussian_elimination()
  print("SOE: ")
  print(soe)
  print("SOLUTION: ")
  print(soe.backward_substitution())
except:
  pass

soe = SystemOfEquations([[4,2,2],[0,2,3],[1,1,1]], [6,6,6])
soe.plain_gaussian_elimination()
print("SOE WITH A ZERO: ")
print(soe)
print("SOLUTION: ")
print(soe.backward_substitution())

print("=== PARTIAL GAUSSIAN ===")
soe = SystemOfEquations([[1,2,3],[4,2,2],[1,1,1]], [6,6,6])
soe.partial_pivoting_gaussian_elimination()
print("SOE: ")
print(soe)
print("SOLUTION: ")
print(soe.backward_substitution())

soe = SystemOfEquations([[0,2,3],[4,2,2],[1,1,1]], [6,6,6])
soe.partial_pivoting_gaussian_elimination()
print("SOE WITH A ZERO AT FIRST: ")
print(soe)
print("SOLUTION: ")
print(soe.backward_substitution())

soe = SystemOfEquations([[4,2,2],[0,2,3],[1,1,1]], [6,6,6])
soe.partial_pivoting_gaussian_elimination()
print("SOE WITH A ZERO: ")
print(soe)
print("SOLUTION: ")
print(soe.backward_substitution())

