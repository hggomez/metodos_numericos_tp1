def solve_system(matrix: [[float]], independent_terms: [float]):
  SystemOfEquations(matrix, 

class SystemOfEquations:
  def __init__(self, matrix: [[float]], independent_terms: [float]):
    self.augmented_matrix = matrix
    for idx in range(len(matrix)):
      self.augmented_matrix[idx].append(independent_terms[idx])

  def plain_gaussian_elimination(self):
      
