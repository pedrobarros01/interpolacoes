import sympy as sp

def lagrange_interpolation(x, y, xi):
    n = len(x)
    result = 0.0

    for i in range(n):
        term = y[i]
        for j in range(n):
            if i != j:
                term *= (xi - x[j]) / (x[i] - x[j])
        result += term

    return result


def lagrange_interpolation_polynomial(x, y):
    n = len(x)
    xi = sp.symbols('x')
    polynomial = 0

    for i in range(n):
        term = y[i]
        for j in range(n):
            if i != j:
                term *= (xi - x[j]) / (x[i] - x[j])
        polynomial += term

    return sp.expand(polynomial)
# Exemplo de uso
x1 = [0, 10, 20]
y1 = [0, 6, 4]
valor1 = 5
result = lagrange_interpolation(x1, y1, valor1)
print(f'O valor interpolado em {valor1} é aproximadamente {result}')
interpolation_polynomial = lagrange_interpolation_polynomial(x1, y1)
print("Polinômio interpolador de Lagrange:")
print(interpolation_polynomial)
