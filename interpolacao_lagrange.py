import sympy as sp

def lagrange_interpolation(x, y, xi, num_decimal_places=2):
    n = len(x)
    result = 0.0

    for i in range(n):
        term = y[i]
        for j in range(n):
            if i != j:
                term *= (xi - x[j]) / (x[i] - x[j])
        result += term

    return round(result, num_decimal_places)

def lagrange_interpolation_polynomial(x, y, num_decimal_places=2):
    n = len(x)
    xi = sp.symbols('x')
    polynomial = 0

    for i in range(n):
        term = y[i]
        for j in range(n):
            if i != j:
                term *= (xi - x[j]) / (x[i] - x[j])
        polynomial += term

    polynomial = sp.expand(polynomial)
    return sp.N(polynomial, num_decimal_places)

# Letra A 2° grau
x2 = [20, 25, 30]
y2 = [0.9982, 0.9971, 0.9957]
valor = 27
decimal_places = 6 

result_2_grau = lagrange_interpolation(x2, y2, valor, decimal_places)
print(f'O valor interpolado {valor} em 2° grau é aproximadamente {result_2_grau}')

interpolation_polynomial = lagrange_interpolation_polynomial(x2, y2, decimal_places)
print("Polinômio interpolador de Lagrange de 2° grau:")
print(interpolation_polynomial)
print('========================')
x3 = [20, 25, 30, 35]
y3 = [0.9982, 0.9971, 0.9957, 0.9941]
decimal_places = 9
result_3_grau = lagrange_interpolation(x3, y3, valor, decimal_places)
print(f'O valor interpolado {valor} em 3° grau é aproximadamente {result_3_grau}')

interpolation_polynomial = lagrange_interpolation_polynomial(x3, y3, decimal_places)
print("Polinômio interpolador de Lagrange de 3° grau:")
print(interpolation_polynomial)