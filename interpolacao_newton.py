from sympy import symbols, expand
def interpolacao_newton(x, y, valor):
    n = len(x)
    f = [[0] * n for _ in range(n)]

    for i in range(n):
        f[i][0] = y[i]

    for j in range(1, n):
        for i in range(n - j):
            f[i][j] = (f[i + 1][j - 1] - f[i][j - 1]) / (x[i + j] - x[i])

    resultado = f[0][0]
    termo = 1

    for i in range(1, n):
        termo *= (valor - x[i - 1])
        resultado += f[0][i] * termo

    return resultado



def polinomio_interpolador_newton(x, y):
    n = len(x)
    f = [[0] * n for _ in range(n)]

    for i in range(n):
        f[i][0] = y[i]

    for j in range(1, n):
        for i in range(n - j):
            f[i][j] = (f[i + 1][j - 1] - f[i][j - 1]) / (x[i + j] - x[i])

    x_symbol = symbols('x')
    polinomio = f[0][0]
    termo = 1

    for i in range(1, n):
        termo *= (x_symbol - x[i - 1])
        polinomio += f[0][i] * termo

    return expand(polinomio)




# Exemplo de uso
x = [0, 10, 20]
y = [0, 6, 4]
valor_interpolar = 5

resultado = interpolacao_newton(x, y, valor_interpolar)
print(f'O valor interpolado em {valor_interpolar} é aproximadamente {resultado}')

polinomio = polinomio_interpolador_newton(x, y)
print(f"O polinômio interpolador de Newton é: {polinomio}")
