import sympy as sp

def fat(n):
    resultado=1
    for num in range(1,n+1):
        resultado *= num
    return resultado


def printar_tabela_delta(x_print, y_ptint):
    n_print = len(x_print)
    for i in range(n_print):
        print(x_print[i], end = "\t")
        for j in range(n_print - i):
            print(y_ptint[i][j], end = "\t")
        print("")

def gregory_newton_interpolation(x1, y1, xi, num_decimal_places=2):
    n = len(x1)
    if len(y1) != n:
        raise ValueError("As listas x e y devem ter o mesmo tamanho")
    
    # Calcule o passo h
    h = x1[1] - x1[0]
    
    y = [[0]*n for _ in range(n) ]
    
    for i in range(0, n):
        y[i][0] = y1[i]
    
    for i in range(1, n):
        for j in range(n - i):
            y[j][i] = y[j + 1][i - 1] - y[j][i - 1]

    

    sum = y[0][0]
    for ordem in range(1, n):
        delta = y[0][ordem]
        fatorial = fat(ordem)
        valor = 1
        for i in range(0, ordem):
            valor *= (xi - x1[i])
        sum += (valor * delta) / (fatorial * pow(h, ordem))
    return round(sum, num_decimal_places), y

def gregory_newton_interpolation_equation(x1, y1, num_decimal_places=2):
    n = len(x1)
    if len(y1) != n:
        raise ValueError("As listas x e y devem ter o mesmo tamanho")
    xi = sp.symbols('x')
    # Calcule o passo h
    h = x1[1] - x1[0]
    
    y = [[0]*n for _ in range(n) ]
    
    for i in range(0, n):
        y[i][0] = y1[i]
    
    for i in range(1, n):
        for j in range(n - i):
            y[j][i] = y[j + 1][i - 1] - y[j][i - 1]

    sum = y[0][0]
    for ordem in range(1, n):
        delta = y[0][ordem]
        fatorial = fat(ordem)
        valor = 1
        for i in range(0, ordem):
            valor *= (xi - x1[i])
        sum += (valor * delta) / (fatorial * pow(h, ordem))
    
    polinomio = sp.expand(sum)
    return sp.N(polinomio, num_decimal_places)



if __name__ == '__main__':
    x = [0, 10, 20]
    y = [0, 6, 4]
    valor_interpolar = 5
    interpolated_value, table_dely = gregory_newton_interpolation(x, y, valor_interpolar)
    print(f"Para xi = {valor_interpolar}, o valor interpolado é {interpolated_value}")
    printar_tabela_delta(x, table_dely)
    polinomio = gregory_newton_interpolation_equation(x, y)
    print(f"A equação;")
    print(polinomio)