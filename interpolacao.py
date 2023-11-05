import sympy as sp
import math

class Interpolacao:
    def __init__(self, x: list[float | int], y: list[float | int], valor_interpolar: float | int, decimal_places=4) -> None:
        self.x = x
        self.y = y
        self.ordem = len(self.x) - 1
        self.valor_interpolar = valor_interpolar
        self.num_decimal_places = decimal_places
        self.table_dely = None
        self.polinomio = None

    def __f_derivado(self, x: int | float):
        return math.exp(x) * (x**2) + 12*math.exp(x) * x + 30*math.exp(x)  

    
    def __phi(self, x_est: int | float):
        phi = 1
        for i in range(0, len(self.x)):
            phi *= (x_est - self.x[i])
        return phi
    
    def erro_com_derivado(self):
        num = self.__phi(self.valor_interpolar) / (self.__fat(self.ordem + 1))
        intervalo_0 = self.__f_derivado(self.x[0])
        intervalo_1 = self.__f_derivado(self.x[len(self.x) - 1])
        lista = [intervalo_0, intervalo_1]
        maximo = max(lista)
        maximo = abs(maximo)
        return num * maximo


    def __f(self, x: int | float):
        return x
        

    
    def __polynomial(self, x: float | int):
        x_var = sp.symbols('x')
        pol = sp.lambdify((x_var), self.polinomio, 'numpy')
        resultado = pol(x)
        return resultado

    def erro_sem_derivada(self, x: int | float):
        return abs(self.__f(x) - self.__polynomial(x))   

    def __fat(self, n: int | float) -> int | float:
        resultado=1
        for num in range(1,n+1):
            resultado *= num
        return resultado
    
    def printar_tabela_delta(self):
        n_print = len(self.x)
        for i in range(n_print):
            print(self.x[i], end = "\t")
            for j in range(n_print - i):
                print(self.table_dely[i][j], end = "\t")
            print("")

    def gregory_newton_interpolation(self):
        n = len(self.x)
        if len(self.y) != n:
            raise ValueError("As listas x e y devem ter o mesmo tamanho")
        
        # Calcule o passo h
        h = self.x[1] - self.x[0]
        
        y = [[0]*n for _ in range(n) ]
        
        for i in range(0, n):
            y[i][0] = self.y[i]
        
        for i in range(1, n):
            for j in range(n - i):
                y[j][i] = y[j + 1][i - 1] - y[j][i - 1]

        

        sum = y[0][0]
        for ordem in range(1, n):
            delta = y[0][ordem]
            fatorial = self.__fat(ordem)
            valor = 1
            for i in range(0, ordem):
                valor *= (self.valor_interpolar - self.x[i])
            sum += (valor * delta) / (fatorial * pow(h, ordem))
        self.table_dely = y
        return round(sum, self.num_decimal_places)

    def gregory_newton_interpolation_equation(self):
        n = len(self.x)
        if len(self.y) != n:
            raise ValueError("As listas x e y devem ter o mesmo tamanho")
        xi = sp.symbols('x')
        # Calcule o passo h
        h = self.x[1] - self.x[0]
        
        y = [[0]*n for _ in range(n) ]
        
        for i in range(0, n):
            y[i][0] = self.y[i]
        
        for i in range(1, n):
            for j in range(n - i):
                y[j][i] = y[j + 1][i - 1] - y[j][i - 1]

        sum = y[0][0]
        for ordem in range(1, n):
            delta = y[0][ordem]
            fatorial = self.__fat(ordem)
            valor = 1
            for i in range(0, ordem):
                valor *= (xi - self.y[i])
            sum += (valor * delta) / (fatorial * pow(h, ordem))
        
        polinomio = sp.expand(sum)
        self.polinomio = polinomio
        return sp.N(polinomio, self.num_decimal_places)
    
    def interpolacao_newton(self):
        n = len(self.x)
        f = [[0] * n for _ in range(n)]

        for i in range(n):
            f[i][0] = self.y[i]

        for j in range(1, n):
            for i in range(n - j):
                f[i][j] = (f[i + 1][j - 1] - f[i][j - 1]) / (self.x[i + j] - self.x[i])

        resultado = f[0][0]
        termo = 1

        for i in range(1, n):
            termo *= (self.valor_interpolar - self.x[i - 1])
            resultado += f[0][i] * termo

        return round(resultado, self.num_decimal_places)

    def polinomio_interpolador_newton(self):
        n = len(self.x)
        f = [[0] * n for _ in range(n)]

        for i in range(n):
            f[i][0] = self.y[i]

        for j in range(1, n):
            for i in range(n - j):
                f[i][j] = (f[i + 1][j - 1] - f[i][j - 1]) / (self.x[i + j] - self.y[i])

        x_symbol = sp.symbols('x')
        polinomio = f[0][0]
        termo = 1

        for i in range(1, n):
            termo *= (x_symbol - self.x[i - 1])
            polinomio += f[0][i] * termo
        resultado = sp.expand(polinomio)
        self.polinomio = resultado
        return sp.N(resultado, self.num_decimal_places)
    
    def lagrange_interpolation(self):
        n = len(self.x)
        result = 0.0

        for i in range(n):
            term = self.y[i]
            for j in range(n):
                if i != j:
                    term *= (self.valor_interpolar - self.x[j]) / (self.x[i] - self.x[j])
            result += term

        return round(result, self.num_decimal_places)

    def lagrange_interpolation_polynomial(self):
        n = len(self.x)
        xi = sp.symbols('x')
        polynomial = 0

        for i in range(n):
            term = self.y[i]
            for j in range(n):
                if i != j:
                    term *= (xi - self.x[j]) / (self.x[i] - self.x[j])
            polynomial += term

        polynomial = sp.expand(polynomial)
        self.polinomio = polynomial
        return sp.N(polynomial, self.num_decimal_places)


if __name__ == '__main__':
        # Exemplo de uso
    x = [5, 10, 15]
    y = [0.9998, 0.9997, 0.9991]
    valor_interpolar = 13
    interpolacao = Interpolacao(x, y, valor_interpolar)
    print(f'gregory: {interpolacao.gregory_newton_interpolation()}')
    print(f'lagrange: {interpolacao.lagrange_interpolation()}')
    print(f'newton: {interpolacao.interpolacao_newton()}')









    