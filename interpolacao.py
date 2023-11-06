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
    '''

    função quu recebe o menor valor de x ou o maior valor de x e retorna o valor de f(x) derivado a self.ordem + 1 vezes
    '''
    def __f_derivado(self, x: int | float):
        return math.exp(x) * (x**2) + 12*math.exp(x) * x + 30*math.exp(x)  

    
    '''
    função que realiza o phi da formula de erro
    o phi é justamente voce pega o valor_interpolar e diminui com todos os x e no final multiplica-los
    '''
    def __phi(self, x_est: int | float):
        phi = 1
        for i in range(0, len(self.x)):
            phi *= (x_est - self.x[i])
        return phi
    
    '''
    função que realiza o erro da derivada ou o limitante superior
    '''
    def erro_com_derivado(self):
        num = self.__phi(self.valor_interpolar) / (self.__fat(self.ordem + 1))
        intervalo_0 = self.__f_derivado(self.x[0])
        intervalo_1 = self.__f_derivado(self.x[len(self.x) - 1])
        lista = [intervalo_0, intervalo_1]
        maximo = max(lista)
        maximo = abs(maximo)
        return num * maximo


    '''
    função que pega o f(x) da questao e retorna seu y de acordo com o valor_interpolar
    '''
    def __f(self, x: int | float):
        return x
        

    '''
        função que pega o polinomio descoberto e pega o valor de y de acordo com o valor_interpolar
    '''
    def __polynomial(self, x: float | int):
        x_var = sp.symbols('x')
        pol = sp.lambdify((x_var), self.polinomio, 'numpy')
        resultado = pol(x)
        return resultado

    '''
        função que faz o erro sem derivada
    '''
    def erro_sem_derivada(self, x: int | float):
        return abs(self.__f(x) - self.__polynomial(x))   

    '''
    função que faz um fatorial
    '''
    def __fat(self, n: int | float) -> int | float:
        resultado=1
        for num in range(1,n+1):
            resultado *= num
        return resultado
    
    '''função que printa a tabela delta'''
    def printar_tabela_delta(self):
        n_print = len(self.x)
        for i in range(n_print):
            print(self.x[i], end = "\t")
            for j in range(n_print - i):
                print(self.table_dely[i][j], end = "\t")
            print("")

    '''
    função que utiliza o gregory_newton
    '''
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

    '''
    função que faz o polinomio de gregory_newton
    '''
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
    
    '''
    função que utiliza a tecnica de newton
    '''
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

    '''função que faz o polinomio de newton'''
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
    
    '''função que utiliza a tecnica de lagrange'''
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

    '''
    função que faz o polinomio de lagrange
    '''
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
        self.polinomio = sp.N(polynomial, self.num_decimal_places)
        return sp.N(polynomial, self.num_decimal_places)


if __name__ == '__main__':
        # Exemplo de uso
    x = [5, 10, 15]  # lista de x  
    y = [0.9998, 0.9997, 0.9991] # lista de y
    valor_interpolar = 13 #valor que sera interpolado
    interpolacao = Interpolacao(x, y, valor_interpolar)
    print(interpolacao.lagrange_interpolation_polynomial())
    print(f'gregory: {interpolacao.gregory_newton_interpolation()}')
    print(f'lagrange: {interpolacao.lagrange_interpolation()}')
    print(f'newton: {interpolacao.interpolacao_newton()}')









    