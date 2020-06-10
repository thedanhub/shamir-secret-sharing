'''
An implementation of the Shamir secret sharing scheme
https://en.wikipedia.org/wiki/Shamir%27s_Secret_Sharing
'''

import secrets
import numpy as np


class ShSeSh:
    '''
        k = the key
        p = prime number for modulo operations in the set Z_p
        w = number of shares
        t = minimum number of shares required to reconstruct the key k
    '''

    def __init__(self, k, p, w, t):
        self.k = k
        self.p = p
        self.w = w
        self.t = t

    '''
        SHARES GENERATION
    '''

    '''
        choose_x chooses w distinct, non-zero elements of Z_p
        these are the public x_i values, 1 <= i <= w
    '''

    def choose_x(self):
        x = []
        for element in range(0, self.w):
            new_x = secrets.randbelow(self.p)
            while new_x == 0 or new_x in x:  # elements must be non-zero and distinct
                new_x = secrets.randbelow(self.p)
            x.append(new_x)
        self.x = x
        return self.x

    '''
        choose_a chooses t-1 elements of Z_p
        these are the random a_i values, 1 <= i <= t-1
    '''

    def choose_a(self):
        a = []
        for element in range(0, self.t - 1):
            a.append(secrets.randbelow(self.p))
        self.a = a
        return self.a

    '''
        generate_shares creates the secret shares to distribute
        these are the y_i values, 1 <= i <= w
    '''

    def generate_shares(self):
        y = []
        for idx, element in enumerate(range(0, self.w)):
            ax = np.poly1d(self.a + [self.k])  # ax is the polynomial we create before computing y_i
            y.append(ax(self.x[idx]) % self.p)  # evaluate the polynomial to generate y_i values
        self.y = y
        return self.y

    '''
        KEY RECONSTRUCTION
    '''

    '''
        The following two functions are used to calculate the inverse modulo p
        of a number (we need it when we calculate the b values, to handle
        the denominator)
        
        Adapted from https://en.wikibooks.org/wiki/Algorithm_Implementation/Mathematics/Extended_Euclidean_algorithm
    '''

    def egcd(self, number1, number2):
        if number1 == 0:
            return (number2, 0, 1)
        else:
            g, y, x = self.egcd(number2 % number1, number1)
            return g, x - (number2 // number1) * y, y

    def modinv(self, number, m):
        g, x, y = self.egcd(number, m)
        if g != 1:
            raise Exception('modular inverse does not exist')
        else:
            return x % m

    '''
        calculate_b is needed to reconstruct the key from the shares
        Part of the Lagrange interpolation formula
    '''

    def calculate_b(self, x, x_values):
        acc = 1
        for element in x_values:
            acc = acc * element
            if (
                    element - x != 0):  # ignore the case when we are using x, we need to focus on all of the other values only
                if (element - x >= 0):
                    inverse = self.modinv(element - x, self.p)
                else:
                    inverse = self.modinv(self.p - abs(element - x),
                                          self.p)  # https://math.stackexchange.com/questions/355066/find-the-inverse-modulo-of-a-number-got-a-negative-result
                acc = acc * inverse
        acc = acc / x  # we need to multiply all the values but the current one. Because we originally multiplied times x too, we divide by x
        return acc % self.p

    '''
        Take the values of x and the values of y as input
        (i.e. the coordinates of the points on the plane)
        and reconstruct the key by summing the products of b_i * y_i mod p
    '''

    def reconstruct_key(self, x, y):
        if (len(y) < self.t):
            print("Not enough shares. Key reconstruction won't be possible.")
        else:
            b = []
            k = 0
            for index in range(len(y)):
                b.append(self.calculate_b(x[index], x))
                k = k + b[index] * y[index]
            return k % self.p


#x = ShSeSh(13, 17, 5, 3)
#x = ShSeSh(1234, 1613, 6, 3)
x = ShSeSh(1234, 31847, 10, 5)
print("x = " + str(x.choose_x()))
print("a = " + str(x.choose_a()))
print("y = " + str(x.generate_shares()))
#print(x.reconstruct_key([1, 3, 5], [8, 10, 11]))
print(x.reconstruct_key([413, 432, 451, 470, 489], [25439, 14847, 24780, 5910, 12734]))



