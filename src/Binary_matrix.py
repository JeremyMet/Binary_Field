
# Binary Matrix Arithmetic
# Author : Métairie Jérémy
# Date : 20/02/2018

import copy

class Binary_matrix(object):


    def __init__(self, a):
        self.nb_rows = len(a) ;
        self.nb_cols = len(a[0]) ;
        self.shape = (self.nb_rows, self.nb_cols) ;
        self.mat = list(a) ;

    def __add__(self, other):
        if self.nb_rows != other.nb_rows or self.nb_cols != other.nb_cols:
            raise("Matrice dimensions do not match.") ;
        ret = [[0 for _ in range(self.nb_cols)] for _ in range(self.nb_rows)] ;
        for i in range(self.nb_rows):
            for j in range(self.nb_cols):
                ret[i][j] = self.mat[i][j]^other.mat[i][j] ;
        return Binary_matrix(ret) ;

    def __sub__(self, other):
        return self+other ;

    def __mul__(self, other):
        if self.nb_cols != other.nb_rows:
            raise ("Matrice Dimensions do not match.");
        ###
        ret = [[0 for _ in range(other.nb_cols)] for _ in range(self.nb_rows)];
        for i in range(other.nb_cols):
            for j in range(self.nb_rows):
                tmp = 0;
                for k in range(self.nb_cols):
                    tmp = tmp ^ (self.mat[j][k] & other.mat[k][i]);
                ret[j][i] = tmp ;
        return Binary_matrix(ret) ;

    def __pow__(self, other):
        if self.nb_cols != self.nb_rows:
            raise("Matrix should be square.") ;
        ret = Binary_matrix.eye(self.nb_cols) ;
        buf = copy.deepcopy(self) ;
        if other < 0:
            buf = ~buf ;
            other=-other ;
        while(other > 0):
            if other&1:
                ret = ret*buf ;
            buf = buf*buf ;
            other = other >> 1 ;
        return ret ;

    def __invert__(self):
        if not(Binary_matrix.determinant(self)):
            raise("Matrix is not invertible.") ;
        return Binary_matrix.gauss_elimination(self, "inv") ;

    def __eq__(self, other):
        if self.nb_cols != other.nb_cols or self.nb_rows != other.nb_rows:
            False ;
        for i in range(self.nb_rows):
            for j in range(self.nb_cols):
                if self.mat[i][j] != other.mat[i][j]:
                    return False ;
        return True ;

    def __ne__(self, other):
        return not(self==other) ;


    def __str__(self):
        return str(self.mat) ;

    def transpose(self):
        ret = [[0 for _ in range(self.nb_rows)] for _ in range(self.nb_cols)];
        for i in range(self.nb_rows):
            for j in range(self.nb_cols):
                ret[j][i] = self.mat[i][j];
        return Binary_matrix(ret);

    @classmethod
    def gauss_elimination(cls, M, inv = ''):
        if inv:
            if inv=="inv" and M.nb_cols != M.nb_rows:
                raise("Matrix should be square.")
            if inv=="ker":
                M = M.transpose() ;
            id = Binary_matrix.eye(M.nb_rows);
        # Run the Gaussian Elimination Algorithm
        M = copy.deepcopy(M);
        internal_cpt = 0;
        for i in range(M.nb_cols):
            found_pivot = False ;
            pivot_value = -1 ;
            # Look for a pivot
            for j in range(internal_cpt, M.nb_rows):
                if M.mat[j][i]:
                    found_pivot = True ;
                    pivot_value = j ;
                    break ;
            # Once pivot has been found, exchange rows
            if found_pivot:
                tmp = M.mat[pivot_value] ;
                M.mat[pivot_value] = M.mat[internal_cpt] ;
                M.mat[internal_cpt] = tmp ;
                if inv:
                    tmp = id.mat[pivot_value];
                    id.mat[pivot_value] = id.mat[internal_cpt];
                    id.mat[internal_cpt] = tmp;
            # And add this row to others if required
                for j in range(internal_cpt+1, M.nb_rows):
                    if M.mat[j][i]:
                        M.mat[j] = Binary_matrix.add_vector(M.mat[j], M.mat[internal_cpt]) ;
                        if inv:
                            id.mat[j] = Binary_matrix.add_vector(id.mat[j], id.mat[internal_cpt]);
                internal_cpt += 1;

        ## Solving the system (inversion).
        if inv=="inv":
            upper = M.nb_rows ;
            for i in reversed(range(upper)):
                if M.mat[i][i]:
                    for j in reversed(range(i)):
                        if (M.mat[j][i]):
                            M.mat[j] = Binary_matrix.add_vector(M.mat[j], M.mat[i]) ;
                            id.mat[j] = Binary_matrix.add_vector(id.mat[j], id.mat[i]);

        if inv:
            tmp = Binary_matrix.clean_matrix(M, id) ;
            return tmp[1] ;
        else:
            return M ;



    @classmethod
    def clean_matrix(cls, a, b=None):
        # "Clean the matrix" (i.e. move the 0 rows at its bottom)
        # If b matrix is provided, the permutation applied to a is also applied to b.
        internal_cpt = 0;
        z_cpt = a.nb_rows - 1;
        z = Binary_matrix.zero(1, a.nb_cols);
        a = copy.deepcopy(a) ;
        if b:
            b = copy.deepcopy(b) ;
        while (internal_cpt < z_cpt):
            if Binary_matrix([a.mat[internal_cpt]]) == z:
                tmp = a.mat[internal_cpt];
                a.mat[internal_cpt] = a.mat[z_cpt];
                a.mat[z_cpt] = tmp;
                z_cpt -= 1;
                if b:
                    tmp = b.mat[internal_cpt];
                    b.mat[internal_cpt] = b.mat[z_cpt];
                    b.mat[z_cpt] = tmp;
            if Binary_matrix([a.mat[internal_cpt]]) != z:
                internal_cpt += 1
        if b:
            return (a,b)
        else:
            return a ;

    @classmethod
    def add_vector(cls, a, b):
        if len(a) != len(b):
            raise("Vector dimensions do not match.")
        ret = [] ;
        for i in range(len(a)):
            ret.append(a[i]^b[i]) ;
        return ret ;

    @classmethod
    def determinant(cls, a):
        if a.nb_rows != a.nb_cols:
            raise("Matrix should be square.")
        a = Binary_matrix.gauss_elimination(a) ;
        ret = 1 ;
        for i in range(a.nb_rows):
            ret = ret & a.mat[i][i] ;
        return ret ;


    @classmethod
    def rank(cls, a):
        a = Binary_matrix.gauss_elimination(a);
        a = Binary_matrix.clean_matrix(a) ;
        z = Binary_matrix.zero(1, a.nb_cols) ;
        i = 0
        while(i < a.nb_rows and Binary_matrix([a.mat[i]])!= z):
            i+= 1 ;
        return i ;

    @classmethod
    def ker(cls, a):
        rk = Binary_matrix.rank(a);
        a = Binary_matrix.gauss_elimination(a, "ker") ;
        ker_dim = a.nb_rows-rk ;
        ker_vec = [] ;
        for i in range(a.nb_rows-ker_dim, a.nb_rows):
            ker_vec.append(a.mat[i]) ;
        if not(ker_vec):
            ker_vec = Binary_matrix.zero(a.nb_rows, a.nb_cols) ;
        else:
            ker_vec = Binary_matrix(ker_vec) ;
        return ker_vec ;

    @classmethod
    def eye(cls, n):
        delta = lambda i,j : 1 if (i==j) else 0 ;
        ret = [[delta(i,j) for i in range(n)] for j in range(n)] ;
        return Binary_matrix(ret) ;

    @classmethod
    def zero(cls, n, m):
        ret = [[0 for i in range(m)] for j in range(n)];
        return Binary_matrix(ret);

    @classmethod
    def solving_system(cls, a, v):
        a = Binary_matrix.gauss_elimination(a, True) ;
        return a*v ;





if __name__ == "__main__":
    a = [[1,0],[1,1],[0,0]] ;
    b = [[1,0],[1,1], [1,1]] ;
    m_a = Binary_matrix(a) ;
    m_b = Binary_matrix(b);
    print(m_a+m_b) ;
    print(m_a-m_b) ;
    print(m_a*m_a.transpose()+m_b*m_b.transpose()) ;
    print(m_a.transpose()) ;
    c = [[1,1,1],[1,0,1],[0,1,1]] ;
    m_c = Binary_matrix(c) ;
    print("Binary Inverse")
    print("Should be equal to ID") ;
    print(m_c*~m_c) ;
    v = Binary_matrix([[1,1,1]]) ;
    solution = Binary_matrix.solving_system(m_c, v.transpose()) ;
    print(m_c*solution) ;
    print(Binary_matrix.determinant(m_c)) ;
    print(Binary_matrix.eye(5)) ;
    print(m_c) ;
    print(m_c*m_c*m_c*m_c) ;
    print(m_c**4)
    print(m_c == m_c)
    print(m_c == m_c*Binary_matrix.eye(3))
    d =[[1,1,0,0],[0,1,1,1],[0,1,0,1]]
    m_d = Binary_matrix(d) ;
    print(Binary_matrix.rank(m_a)) ;
    print('Ker de m_d')
    print(Binary_matrix.ker(m_d));
    e = [[0, 0, 0, 0, 0, 0, 1], [0, 1, 0, 0, 1, 1, 0], [0, 1, 1, 0, 0, 0, 1], [0, 0, 0, 1, 1, 0, 0], [0, 0, 1, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 1, 1, 0, 0]] ;
    m_e = Binary_matrix(e) ;
    print("Value of m_e  "+str(m_e)) ;
    print(Binary_matrix.ker(m_e));
    print("Simple Example")
    M = Binary_matrix([[1, 0, 1], [0, 1, 1]]);
    print(M) ;
    M2 = M*M.transpose() ;
    print(M2) ;
    M2_inv = M2**-1 ;
    print(M2*M2_inv) ;





