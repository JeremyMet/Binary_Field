from GF2 import GF2 ;
from Binary_matrix import Binary_matrix ;


# Berlekamp Factorization Algorithm for Binary Polynomials.
# Author : Métairie Jérémy
# Date : 20/02/2018

class Berlekamp():

    # convert a binary vector [a_0, a_1 ... a_{n-1}] to the integer a_0+a_1*2+...+a_{n-1}*2^{n-1}
    @classmethod
    def from_vec_to_pol(cls, vec):
        ret = 0 ;
        for i in range(len(vec)):
            if (vec[i]%2 != 0):
                ret+=2**i ;
        return ret ;

    @classmethod
    def square_factors(cls, F):
        dF = F.derivate() ;
        return GF2(GF2.gcd(F.val, dF.val)) ;


    @classmethod
    ## Get root if F is written as G^(2n) (i.e. derivate(F) = 0).
    def get_root(cls, F, deg = None):
        F = F.val ;
        ret = 0 ;
        cpt = 0 ;
        while(F>0):
            if (F&1):
                ret = ret^(1<<cpt) ;
            cpt+=1 ;
            F = F >> 2 ;
        ret = GF2(ret) ;
        if not(deg):
            deg = 2 ;
        else:
            deg = deg << 1 ;
        if ret.derivate() == GF2(0):
            return Berlekamp.get_root(ret, deg) ;
        else:
            return (ret, deg) ;


    ## The factorization algorithm is based on the Berlekamp algorithm (Elwyn Berlekamp in 1967)
    @classmethod
    def Berlekamp_Alg(cls, F):
        F = F.val ;
        GF2.set_polynomial(F) ;
        deg = GF2.F_DEG ;
        mat = [] ;
        for i in range(deg):
            tmp = GF2.reduce(1 << (2*i)) ;
            vec = [0]*deg ;
            for j in range(deg):
                if (tmp&1):
                    vec[j] = 1 ;
                tmp = tmp >> 1 ;
            mat.append(vec) ;
        mat = Binary_matrix(mat) ;
        mat = mat.transpose() + Binary_matrix.eye(deg) ;
        mat = Binary_matrix.ker(mat) ;

        for i in range(mat.nb_rows):
            v_0 = mat.mat[i] ;
            v_1 = list(v_0) ;
            v_1[0] = v_1[0]^1 ;
            tmp_0 = GF2.gcd(F, Berlekamp.from_vec_to_pol(v_0));
            tmp_1 = GF2.gcd(F, Berlekamp.from_vec_to_pol(v_1));
            if tmp_0 != 1 and tmp_0 != F: return (GF2(F)/GF2(tmp_0), GF2(tmp_0)) ;
            if tmp_1 != 1 and tmp_1 != F: return (GF2(F) / GF2(tmp_1), GF2(tmp_1));
        return False ;

    @classmethod
    def Factorize(cls, F):
        ret = {} ;
        remain = True ;
        ret[F.val] = 1 ;
        while(remain):
            keys = list(ret.keys()) ;
            remain = False ;
            for k in keys:
                GF2k = GF2(k) ;
                # Squares management
                if GF2k.derivate() == GF2(0):
                    squares = Berlekamp.get_root(GF2k) ;
                    t0 = squares[0] ; t1 = squares[1] ;
                    if t0.val in ret:
                        ret[t0.val] += t1*ret[k] ;
                    else:
                        ret[t0.val] = t1*ret[k] ;
                    del ret[k] ;
                    remain = True ;
                # Squares management (with odd powers)
                else:
                    squares = Berlekamp.square_factors(GF2k) ;
                    if squares != GF2(1):
                        remain = True ;
                        free_squares = GF2k/squares ;
                        if free_squares.val in ret:
                            ret[free_squares.val] += ret[k] ;
                        else:
                            ret[free_squares.val] = ret[k];
                        if squares.val in ret:
                            ret[squares.val] += ret[k] ;
                        else:
                            ret[squares.val] = ret[k];
                        del ret[k] ;
                    else:
                        tmp = Berlekamp.Berlekamp_Alg(GF2k) ;
                        if tmp:
                            t0 = tmp[0] ;
                            t1 = tmp[1] ;
                            if t0.val in ret:
                                ret[t0.val] += ret[k] ;
                            else:
                                ret[t0.val] = ret[k] ;
                            if t1.val in ret:
                                ret[t1.val] += ret[k] ;
                            else:
                                ret[t1.val] = ret[k] ;
                        if tmp:
                            remain = True;
                            del ret[k]
        return ret ;


if __name__ == "__main__":
    print("Factorization with the Berlekamp Algorithm .")
    GF2.set_polynomial(0) ;
    a = GF2("x^3+x+1") ;
    b = GF2("x^4+x+1") ;
    c = GF2("x^5+x^2+1") ;
    d = a*a*a*b*b*c ;
    d = d**129 ;
    print("d: " + str(d));
    prod = Berlekamp.Factorize(d) ;
    print(" Factorisation ...") ;
    for k, v in prod.items():
        print("(key: "+str(GF2(k))+", pow: "+str(v)+")") ;


