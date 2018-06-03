
# GF2 Class for Binary Field Manipulation
# Author : Métairie Jérémy
# Date : 20/02/2018

import copy ;
from ParsingTree import ParsingTree ;

class GF2(object):

    F = 19 ; # x^4+x+1
    F_DEG = 4 ;

    def __init__(self, a):
        if type(a) == str:
            self.val = GF2.build_from_string(a) ;
        else:
            self.val = a ;
        self.deg = GF2.degree(self.val) ;

    def __add__(self, other):
        return GF2(self.val ^ other.val) ;

    def __eq__(self, other):
        return self.val == other.val ;

    def __ne__(self, other):
        return not(self==other) ;

    def __sub__(self, other):
        return GF2(self.val ^ other) ;

    def __mul__(self, other):
        a = self.val ;
        b = other.val ;
        ret = 0 ;
        while(a>0):
            if (a&1):
                ret = ret^b ;
                tmp = ret ^ GF2.F ;
            a = a >> 1 ;
            b = b << 1 ;
            if (GF2.F != 0):
                tmp = b ^ GF2.F ;
                if tmp < b :
                    b = tmp ;
        return GF2(ret) ;

    def __mod__(self, other):
        a = self.val ;
        if type(other) == str:
            other = GF2.build_from_string(a) ;
        other = other.val ;
        if not(other):
            raise("Div by zero.") ;
        deg_other = GF2.degree(other) ;
        while(GF2.degree(a) >= deg_other):
            tmp = GF2.degree(a) - deg_other ;
            a = a ^ (other << tmp) ;
        return GF2(a) ;

    def __invert__(self):
        r0 = GF2.F; r1 = self.val ;
        if not(r1):
            raise("Div by zero.") ;
        u0 = 0 ; u1 = 1 ;
        while(r0 != 1):
            tmp = GF2.degree(r0) - GF2.degree(r1) ;
            r0 = r0 ^ (r1 << tmp) ;
            u0 = u0 ^ (u1 << tmp) ;
            if GF2.degree(r0) < GF2.degree(r1):
                r0, r1 = r1, r0 ;
                u0, u1 = u1, u0 ;
        return GF2(GF2.reduce(u0)) ;

    def __floordiv__(self, other):
        ret = 0 ;
        if not(other.val):
            raise("Div by zero.") ;
        a = self.val ; deg_a = GF2.degree(a) ;
        b = other.val ; deg_b = GF2.degree(b) ;
        while(deg_a >= deg_b):
            tmp = deg_a-deg_b ;
            ret = ret^(1<<tmp) ;
            a = a^(b<<tmp) ;
            deg_a = GF2.degree(a) ;
        ret = GF2(ret) ;
        return ret ;

    def __truediv__(self, other):
        print(GF2.gcd(other.val, GF2.F)) ;
        if GF2.gcd(other.val, GF2.F) != 1:
            raise("Element is not invertible.")
        return self*~other ;

    def __str__(self):
        ret = "" ;
        a = self.val ;
        cpt = 0 ;
        if (a==0):
            return "0" ;
        while(a>0):
            if (a&1):
                if ret: ret += "+";
                if (cpt == 0):
                    ret+="1" ;
                elif (cpt ==1 ):
                    ret+="x" ;
                else:
                    ret+="x^"+str(cpt) ;
            cpt+=1 ;
            a = a >> 1 ;
        return ret ;

    def __pow__(self, power):
        buffer = copy.copy(self) ;
        if power < 0:
            power = -power ;
            buffer = ~buffer ;
        ret = GF2(1) ;
        while(power>0):
            if power&1:
                ret = ret*buffer ;
            buffer = buffer*buffer ;
            power = power >> 1 ;
        return ret ;

    def derivate(self):
        a = self.val ;
        a = a >> 1 ;
        i = 1 ;
        while((1<< i)<=a):
            if (a & (1 << i)):
                a = a^(1<<i) ;
            i+=2 ;
        return GF2(a) ;

    @classmethod
    def build_from_string(cls, s):
        ret = 0 ;
        if "(" in s or "*" in s or "/" in s:
            tmp_tree = ParsingTree(s) ;
            ret = tmp_tree.eval(GF2).val ;
        else:
            s = s.split("+") ;
            for c in s:
                if '^' in c:
                    c = c.split("^") ;
                    ret+= 2**int(c[1]) ;
                else:
                    if c == '1':
                        ret += 1 ;
                    elif c == 'x':
                        ret+= 2 ;
        return ret ;

    @classmethod
    def degree(cls, a):
        dg = 0 ;
        a = a >> 1 ;
        while(a>0):
            dg+=1 ;
            a = a >> 1 ;
        return dg ;

    @classmethod
    def reduce(cls, a):
        if type(a) == str:
            a = GF2.build_from_string(a) ;
        while(GF2.degree(a) >= GF2.F_DEG):
            tmp = GF2.degree(a) - GF2.F_DEG ;
            a = a ^ (GF2.F << tmp) ;
        return a ;

    @classmethod
    def set_polynomial(cls, F):
        if type(F) == str:
            GF2.F = GF2.build_from_string(F) ;
        else:
            GF2.F = F ;
        GF2.F_DEG = GF2.degree(GF2.F) ;

    @classmethod
    def gcd(cls, a, b):
        if type(a) == str:
            a = GF2.build_from_string(a) ;
        if type(b) == str:
            b = GF2.build_from_string(b);
        if a == 0 or b == 0:
            return 0 ;
        r0 = a ; r1 = b ;
        while(r1!=0):
            tmp = GF2.degree(r0) - GF2.degree(r1);
            if tmp < 0:
                r0, r1 = r1, r0;
                tmp = -tmp ;
            r0 = r0 ^ (r1 << tmp) ;
            if r0 == 0:
                r0, r1 = r1, r0 ;
        return r0 ;

    @classmethod
    def derive(cls, a):
        pass



if __name__ == "__main__":
    ## Set Finite Field Polynomial (for reduction)
    print("- Product of the finite field defined by x^3+x+1")
    GF2.set_polynomial("(x^2+1)*x+1") ;
    a = GF2("(x+1)*x+x/x") ;
    b = a*a ;
    print("a: "+str(a)) ;
    print("b: "+str(b)) ;
    print(a.deg) ; # get the degree
    print(a**-1) ; # compute a^{-1}
    print(~a) ;  # compute the inverse (based on Euclidean Algorithm)
    print(GF2.gcd("x^3+x^2+x+1", "x^3+x")) ;



    ## By Setting the F to zero, one can perform polynomial operations

    print("- Polynomial Arithmetic")
    GF2.set_polynomial(0) ;
    print("a*b:"+str(a*b)) ;
    print("(a*b)%(x^2+x+1): "+str((a*b)%GF2("x^2+x+1"))) ;
    print("x^10%(x^3+1): "+str(GF2("x^10")%GF2("x^3+1"))) ;

    print(GF2("x^5+x+1")//GF2("x^3+1")) ;
    print(GF2("x^5+x+1") % GF2("x^3+1"));

    c = GF2("x^7+x^5+x^4+x^3+x^2+1") ;
    print(c.derivate()) ;

    print(" a == a ? :"+str(a==a)) ;
    print(" a != a ? :" + str(a != a));

    GF2.set_polynomial("x^4+x+1") ;


    print(GF2(GF2.gcd("x^6+x^5+x^4+x^2+x+1", "x^5+x^4+x^2+x"))) ;
    print(GF2(GF2.gcd("x^3+x","x^3+x+1")));

