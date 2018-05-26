# Binary_Field

The project includes several classes.
1. GF2: Allows to perform polynomial/finite field arithmetics in characteristic 2.
2. Binary_Matrix: Allows to perform matrix arithmetic in characteristic 2.
3. Berlekamp: Allows to factorize binary polynomial thanks to the Berlekamp algorithm.
4. Mastrovito: A FPGA Design Tool to generate Mastrovito Matrices.


## GF2 Class

The class can be used to perform arithmetic over any binary field.
The "reduction" polynomial (say __F__)  can be set through the *set_polynomial* method.

```python
GF2.set_polynomial("x^3+x+1") 
```

Hereafter, one can create GF2 elements with the GF2 class constructor,

```python
a = GF2("x+1") ;
b = GF2("x^2") ;
```
and can run products, additions, modular exponentiations, reductions and inversions in a very natural and intuitive way. For instance, 
```python
c = a*b ; # (x+1)*(x^2) mod (x^3+x+1) = (x^3+x^2) mod (x^3+x+1) = x^2+x+1
```
should provide the correct output.

Moreover, setting the __F__ to 0 allows the user to perform calculations over the polynomial ring (mod 2).
Note the constructor should be provided "expanded form" of polynomial since there is no algorithm to carry out more complex input (for instance, *GF2("(x+1)*(x+1)")* is not allowed, one should first define *A=GF2("x+1")* and then compute *A\*A*).

Notice the difference between *A/B* and *A//B*. The former actually computes *A\*B^{-1}* if the inverse of B exists while the latter computes A/B interpreted as polynomials *(i.e. (x^3+x+1)//x == x^2+1)* .

Numerous examples are given in the "main body" and should show the user how to handle the class functionalities.

## Binary Matrix Class

This class allows to perform arithmetic over binary matrices. One can perform multiplications, inversions, evaluate the determinant, get a triangular form (based on the Gauss Elimination method) as well as computing kernels of linear applications. 

To create a matrix, simply write it as a list of "rows".
```python
M = Binary_matrix([[1,0,1], [0,1,1]]) ; 
```

Once matrix objects are created, one can manipulated them as if they were simple integers
```python
M2 = M*M.transpose() ; 
M3 = M**5 ;
M4 = M3+M ; 
```

If a matrix is invertible (i.e. determinant differs from 0), one can compute its inverse thanks to the ~ symbol.

```python
M2_inv = ~M_inv ; 
```

One can check that 

```python
Id = M2*M_inv ; 
```
is equal to the identity matrix indeed. The identity matrix can be instantiated thanks to the class method *Binary_matrix.eye(n)* where n is the dimension of the identity matrix (i.e. square matrix with n^2 elements). The zero matrix can be created as well through the *Binary_matrix.zero(m,n)* method which spans a (m,n) matrix filled up with 0.

Kernel of a matrix can be computed with the class method *Binary_matrix.ker* which returns a matrix whose rows are those which span the ker seen as a vector space.

Note the code is not fully optimized as some code parts are run several times while results could have been stored. Nonetheless, this code is still very usable and somehow easy to understand.

Numerous examples are given in the "main body" and should demonstrate to the user how to handle the class functionalities.

## Berlekamp Algorithm

This class includes tool to factorize binary polynomials (i.e. whose coefficients belong to __GF(2)__).
Actually, only one method should be called (namely *Berlekamp.Factorize*) as the others could be considered as "private" (see the main body for example). This class relies on the Berlekamp algorithm which mainly transforms the factorization problem into a linear algebra one (i.e. finding the kernel of a linear application). An example of the Factorisation method usage is given below.

```python
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
```
## Mastrovito Class (FPGA Design Tool)

For now, this class allows to generate a *.vhd* file which matches the Mastrovito Matrix structure one would generate for its FPGA design. This matrix is used to perform a finite field multiplication in one clock cycle (i.e. pure combinatatioral circuit) and could easily be "plugged" in your FPGA architecture. For the *x^4+x+1* reduction polynomial, the Mastrovito matrix would look like the following.

```vhld
c(0) <= (b(0) and (a(0))) xor (b(1) and (a(3))) xor (b(2) and (a(2))) xor (b(3) and (a(1))) ; 
c(1) <= (b(0) and (a(1))) xor (b(1) and (a(0) xor a(3))) xor (b(2) and (a(2) xor a(3))) xor (b(3) and (a(1) xor a(2))) ; 
c(2) <= (b(0) and (a(2))) xor (b(1) and (a(1))) xor (b(2) and (a(0) xor a(3))) xor (b(3) and (a(2) xor a(3))) ; 
c(3) <= (b(0) and (a(3))) xor (b(1) and (a(2))) xor (b(2) and (a(1))) xor (b(3) and (a(0) xor a(3))) ; 

```

## Setup

Copy the source files into your project folder and import the desired classes. For instance,
```python
from GF2 import GF2
```


## Note

The code has not been extensively tested so bugs may remain. If you meet any, don't hesite to contact me. 
