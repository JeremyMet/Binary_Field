# Binary_Field

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

Kernel of a matrix can be computed with the class method *Binary_matrix.ker* which returns a matrix whose rows are those which span the ker seen as a vectorial space.

Note the code is not fully optimized as some code parts are run several times while results could have been stored. Nonetheless, this code is still very usable and somehow easy to understand.

Numerous examples are given in the "main body" and should demonstrate to the user how to handle the class functionalities.

## Berlekamp Algorithm

This class includes tool to factorize binary polynomials (i.e. whose coefficients belong to __GF(2)__).
Actually, only one method should be called (namely *Berlekamp.Factorize*) as the others could be considered as "private" (see the main body for example). This class relies on the Berlekamp algorithm which mainly transforms the factorization problem into a linear algebra one.
