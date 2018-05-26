from GF2 import GF2 ;


class Mastrovito(object):

    xor_str = "xor" ;
    and_str = "and"

    @classmethod
    def to_vhd(cls, F):
        F = GF2.build_from_string(F) ;
        F_DEG = GF2.degree(F) ;
        current_column = [(1 << i) for i in range(F_DEG)] ;
        output_string = ["c("+str(i)+") <= " for i in range(F_DEG)] ;
        ## First, get indexes non zero element of F (it will allow faster reduction then)
        cpy_F = F ; cpt = 0 ; F_idx = [] ;
        while(cpy_F > 0):
            if cpy_F&1:
                F_idx.append(cpt) ;
            cpt += 1;
            cpy_F = cpy_F >> 1 ;
        F_idx = F_idx[:-1] ;
        F_idx_len = len(F_idx) ;
        ## Thus, build the matrix
        for i in range(F_DEG):
            ## Updating output string
            for j in range(F_DEG):
                if current_column[j]:
                    output_string[j] += "(b("+str(i)+") "+Mastrovito.and_str+" ("
                    for k in range(F_DEG):
                        if (current_column[j] & (1 << k)):
                            output_string[j] += "a("+str(k)+") "+Mastrovito.xor_str+" " ;
                output_string[j] = output_string[j][:-5]+"))" ;
                if i+1 != F_DEG:
                    output_string[j] += " "+Mastrovito.xor_str+" " ; 
            ## Shift
            current_column.insert(0, 0) ;
            val = current_column.pop(F_DEG) ;
            ## Reduction
            for idx in F_idx:
                current_column[idx] = current_column[idx]^val ;
        ## Convert output_string into full string
        ret = "" ;
        for s in output_string:
            ret+=s+" ; \n" ;
        return ret[:-1]


if __name__ == "__main__":
    pol = "x^4+x+1" ;
    print(Mastrovito.to_vhd(pol)) ;