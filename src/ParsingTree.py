
class ParsingTree(object):


    ## Encoding transformation rule
    ## A*(B+C) => A*B + A*C
    ## (A*B)*C => (A*C)*B

    operators = ["+", "-", "*", "/"] ; # sort in priority order.

    def __init__(self, str_expression = None):
        self.node = "" ;
        self.right = None ;
        self.left = None ;
        if str_expression:
            ParsingTree.parse(self, str_expression) ;

    def eval(self, _class = None):
        return ParsingTree.evaluate(self, _class) ;


    @classmethod
    def split(cls, str_expression, op, max_split = None):
        bracket_cpt = 0 ;
        ret = [] ;
        for i in range(len(str_expression)):
            if (str_expression[i] == '('):
                bracket_cpt += 1 ;
            elif (str_expression[i] == ')'):
                bracket_cpt -= 1;
            elif (bracket_cpt == 0 and str_expression[i] == op):
                ret.append(str_expression[:i]) ;
                if max_split != None and max_split == len(ret):
                    ret.append(str_expression[i+1:]) ;
                    break ;
        if not(len(ret)):
            ret.append(str_expression) ;
        return ret ;

    @classmethod
    def remove_brackets(cls, str_expression):
        if (str_expression[0] == '(' and str_expression[-1] == ')'):
            return str_expression[1:-1] ;
        else:
            return "" ;

    @classmethod
    def parse(cls, root, str_expression):
        for op in ParsingTree.operators:
            tmp = ParsingTree.split(str_expression, op, 1) ;
            if (len(tmp) == 2): break ;
        if (len(tmp) == 2):
            left = tmp[0] ;
            right = tmp[1] ;
            root.node = op ;
            root.left = ParsingTree() ; ParsingTree.parse(root.left, left) ;
            root.right = ParsingTree() ; ParsingTree.parse(root.right, right);
        else:
            tmp_expression = ParsingTree.remove_brackets(str_expression) ;
            if (len(tmp_expression) == 0):
                    root.node = str_expression ;
            else:
                ParsingTree.parse(root, tmp_expression) ;

    @classmethod
    def evaluate(cls, root, _class = None):
        if (root.left == None or root.right == None):
            if _class:
                return _class(root.node) ;
            else:
                return int(root.node) ;
        elif root.node == "+":
            return ParsingTree.evaluate(root.left, _class)+ParsingTree.evaluate(root.right, _class) ;
        elif root.node == "-":
            return ParsingTree.evaluate(root.left, _class) - ParsingTree.evaluate(root.right, _class);
        elif root.node == "*":
            return ParsingTree.evaluate(root.left, _class) * ParsingTree.evaluate(root.right, _class);
        elif root.node == "/":
            return ParsingTree.evaluate(root.left, _class) / ParsingTree.evaluate(root.right, _class);





if __name__ == "__main__":
    a = "((2*1+2)*5+1)" ;
    exp = ParsingTree(a) ;
    print(exp.eval()) ;