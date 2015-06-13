# Index of functions in my_r_func.R

ParseGSE("GSEnumber") 
    return list("exp".matrix, "pat".matrix, "prob".matrix)

GroupPat(pat.matrix) 
    return same matrix with one coloum "group" contains 'G0','G1'

Preproc(exp.matrix), 
    return same matrix.

MatrixTest(x,y,manner=1(if x vs y) /2(if x~y) ), 
    return matrix. probs vs. categories(t.adj.p,w.p.adj)

SelectThres(matrix, p.val-colname, thres=0.05)
    return matrix with p.val < thres

ProbsToGene(list(list(exp.matrix, prob.matrix),
                 list(exp.matrix, prob.matrix)))
    return list(exp.dataframe(add "geneID" coloum to the last),
                exp.dataframe(second matrix))

SelectSpec(list.contains.expression.matrix.as.above),
           "colname which contains geneID(e.x. 'geneID')")
    return list("mat".gene distribution matrix,
                "summ".summary matrix)
	
S1Parse("GSEnumber")
    return list("exp".matrix, "prob".matrix, "pat".matrix,
                "test".matrix, "sig".matrix)

S2MultiCanc(list(cancer.name=cancer.obj, caner.name=cancer.obj2))
    return list("mat".matrix, "summ".matrix, "cancsum".matrix)
