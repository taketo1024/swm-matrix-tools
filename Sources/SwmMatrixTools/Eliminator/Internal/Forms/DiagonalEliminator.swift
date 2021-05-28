//
//  DiagonalEliminator.swift
//  SwiftyMath
//
//  Created by Taketo Sano on 2017/11/08.
//  Copyright © 2017年 Taketo Sano. All rights reserved.
//

import SwmCore

internal final class DiagonalEliminator<R: EuclideanRing>: MatrixEliminator<R> {
    override var form: MatrixEliminationForm {
        .Diagonal
    }
    
    override func isDone() -> Bool {
        data.entries.allSatisfy { (i, j, a) in
            (i == j) && a.isNormalized
        }
    }
    
    override func iteration() {
        subrun(RowEchelonEliminator<R>.self)
        
        if isDone() {
            return
        }
        
        subrun(RowEchelonEliminator<R>.self, transpose: true)
    }
}
