//
//  DiagonalEliminator.swift
//  SwiftyMath
//
//  Created by Taketo Sano on 2017/11/08.
//  Copyright © 2017年 Taketo Sano. All rights reserved.
//

import SwmCore

internal final class DiagonalEliminator<R>: MatrixEliminator<R>
where R: EuclideanRing & ComputationalRing {
    override var form: MatrixEliminationForm {
        .Diagonal
    }
    
    override func isDone() -> Bool {
        data.rows.enumerated().allSatisfy { (i, row) in
            if let head = row.head {
                return row.isSingle && head.col == i && head.value.isNormalized
            } else {
                return true
            }
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
