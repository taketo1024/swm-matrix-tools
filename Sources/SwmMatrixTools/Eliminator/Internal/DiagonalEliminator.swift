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
        data.rows.enumerated().allSatisfy { (i, row) in
            row.isEmpty || row.head.flatMap { ptr -> Bool in
                let (j, a) = ptr.element
                return i == j && a.isNormalized && !ptr.hasNext
            } ?? false
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
