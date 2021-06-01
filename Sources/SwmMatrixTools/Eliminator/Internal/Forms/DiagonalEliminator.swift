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
            row.isEmpty || row.headPointer.flatMap { ptr -> Bool in
                let node = ptr.pointee
                let (j, a) = node.element
                return i == j && a.isNormalized && !node.hasNext
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
