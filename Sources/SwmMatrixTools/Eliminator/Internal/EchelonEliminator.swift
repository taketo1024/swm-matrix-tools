//
//  EchelonEliminator.swift
//  SwiftyMath
//
//  Created by Taketo Sano on 2017/11/08.
//  Copyright © 2017年 Taketo Sano. All rights reserved.
//

import SwmCore

internal final class RowEchelonEliminator<R: EuclideanRing>: MatrixEliminator<R> {
    var currentRow = 0
    var currentCol = 0
    
    override var form: MatrixEliminationForm {
        .RowEchelon
    }
    
    override func isDone() -> Bool {
        currentRow >= data.size.rows || currentCol >= data.size.cols
    }
    
    @_specialize(where R == 𝐙)
    override func iteration() {
        
        // find pivot point
        let colEntries = data.colEntries(withHeadInCol: currentCol)
        guard let pivot = findPivot(in: colEntries) else {
            currentCol += 1
            return
        }
        
        log("Pivot: \((pivot.row, currentCol)), \(pivot.value)")
        
        let again = eliminate(colEntries, byPivot: pivot)
        if again {
            return
        }
        
        // final step
        if !pivot.value.isNormalized {
            apply(.MulRow(at: pivot.row, by: pivot.value.normalizingUnit))
        }
        
        if pivot.row != currentRow {
            apply(.SwapRows(pivot.row, currentRow))
        }
        
        currentRow += 1
        currentCol += 1
    }
    
    @_specialize(where R == 𝐙)
    private func findPivot(in candidates: [ColEntry<R>]) -> ColEntry<R>? {
        candidates.min { (c1, c2) in
            let (i1, i2) = (c1.row, c2.row)
            let (d1, d2) = (c1.value.euclideanDegree, c2.value.euclideanDegree)
            return d1 < d2 || (d1 == d2 && data.rowWeight(i1) < data.rowWeight(i2))
        }
    }
    
    private func eliminate(_ colEntries: [ColEntry<R>], byPivot pivot: ColEntry<R>) -> Bool {
        var again = false
        let targets = colEntries
            .exclude{ $0.row == pivot.row }
            .map { (i, a) -> (Int, R) in
                let (q, r) = a /% pivot.value
                if !r.isZero {
                    again = true
                }
                return (i, -q)
            }
        
        addRow(at: pivot.row, to: targets)
        
        return again
    }
    
    private func addRow(at i0: Int, to: [(Int, R)]) {
        data.addRow(at: i0, to: to)
        append(to.map{ (i, r) in
            .AddRow(at: i0, to: i, mul: r)
        })
    }
}
