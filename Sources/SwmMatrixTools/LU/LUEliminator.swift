//
//  File.swift
//  
//
//  Created by Taketo Sano on 2021/06/02.
//

import SwmCore

internal class LUEliminator<R: Ring>: MatrixEliminator<R> {
    private var currentRow = 0
    private var currentCol = 0
    private var tmpL: [[RowEntry<R>]]
    
    required init(data: MatrixEliminationData<R>) {
        self.tmpL = Array(repeating: [RowEntry<R>].empty, count: data.size.rows)
        super.init(data: data)
    }
    
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
        
        if !pivot.value.isInvertible {
            abort()
            return
        }
        
        if pivot.row != currentRow {
            apply(.SwapRows(currentRow, pivot.row))
            tmpL.swapAt(currentRow, pivot.row)
            return
        }
        
        log("Pivot: \((pivot.row, currentCol)), \(pivot.value)")
        
        eliminate(colEntries, byPivot: pivot)
        
        currentRow += 1
        currentCol += 1
    }
    
    @_specialize(where R == 𝐙)
    private func findPivot(in candidates: [ColEntry<R>]) -> ColEntry<R>? {
        candidates.min { (c1, c2) in
            let (i1, i2) = (c1.row, c2.row)
            let (d1, d2) = (c1.value.isInvertible ? 0 : 1, c2.value.isInvertible ? 0 : 1)
            return d1 < d2 || (d1 == d2 && data.rowWeight(i1) < data.rowWeight(i2))
        }
    }
    
    private func eliminate(_ colEntries: [ColEntry<R>], byPivot pivot: ColEntry<R>) {
        let i0 = pivot.row
        let e = pivot.value.inverse!
        
        let targets = colEntries[1...]
            .map { (i, a) -> (Int, R) in
                (i, -a * e)
            }
        
        addRow(at: pivot.row, to: targets)
        
        tmpL[i0].append((i0, .identity))
        for (i, a) in targets {
            tmpL[i].append( (i0, -a) )
        }
    }
    
    fileprivate func addRow(at i: Int, to: [(Int, R)]) {
        data.addRow(at: i, to: to)
        append(to.map{ (i1, r) in
            .AddRow(at: i, to: i1, mul: r)
        })
    }
    
    // MARK: PQLU results
    
    func PQLU<M: MatrixImpl>(_ type: M.Type) -> LUFactorizer<M>.Result? where M.BaseRing == R {
        !aborted ? (P, Q, L(M.self), U(M.self)) : nil
    }
    
    private var P: Permutation<anySize> {
        let n = data.size.rows
        var indices = Array(0 ..< n)
        for case let .SwapRows(i, j) in rowOps {
            indices.swapAt(i, j)
        }
        return Permutation(length: n, indices: indices).inverse!
    }
    
    private var Q: Permutation<anySize> {
        let m = data.size.cols
        let indices = data.headEntries.map { $0.col }
        return Permutation(length: m, indices: indices, fillRemaining: true).inverse!
    }
    
    private func L<M: MatrixImpl>(_ type: M.Type) -> M where M.BaseRing == R {
        let n = data.size.rows
        let r = data.countNonZeroRows
        return .init(
            size: (n, r),
            entries: tmpL.enumerated().flatMap{ (i, row) in
                row.map { (j, a) in (i, j, a) }
            }
        )
    }
    
    private func U<M: MatrixImpl>(_ type: M.Type) -> M where M.BaseRing == R {
        let m = data.size.cols
        let r = data.countNonZeroRows
        return data.resultAs(M.self).permuteCols(by: Q).submatrix(rowRange: 0 ..< r, colRange: 0 ..< m)
    }
}
