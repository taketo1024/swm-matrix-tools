//
//  File.swift
//  
//
//  Created by Taketo Sano on 2021/06/02.
//

import SwmCore

public final class LUFactorizer<M: LUFactorizable> {
    public typealias R = M.BaseRing
    private let eliminator: LUEliminator
    
    public init(_ A: M) {
        self.eliminator = LUEliminator(data: MatrixEliminationData(A))
    }

    public convenience init<n, m>(_ A: MatrixIF<M, n, m>) {
        self.init(A.impl)
    }
    
    public func run() {
        eliminator.run()
    }
    
    public var isDone: Bool {
        eliminator.isDone()
    }
    
    // MARK: PQLU results
    
    public var result: M.RawLUFactorizationResult {
        assert(isDone)
        return (P, Q, L, U)
    }
    
    public var P: Permutation<anySize> {
        let data = eliminator.data
        let n = data.size.rows
        var indices = Array(0 ..< n)
        for case let .SwapRows(i, j) in eliminator.rowOps {
            indices.swapAt(i, j)
        }
        return Permutation(length: n, indices: indices).inverse!
    }
    
    public var Q: Permutation<anySize> {
        let data = eliminator.data
        let m = data.size.cols
        let indices = data.headEntries.map { $0.col }
        return Permutation(length: m, indices: indices, fillRemaining: true).inverse!
    }
    
    public var L: M {
        let data = eliminator.data
        let n = data.size.rows
        let r = data.countNonZeroRows
        return .init(
            size: (n, r),
            entries: eliminator.tmpL.enumerated().flatMap{ (i, row) in
                row.map { (j, a) in (i, j, a) }
            }
        )
    }
    
    public var U: M {
        let data = eliminator.data
        let m = data.size.cols
        let r = data.countNonZeroRows
        return data.resultAs(M.self).permuteCols(by: Q).submatrix(rowRange: 0 ..< r, colRange: 0 ..< m)
    }
    
    private class LUEliminator: MatrixEliminator<R> {
        private var currentRow = 0
        private var currentCol = 0
        fileprivate var tmpL: [[RowEntry<R>]]
        
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
        
        private func addRow(at i: Int, to: [(Int, R)]) {
            data.addRow(at: i, to: to)
            append(to.map{ (i1, r) in
                .AddRow(at: i, to: i1, mul: r)
            })
        }
    }
}
