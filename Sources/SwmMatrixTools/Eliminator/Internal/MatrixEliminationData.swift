//
//  RowAlignedMatrixData.swift
//  SwiftyMath
//
//  Created by Taketo Sano on 2019/11/02.
//

import SwmCore

internal final class MatrixEliminationData<R: Ring> {
    typealias Row = LinkedList<RowEntry<R>>
    
    private(set) var size: MatrixSize
    private(set) var rows: [Row]
    private(set) var rowWeights: [Double]
    private var tracker: RowHeadTracker!

    init<S: Sequence>(size: MatrixSize, entries: S) where S.Element == MatrixEntry<R> {
        self.size = size
        self.rows = []
        self.rowWeights = []
        self.tracker = RowHeadTracker([])
        self.setup(size: size, entries: entries)
    }
    
    convenience init<Impl, n, m>(_ A: MatrixIF<Impl, n, m>) where Impl.BaseRing == R {
        self.init(size: A.size, entries: A.nonZeroEntries)
    }
    
    private func setup<S>(size: MatrixSize, entries: S)
    where S: Sequence, S.Element == MatrixEntry<R> {
        self.size = size
        
        let group = entries.group{ c in c.row }
        self.rows = (0 ..< size.rows).map { i in
            if let list = group[i] {
                let sorted = list.map{ c in RowEntry(c.col, c.value) }.sorted{ $0.col }
                return Row(sorted)
            } else {
                return Row()
            }
        }
        
        self.rowWeights = rows.map { row in
            row.sum{ weight(of: $0.value) }
        }
        self.tracker = RowHeadTracker(headEntries)
    }
    
    @inlinable
    func row(_ i: Int) -> Row {
        rows[i]
    }
    
    @inlinable
    func weight(of a: R) -> Double {
        a.computationalWeight
    }
    
    @inlinable
    func rowWeight(_ i: Int) -> Double {
        rowWeights[i]
    }
    
    var allEntries: [MatrixEntry<R>] {
        rows.enumerated().flatMap { (i, row) in
            row.map { (j, a) in (i, j, a) }
        }
    }
    
    var headEntries: [MatrixEntry<R>] {
        rows.enumerated().compactMap { (i, row) in
            row.headElement.map { (j, a) in
                (i, j, a)
            }
        }
    }
    
    // rename: colEntries(withRowHeadInCol j: Int)
    func colEntries(withHeadInCol j: Int) -> [ColEntry<R>] {
        tracker
            .rows(withHeadInCol: j)
            .map{ i in (i, row(i).headElement!.value) }
    }
    
    func colEntries(in j0: Int, aboveRow i0: Int) -> [ColEntry<R>] {
        (0 ..< i0).compactMap { i -> ColEntry<R>? in
            if let a = find(i, j0).hit?.pointee.element.value {
                return (i, a)
            } else {
                return nil
            }
        }
    }
    
    func find(_ i: Int, _ j: Int) -> (hit: Row.NodePointer?, prev: Row.NodePointer?) {
        rows[i].find({ e in e.col == j}, while: { e in e.col <= j})
    }
    
    func transpose() {
        setup(
            size: (size.cols, size.rows),
            entries: allEntries.map { (i, j, a) in (j, i, a) }
        )
    }
    
    func resultAs<Impl, n, m>(_ type: MatrixIF<Impl, n, m>.Type) -> MatrixIF<Impl, n, m> where Impl.BaseRing == R {
        .init(size: size) { setEntry in
            for (i, row) in rows.enumerated() {
                for e in row {
                    setEntry(i, e.col, e.value)
                }
            }
        }
    }
    
    // only for debug
    var density: Double {
        let nnz = rowWeights.sum()
        return nnz == 0 ? 0 : Double(nnz) / Double(size.rows * size.cols)
    }
    
    @discardableResult
    func apply(_ s: RowElementaryOperation<R>) -> Self {
        switch s {
        case let .AddRow(i, j, r):
            addRow(at: i, to: j, multipliedBy: r)
        case let .MulRow(i, r):
            multiplyRow(at: i, by: r)
        case let .SwapRows(i, j):
            swapRows(i, j)
        }
        return self
    }

    @discardableResult
    func applyAll<S>(_ seq: S) -> Self
    where S: Sequence, S.Element == RowElementaryOperation<R> {
        for s in seq {
            apply(s)
        }
        return self
    }

    @discardableResult
    func addRow(at i1: Int, to i2: Int, multipliedBy r: R) -> Self {
        if row(i1).isEmpty {
            return self
        }
        
        let j2 = row(i2).headElement?.col
        let w = addRow(row(i1), into: row(i2), multipliedBy: r)
        
        rowWeights[i2] += w
        tracker.headMoved(inRow: i2, fromCol: j2, toCol: row(i2).headElement?.col)
        
        return self
    }
    
    @discardableResult
    func batchAddRow(at i1: Int, to rows: [Int], multipliedBy rs: [R]) -> Self {
        let from = row(i1)
        if from.isEmpty {
            return self
        }
        
        let oldCols = rows.map{ i in
            row(i).headElement?.col
        }
        
        let weights = zip(rows, rs)
            .map{ (i2, r) in (row(i2), r)}
            .parallelMap { (to, r) in
                addRow(from, into: to, multipliedBy: r)
            }
        
        for (i, w) in zip(rows, weights) {
            rowWeights[i] += w
        }
        for (i, j) in zip(rows, oldCols) {
            tracker.headMoved(inRow: i, fromCol: j, toCol: row(i).headElement?.col)
        }
        
        return self
    }
    
    @discardableResult
    @_specialize(where R == 𝐙)
    @_specialize(where R == 𝐐)
    @_specialize(where R == 𝐅₂)
    
    private func addRow(_ from: Row, into to: Row, multipliedBy r: R) -> Double {
        if from.isEmpty {
            return 0
        }

        var w = 0.0
        
        let fromHeadCol = from.headElement!.col
        if to.isEmpty || fromHeadCol < to.headElement!.col {
            
            // from: ●-->○-->○----->○-------->
            //   to:            ●------->○--->
            //
            //   ↓
            //
            // from: ●-->○-->○----->○-------->
            //   to: ●--------->○------->○--->
            
            to.insertHead( (fromHeadCol, .zero) )
        }
        
        var fromItr = from.makeIterator()
        var toPtr = to.headPointer!
        var toPrevPtr = toPtr
        
        while let (j1, a1) = fromItr.next() {
            // At this point, it is assured that `from.value.col >= to.value.col`.
            // Proceed `to` so that it comes closed to `from`.
            
            // from: ------------->●--->○-------->
            //   to: -->○----->●------------>○--->
            //
            //   ↓
            //
            // from: ------------->●--->○-------->
            //   to: -->○----->●------------>○--->

            while let next = toPtr.pointee.next, next.pointee.element.col <= j1 {
                (toPrevPtr, toPtr) = (toPtr, next)
            }
            
            let (j2, a2) = toPtr.pointee.element
            
            if j1 == j2 {
                //                     j1 = j2
                // from: ------------->●--->○-------->
                //   to: -->○--------->●-------->○--->

                let b = a2 + r * a1
                
                if b.isZero && toPtr != toPrevPtr {
                    toPtr = toPrevPtr
                    toPtr.pointee.dropNext()
                } else {
                    toPtr.pointee.element.value = b
                }
                
                w += weight(of: b) - weight(of: a2)
                
            } else {
                //                 j2  j1
                // from: ------------->●--->○-------->
                //   to: -->○----->●------------>○--->

                let b = r * a1
                toPtr.pointee.insertNext( RowEntry(j1, b) )
                (toPrevPtr, toPtr) = (toPtr, toPtr.pointee.next!)
                
                w += weight(of: b)
            }
        }
        
        if to.headElement!.value.isZero {
            to.dropHead()
        }
        
        return w
    }
    
    @discardableResult
    @_specialize(where R == 𝐙)
    @_specialize(where R == 𝐐)
    @_specialize(where R == 𝐅₂)
    
    func multiplyRow(at i: Int, by r: R) -> Self {
        row(i).modifyEach { e in
            e.value = r * e.value
        }
        return self
    }
    
    @discardableResult
    func swapRows(_ i1: Int, _ i2: Int) -> Self {
        let j1 = row(i1).headElement?.col
        let j2 = row(i2).headElement?.col
        
        rows.swapAt(i1, i2)
        rowWeights.swapAt(i1, i2)
        tracker.swapRows( (i1, j1), (i2, j2) )
        
        return self
    }
    
    private final class RowHeadTracker {
        private var rowHeads: [Int : Set<Int>] // [col : Set<rows>]

        init<S>(_ headEntries: S) where S: Sequence, S.Element == MatrixEntry<R> {
            self.rowHeads = [:]
            for (i, j, _) in headEntries {
                rowHeads[j, default: []].insert(i)
            }
        }
        
        func rows(withHeadInCol j: Int) -> [Int] {
            rowHeads[j, default: []].sorted()
        }
        
        func swapRows(_ e1: (Int, Int?), _ e2: (Int, Int?)) {
            let (i1, j1) = e1
            let (i2, j2) = e2
            
            if j1 == j2 {
                return
            }
            
            if let j1 = j1 {
                rowHeads[j1]?.remove(i1)
                rowHeads[j1, default: []].insert(i2)
            }
            
            if let j2 = j2 {
                rowHeads[j2]?.remove(i2)
                rowHeads[j2, default: []].insert(i1)
            }
        }
        
        func headMoved(inRow i: Int, fromCol j0: Int?, toCol j1: Int?) {
            if j0 == j1 {
                return
            }
            
            if let j0 = j0 {
                rowHeads[j0]?.remove(i)
            }
            if let j1 = j1 {
                rowHeads[j1, default: []].insert(i)
            }
        }
    }
}

