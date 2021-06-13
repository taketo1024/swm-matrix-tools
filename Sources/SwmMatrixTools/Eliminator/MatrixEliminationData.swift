//
//  RowAlignedMatrixData.swift
//  SwiftyMath
//
//  Created by Taketo Sano on 2019/11/02.
//

import SwmCore

internal struct MatrixEliminationData<R> where R: ComputationalRing {
    private(set) var size: MatrixSize
    private(set) var rows: [Row]
    private(set) var rowWeights: [Double]
    private var tracker: RowHeadTracker!

    init<S: Sequence>(size: MatrixSize, entries: S, transpose: Bool) where S.Element == MatrixEntry<R> {
        self.size = (0, 0)
        self.rows = []
        self.rowWeights = []
        self.tracker = RowHeadTracker([])
        self.setup(size: size, entries: entries, transpose: transpose)
    }
    
    init<Impl: MatrixImpl>(_ A: Impl, transpose: Bool = false) where Impl.BaseRing == R {
        self.init(size: A.size, entries: A.nonZeroEntries, transpose: transpose)
    }
    
    init<Impl, n, m>(_ A: MatrixIF<Impl, n, m>, transpose: Bool = false) where Impl.BaseRing == R {
        self.init(A.impl, transpose: transpose)
    }
    
    private mutating func setup<S>(size: MatrixSize, entries: S, transpose: Bool)
    where S: Sequence, S.Element == MatrixEntry<R> {
        if transpose {
            // TODO generate tranposed rows directly
            setup(
                size: (size.cols, size.rows),
                entries: entries.map { (i, j, a) in (j, i, a) },
                transpose: false
            )
            return
        }
        
        self.size = size
        self.rows = generateRows(entries)
        self.rowWeights = rows.map { row in
            row.sum{ $0.value.computationalWeight }
        }
        self.tracker = RowHeadTracker(headEntries)
    }
    
    private func generateRows<S>(_ entries: S) -> [Row]
    where S: Sequence, S.Element == MatrixEntry<R> {
        let group = entries.group{ c in c.row }
        return (0 ..< size.rows).map { i in
            let sorted = (group[i] ?? []).map{ (_, j, a) in
                RowEntry(j, a)
            }.sorted{ $0.col }
            return Row(sorted)
        }
    }
    
    static func empty(size: MatrixSize) -> Self {
        .init(size: size, entries: [], transpose: false)
    }
    
    @inlinable
    var countNonZeroRows: Int {
        rows.count{ !$0.isEmpty }
    }
    
    @inlinable
    func row(_ i: Int) -> Row {
        rows[i]
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
            row.head.flatMap{ (j, a) in
                (i, j, a)
            }
        }
    }
    
    // rename: colEntries(withRowHeadInCol j: Int)
    func colEntries(withHeadInCol j: Int) -> [ColEntry<R>] {
        tracker
            .rows(withHeadInCol: j)
            .map{ i in (i, row(i).head!.value) }
    }
    
    mutating func transpose() {
        setup(size: size, entries: allEntries, transpose: true)
    }
    
    func resultAs<M: MatrixImpl>(_ type: M.Type) -> M where M.BaseRing == R {
        .init(size: size) { setEntry in
            for (i, row) in rows.enumerated() {
                for e in row {
                    setEntry(i, e.col, e.value)
                }
            }
        }
    }
    
    func resultAs<M, n, m>(_ type: MatrixIF<M, n, m>.Type) -> MatrixIF<M, n, m> where M.BaseRing == R {
        .init(resultAs(M.self))
    }
    
    // only for debug
    var density: Double {
        let nnz = rowWeights.sum()
        return nnz == 0 ? 0 : Double(nnz) / Double(size.rows * size.cols)
    }
    
    @discardableResult
    mutating func apply(_ s: RowElementaryOperation<R>) -> Self {
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
    mutating func applyAll<S>(_ seq: S) -> Self
    where S: Sequence, S.Element == RowElementaryOperation<R> {
        for s in seq {
            apply(s)
        }
        return self
    }

    @discardableResult
    mutating func addRow(at i1: Int, to i2: Int, multipliedBy r: R) -> Self {
        let from = row(i1)
        if from.isEmpty {
            return self
        }
        
        let oldCol = row(i2).head?.col
        let w = Self.addRow(from, into: &rows[i2], multipliedBy: r)
        
        rowWeights[i2] += w
        tracker.headMoved(inRow: i2, fromCol: oldCol, toCol: row(i2).head?.col)
        
        return self
    }
    
    @discardableResult
    mutating func addRow(at i1: Int, to: [(row: Int, multipliedBy: R)]) -> Self {
        let from = row(i1)
        if from.isEmpty {
            return self
        }
        
        let targetRows = to.map{ $0.row }
        let oldCols = targetRows.map{ i in
            row(i).head?.col
        }
        
        var weights: [Double] = []
        rows.withUnsafeMutableBufferPointer { buff in
            weights = to.parallelMap { (i, r) in
                Self.addRow(from, into: &buff[i], multipliedBy: r)
            }
        }
        
        for (i, w) in zip(targetRows, weights) {
            rowWeights[i] += w
        }
        
        for (i, j) in zip(targetRows, oldCols) {
            tracker.headMoved(inRow: i, fromCol: j, toCol: row(i).head?.col)
        }
        
        return self
    }
    
    @discardableResult
    @_specialize(where R == 𝐙)
    @_specialize(where R == 𝐐)
    @_specialize(where R == 𝐅₂)
    
    private static func addRow(_ from: Row, into to: inout Row, multipliedBy r: R) -> Double {
        if from.isEmpty {
            return 0
        }

        var w = 0.0
        
        let fromHeadCol = from.head!.col
        if to.isEmpty || fromHeadCol < to.head!.col {
            
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
        
        to.modify { toPtr in
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
                
                while toPtr.hasNext && toPtr.next.col <= j1 {
                    toPrevPtr = toPtr
                    toPtr.proceed()
                }
                
                let (j2, a2) = (toPtr.col, toPtr.value)
                
                if j1 == j2 {
                    //                     j1 = j2
                    // from: ------------->●--->○-------->
                    //   to: -->○--------->●-------->○--->
                    
                    let b = a2 + r * a1
                    
                    if b.isZero && toPtr != toPrevPtr {
                        toPrevPtr.dropNext()
                        toPtr = toPrevPtr
                    } else {
                        toPtr.value = b
                    }
                    
                    w += b.computationalWeight - a2.computationalWeight
                    
                } else {
                    //                 j2  j1
                    // from: ------------->●--->○-------->
                    //   to: -->○----->●------------>○--->
                    
                    let b = r * a1
                    toPtr.insertNext( (j1, b) )
                    
                    toPrevPtr = toPtr
                    toPtr.proceed()
                    
                    w += b.computationalWeight
                }
            }
        }
        
        if to.head!.value.isZero {
            to.dropHead()
        }
        
        return w
    }
    
    @_specialize(where R == 𝐙)
    @_specialize(where R == 𝐐)
    @_specialize(where R == 𝐅₂)
    
    @discardableResult
    mutating func multiplyRow(at i: Int, by r: R) -> Self {
        rows[i].mapValues { a in
            r * a
        }
        return self
    }
    
    @discardableResult
    mutating func swapRows(_ i1: Int, _ i2: Int) -> Self {
        let j1 = row(i1).head?.col
        let j2 = row(i2).head?.col
        
        rows.swapAt(i1, i2)
        rowWeights.swapAt(i1, i2)
        tracker.swapRows( (i1, j1), (i2, j2) )
        
        return self
    }
    
    @inlinable
    mutating func modifyRow(_ i: Int, _ modifier: (inout Row.NodePointer) -> Void) {
        rows[i].modify(modifier)
    }
    
    struct Row: Sequence {
        typealias Element = RowEntry<R>
        
        private var headIndex: Int
        private var values: [R]
        private var cols: [Int]
        private var nexts: [Int]
        private var free: Set<Int>
        
        init(_ seq: [RowEntry<R>]) {
            self.headIndex = -1
            self.values = []
            self.cols = []
            self.nexts = []
            self.free = []
            
            let count = seq.count

            values.reserveCapacity(count)
            cols.reserveCapacity(count)
            nexts.reserveCapacity(count)
            free.reserveCapacity(count)
            
            for (j, a) in seq {
                values.append(a)
                cols.append(j)
            }
            
            if count > 0 {
                headIndex = 0
                nexts.append(contentsOf: 1 ..< count)
                nexts.append(-1)
            }
        }
        
        @inlinable
        var isEmpty: Bool {
            headIndex == -1
        }
        
        @inlinable
        var isSingle: Bool {
            headIndex != -1 && nexts[headIndex] == -1
        }
        
        @inlinable
        var head: RowEntry<R>? {
            isEmpty ? nil : (cols[headIndex], values[headIndex])
        }
        
        mutating func insertHead(_ element: RowEntry<R>) {
            let old = headIndex // possibly -1
            headIndex = insert(element)
            nexts[headIndex] = old
        }
        
        mutating func insert(_ element: RowEntry<R>, after: Int) {
            let next = insert(element)
            nexts[next] = nexts[after]
            nexts[after] = next
        }

        mutating func dropHead() {
            assert(!isEmpty)
            let next = nexts[headIndex] // possibly -1
            remove(headIndex)
            headIndex = next
        }
        
        mutating func dropNext(_ index: Int) {
            assert(nexts[index] != -1)
            let next = nexts[index]
            nexts[index] = nexts[next]
            remove(next)
        }

        mutating func mapValues(_ map: (R) -> R) {
            for i in 0 ..< values.count where !free.contains(i) {
                values[i] = map(values[i])
            }
        }
        
        private mutating func insert(_ element: Element) -> Int {
            let index: Int
            if let freed = free.popFirst() {
                index = freed
                (cols[index], values[index]) = element
                nexts[index] = -1
            } else {
                index = values.count
                values.append(element.value)
                cols.append(element.col)
                nexts.append(-1)
            }
            return index
        }
        
        private mutating func remove(_ index: Int) {
            free.insert(index)
        }
        
        @inlinable
        mutating func modify(_ modifier: (inout NodePointer) -> Void) {
            var ptr = NodePointer(&self, headIndex)
            modifier(&ptr)
        }
        
        struct NodePointer: Equatable {
            let row: UnsafeMutablePointer<Row>
            var index: Int
            
            @inlinable
            init(_ row: UnsafeMutablePointer<Row>, _ index: Int) {
                assert(index != -1)
                self.row = row
                self.index = index
            }
            
            @inlinable
            var value: R {
                get {
                    row.pointee.values[index]
                } set  {
                    row.pointee.values[index] = newValue
                }
            }
            
            @inlinable
            var col: Int {
                get {
                    row.pointee.cols[index]
                } set {
                    row.pointee.cols[index] = newValue
                }
            }
            
            @inlinable
            var hasNext: Bool {
                row.pointee.nexts[index] != -1
            }
            
            @inlinable
            var next: NodePointer {
                NodePointer(row, row.pointee.nexts[index])
            }
            
            @inlinable
            func insertNext(_ e: Element) {
                row.pointee.insert(e, after: index)
            }
            
            @inlinable
            func dropNext() {
                row.pointee.dropNext(index)
            }
            
            @inlinable
            mutating func proceed() {
                assert(hasNext)
                index = row.pointee.nexts[index]
            }
            
            @inlinable
            static func == (p: Self, q: Self) -> Bool {
                p.index == q.index
            }
        }
        
        func makeIterator() -> ElementIterator {
            ElementIterator(self, headIndex)
        }

        struct ElementIterator: IteratorProtocol {
            private var row: Row
            private var index: Int
            
            fileprivate init(_ row: Row, _ index: Int) {
                self.row = row
                self.index = index
            }

            mutating func next() -> Element? {
                if index < 0 {
                    return nil
                }
                defer{ index = row.nexts[index]}
                return (row.cols[index], row.values[index])
            }
        }
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

