//
//  MatrixPivotFinder.swift
//  SwiftyMath
//
//  Created by Taketo Sano on 2019/10/26.
//

//  Implementation based on:
//
//  "Parallel Sparse PLUQ Factorization modulo p", Charles Bouillaguet, Claire Delaplace, Marie-Emilie Voge.
//  https://hal.inria.fr/hal-01646133/document
//
//  see also:
//
//  SpaSM (Sparse direct Solver Modulo p)
//  https://github.com/cbouilla/spasm

import SwmCore
import Dispatch
import TSCBasic

public enum PivotMode {
    case rowBased, colBased
}

public final class MatrixPivotFinder<R: Ring> {
    typealias Row = MatrixEliminationData<R>.Row
    
    private let data: MatrixEliminationData<R>
    private let sortedRows: [Int]   // indices of non-zero rows, sorted by weight.
    private var pivotRows: Set<Int> // Set(rows)
    private var pivotTable: [Int : Int] // [col : row]
    private var result: [(Int, Int)]

    public let mode: PivotMode
    public var debug: Bool = false
    
    internal init(data: MatrixEliminationData<R>, mode: PivotMode = .rowBased) {
        self.mode = mode
        self.data = data
        self.sortedRows = (0 ..< data.size.rows)
            .exclude{ i in data.row(i).isEmpty }
            .sorted(by: { i in data.rowWeight(i) } )
        self.pivotTable = [:]
        self.pivotRows = []
        self.result = []
        
        pivotTable.reserveCapacity(data.size.cols)
        pivotRows.reserveCapacity(data.size.rows)
    }
    
    public convenience init<Impl: MatrixImpl>(_ A: Impl, mode: PivotMode = .rowBased) where Impl.BaseRing == R {
        let data = MatrixEliminationData(A, transpose: mode == .colBased)
        self.init(data: data, mode: mode)
    }
    
    public convenience init<Impl, n, m>(_ A: MatrixIF<Impl, n, m>, mode: PivotMode = .rowBased) where Impl.BaseRing == R {
        self.init(A.impl, mode: mode)
    }
    
    public var size: MatrixSize {
        mode == .rowBased
            ? data.size
            : (data.size.cols, data.size.rows)
    }
    
    public var pivots: [(Int, Int)] {
        mode == .rowBased
            ? result
            : result.map{ (i, j) in (j, i) }
    }
    
    public var rowPermutation: Permutation<anySize> {
        asPermutation(size.rows, pivots.map{ $0.0 })
    }
    
    public var colPermutation: Permutation<anySize> {
        asPermutation(size.cols, pivots.map{ $0.1 })
    }
    
    private func asPermutation(_ length: Int, _ order: [Int]) -> Permutation<anySize> {
        .init(length: length, indices: order, fillRemaining: true).inverse!
    }
    
    public func run() {
        log("Start pivot search.")
        log("size: \(data.size), density: \(data.density)")
        
        if size.rows < 100 && size.cols < 100 {
            log("before:\n\(data.resultAs(AnySizeMatrix.self).detailDescription)")
        }
        
        findFLPivots()
        findFLColumnPivots()
        findCycleFreePivots()
        
        let pivots = sortPivots()
        
        if pivots.isEmpty {
            log("No pivots found.")
        } else if pivots.allSatisfy({ (j, i) in i == j }) {
            log("No need for permutation.")
        } else {
            self.result = pivots
            
            log("Total: \(pivots.count)\(pivots.allSatisfy({ (i, _) in data.row(i).count == 1 }) ? " perfect" : "") pivots.")
            log("")
            
            if debug && size.rows < 100 && size.cols < 100 {
                let p = asPermutation(data.size.rows, pivots.map{ $0.0 })
                let q = asPermutation(data.size.cols, pivots.map{ $0.1 })
                log("after:\n\(data.resultAs(AnySizeMatrix.self).permute(rowsBy: p, colsBy: q).detailDescription)")
            }
        }
    }
    
    // Faugère-Lachartre pivot search
    private func findFLPivots() {
        var candidates: [Int : Int] = [:] // col -> row
        
        for i in sortedRows {
            let (j, a) = data.row(i).head!.element
            if isCandidate(a) && !candidates.contains(key: j) {
                candidates[j] = i
            }
        }
        
        candidates.forEach{ (j, i) in
            setPivot(i, j)
        }
        
        log("FL-pivots: \(candidates.count)")
    }
    
    private func findFLColumnPivots() {
        let count = pivotTable.count
        var occupiedCols = Set(pivotRows.flatMap{ i in
            data.row(i).map { $0.col }
        })
        
        for i in sortedRows where !pivotRows.contains(i) {
            var currentCols: Set<Int> = []
            var candidates: [RowEntry<R>] = []
            
            for (j, a) in data.row(i) where !occupiedCols.contains(j) {
                currentCols.insert(j)
                if isCandidate(a) {
                    candidates.append((j, a))
                }
            }
            
            if let c = candidates.min(by: { data.weight(of: $0.value) }) {
                setPivot(i, c.col)
                occupiedCols.formUnion(currentCols)
            }
        }
        
        log("FL-col-pivots: \(pivotTable.count - count)")
    }
    
    private func findCycleFreePivots() {
        let count = pivotTable.count
        let atomic = DispatchQueue(label: "atomic", qos: .userInteractive)
        
        let remainingRows = sortedRows.exclude{ i in pivotRows.contains(i) }
        remainingRows.parallelForEach { i in
            while true {
                let copy = atomic.sync { self.pivotTable }
                guard let pivot = self.findCycleFreePivot(inRow: i, pivots: copy) else {
                    break
                }
                
                let done = atomic.sync { () -> Bool in
                    if copy.count != self.pivotTable.count {
                        return false
                    }
                    self.setPivot(pivot.0, pivot.1)
                    return true
                }
                
                if done {
                    break
                }
            }
        }
        
        log("cycle-free-pivots: \(pivotTable.count - count)")
    }
    
    private func findCycleFreePivot(inRow i: Int, pivots: [Int : Int]) -> (Int, Int)? {
        let row = data.row(i)
        
        //       j1
        //  i [  O     O     C  C     ]    O: queued, C: candidates
        //       |     |     ^
        //       V     |     | rmv
        // i2 [  X     o     o     O  ]
        //             |           |
        //             V           |
        // i3 [        X           |  ]
        //                         |
        
        // the following are col-indexed.
        var candidates: Set<Int> = []
        var values: [Int : R] = [:]
        var queue: [Int] = []
        var queued: Set<Int> = []
        
        // initialize
        for (j, a) in row {
            if pivots.contains(key: j) {
                queue.append(j)
                queued.insert(j)
            } else if isCandidate(a) {
                candidates.insert(j)
                values[j] = a
            }
        }
        
        while !queue.isEmpty && !candidates.isEmpty {
            let j1 = queue.removeFirst()
            let i2 = pivots[j1]!
            
            for (j2, _) in data.row(i2) {
                if pivots.contains(key: j2) && !queued.contains(j2) {
                    queue.append(j2)
                    queued.insert(j2)
                } else if candidates.contains(j2) {
                    candidates.remove(j2)
                    if candidates.isEmpty {
                        break
                    }
                }
            }
        }
        
        if let j = candidates.min(by: { j in data.weight(of: values[j]!) }) {
            return (i, j)
        } else {
            return nil
        }
    }
    
    private func sortPivots() -> [(Int, Int)] {
        let tree = Dictionary(keys: pivotTable.keys) { j1 -> [Int] in
            let i1 = pivotTable[j1]!
            return data.row(i1).map{ $0.col }.filter { j2 in
                j1 != j2 && pivotTable.contains(key: j2)
            }
        }
        let sorted = try! topologicalSort(tree.keys.toArray(), successors: { j in tree[j] ?? [] })
        return sorted.map { j in
            (pivotTable[j]!, j)
        }
    }
    
    private func setPivot(_ i: Int, _ j: Int) {
        pivotRows.insert(i)
        pivotTable[j] = i
    }
    
    private func log(_ msg: @autoclosure () -> String) {
        if debug {
            print(msg())
        }
    }
    
    private func isCandidate(_ a: R) -> Bool {
        a == .identity || a == -.identity
    }
}
