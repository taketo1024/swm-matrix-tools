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

public enum PivotMode {
    case rowBased, colBased
}

public final class MatrixPivotFinder {
    typealias Data = [[(Int, Double)]]
    
    public let mode: PivotMode
    public let size: MatrixSize
    private let data: Data
    private let heads: [(Int, Double)]
    private let weights: [Double]
    
    private var pivotRows: Set<Int> // Set(rows)
    private var pivotTable: [Int : Int] // [col : row]
    private var result: [(Int, Int)]

    public var debug: Bool = false
    
    internal init(mode: PivotMode = .rowBased, size: MatrixSize, data: Data, heads: [(Int, Double)], weights: [Double]) {
        self.mode = mode
        self.size = size
        self.data = data
        self.heads = heads
        self.weights = weights
        self.pivotTable = [:]
        self.pivotRows = []
        self.result = []
    }
    
    public convenience init<Impl: MatrixImpl>(_ A: Impl, mode: PivotMode = .rowBased) where Impl.BaseRing: ComputationalRing {
        let n = (mode == .rowBased) ? A.size.rows : A.size.cols
        var data: Data = Array(repeating: [], count: n)
        var heads: [(Int, Double)] = Array(repeating: (0, 0), count: n)
        var weights: [Double] = Array(repeating: 0, count: n)
        
        if mode == .rowBased {
            for (i, j, a) in A.nonZeroEntries {
                let w = a.computationalWeight
                data[i].append((j, w))
                weights[i] += w
                if heads[i].0 > j {
                    heads[i] = (j, w)
                }
            }
        } else {
            for (i, j, a) in A.nonZeroEntries {
                let w = a.computationalWeight
                data[j].append((i, w))
                weights[j] += w
                if heads[j].0 > i {
                    heads[j] = (i, w)
                }
            }
        }
        
        self.init(mode: mode, size: A.size, data: data, heads: heads, weights: weights)
    }
    
    public convenience init<Impl, n, m>(_ A: MatrixIF<Impl, n, m>, mode: PivotMode = .rowBased) where Impl.BaseRing: ComputationalRing {
        self.init(A.impl, mode: mode)
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
        Permutation(length: length, indices: order, fillRemaining: true).inverse!
    }
    
    public func run() {
        log("Start pivot search.")
        log("size: \(size)")
        
        findFLPivots()
        findFLColumnPivots()
        findCycleFreePivots()
        
        let pivots = sortPivots()
        
        if pivots.isEmpty {
            log("No pivots found.")
        } else {
            self.result = pivots
            
            log("Total: \(pivots.count)\(pivots.allSatisfy({ (i, _) in data[i].count == 1 }) ? " perfect" : "") pivots.")
        }
    }
    
    // Faugère-Lachartre pivot search
    private func findFLPivots() {
        var candidates: [Int : Int] = [:] // col -> row
        
        for i in heads.indices {
            let (j, a) = heads[i]
            if a == 1 && (!candidates.contains(key: j) || weights[i] < weights[candidates[j]!]) {
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
            data[i].map { $0.0 }
        })
        
        for i in data.indices where !pivotRows.contains(i) {
            var currentCols: Set<Int> = []
            var candidates: [(Int, Double)] = []
            
            for (j, a) in data[i] where !occupiedCols.contains(j) {
                currentCols.insert(j)
                if a == 1 {
                    candidates.append((j, a))
                }
            }
            
            if let c = candidates.min(by: { $0.1 }) {
                setPivot(i, c.0)
                occupiedCols.formUnion(currentCols)
            }
        }
        
        log("FL-col-pivots: \(pivotTable.count - count)")
    }
    
    private func findCycleFreePivots() {
        let count = pivotTable.count
        let atomic = DispatchQueue(label: "atomic", qos: .userInteractive)

        let remainingRows = data.indices.exclude{ i in pivotRows.contains(i) }
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
        let row = data[i]

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
        var queue: [Int] = []
        var queued: Set<Int> = []

        // initialize
        for (j, a) in row {
            if pivots.contains(key: j) {
                queue.append(j)
                queued.insert(j)
            } else if a == 1 {
                candidates.insert(j)
            }
        }

        while !queue.isEmpty && !candidates.isEmpty {
            let j1 = queue.removeFirst()
            let i2 = pivots[j1]!

            for (j2, _) in data[i2] {
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

        if let j = candidates.anyElement {
            return (i, j)
        } else {
            return nil
        }
    }
    
    private func sortPivots() -> [(Int, Int)] {
        let tree = PlainGraph(structure: Dictionary(keys: pivotTable.keys) { j1 -> [Int] in
            let i1 = pivotTable[j1]!
            return data[i1].reduce(into: []) { (res, next) in
                let j2 = next.0
                if j1 != j2 && pivotTable.contains(key: j2) {
                    res.append(j2)
                }
            }
        })
        let sorted = tree.topologicalSort()
        return sorted.map { v in
            let j = v.id
            return (pivotTable[j]!, j)
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
}
