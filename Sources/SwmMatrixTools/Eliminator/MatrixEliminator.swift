//
//  MatrixEliminator.swift
//  SwiftyMath
//
//  Created by Taketo Sano on 2017/06/09.
//  Copyright © 2017年 Taketo Sano. All rights reserved.
//

import SwmCore

extension MatrixIF where BaseRing: EuclideanRing {
    public func eliminate(form: MatrixEliminationForm = .Diagonal) -> MatrixEliminationResult<Impl, n, m> {
        let (type, transpose) = eliminatorType(form)
        let worker = !transpose
            ? MatrixEliminationWorker(self)
            : MatrixEliminationWorker(self.transposed) // TODO directly pass tranposed entries
        
        let e = type.init(worker: worker)
        e.run()
        
        return !transpose
            ? e.result(as: MatrixEliminationResult.self)
            : e.result(as: MatrixEliminationResult.self).transposed
    }
    
    private func eliminatorType(_ form: MatrixEliminationForm) -> (MatrixEliminator<BaseRing>.Type, Bool) {
        switch form {
        case .RowEchelon:
            return (RowEchelonEliminator.self, false)
        case .ColEchelon:
            return (RowEchelonEliminator.self, true)
        case .RowHermite:
            return (ReducedRowEchelonEliminator.self, false)
        case .ColHermite:
            return (ReducedRowEchelonEliminator.self, true)
        case .Diagonal:
            return (DiagonalEliminator.self, false)
        case .Smith:
            return (SmithEliminator.self, false)
        default:
            return (MatrixEliminator.self, false)
        }
    }
}

public class MatrixEliminator<R: Ring> {
    let worker: MatrixEliminationWorker<R>
    var debug: Bool = false
    var aborted: Bool = false
    
    private var rowOps: [RowElementaryOperation<R>]
    private var colOps: [ColElementaryOperation<R>]
    
    required init(worker: MatrixEliminationWorker<R>) {
        self.worker = worker
        self.rowOps = []
        self.colOps = []
    }
    
    public final func run() {
        log("Start: \(self)")
        
        prepare()
        
        var itr = 0
        while !aborted && !isDone() {
            log("\(self) iteration: \(itr)")
            
            if debug {
                printCurrentMatrix()
            }
            
            iteration()
            
            log("")
            itr += 1
        }
        
        finalize()
        
        log("Done:  \(self), \(rowOps.count + colOps.count) steps")
        
        if debug {
            printCurrentMatrix()
        }
    }
    
    public func result<Impl, n, m>(as: MatrixEliminationResult<Impl, n, m>.Type) -> MatrixEliminationResult<Impl, n, m> where Impl.BaseRing == R {
        .init(
            form: form,
            size: worker.size,
            entries: worker.entries,
            headEntries: worker.headEntries,
            rowOps: rowOps,
            colOps: colOps
        )
    }
    
    public var description: String {
        "\(type(of: self))"
    }
    
    // MARK: Internal methods

    final func subrun(_ type: MatrixEliminator<R>.Type, transpose: Bool = false) {
        let e = type.init(worker: worker)
        if transpose {
            e.transpose()
        }
        
        e.run()
        
        if transpose {
            e.transpose()
        }
        
        rowOps += e.rowOps
        colOps += e.colOps
    }
    
    final func abort() {
        aborted = true
    }
    
    final func apply(_ s: RowElementaryOperation<R>) {
        worker.apply(s)
        append(s)
    }

    final func append(_ s: RowElementaryOperation<R>) {
        rowOps.append(s)
        log("\(s)")
    }
    
    final func append(_ s: ColElementaryOperation<R>) {
        colOps.append(s)
        log("\(s)")
    }
    
    final func append(_ s: [RowElementaryOperation<R>]) {
        rowOps.append(contentsOf: s)
        log(s.map{ "\($0)"}.joined(separator: "\n"))
    }
    
    final func transpose() {
        log("Transpose: \(self)")
        worker.transpose()
        (rowOps, colOps) = (colOps.map{ $0.transposed }, rowOps.map{ $0.transposed })
    }
    
    final func log(_ msg: @autoclosure () -> String) {
        if debug {
            print(msg())
        }
    }
    
    final func printCurrentMatrix() {
        let (size, entries) = (worker.size, worker.entries)
        if size.rows > 100 || size.cols > 100 {
            return
        }
        
        print("\n", AnySizeMatrix(size: size, entries: entries).detailDescription, "\n")
    }

    // MARK: Methods to be overridden
    
    public var form: MatrixEliminationForm { .none }
    func prepare() {}
    func isDone() -> Bool { true }
    func iteration() {}
    func finalize() {}
}
