//
//  File.swift
//  
//
//  Created by Taketo Sano on 2021/06/11.
//

import Foundation
import XCTest
import SwmCore
import SwmEigen
@testable import SwmMatrixTools

class SparseLUFactorizerTests: XCTestCase {
    typealias R = RationalNumber
    typealias M<n: SizeType, m: SizeType> = MatrixIF<DefaultSparseMatrixImpl<R>, n, m>

    func testFactorizerPartialLU() {
        let A: M<_6, _9> = [
            1, 0, 1, 0, 0, 1, 1, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 3, 1, 0, 0, 0, 1, 1,
            0, 1, 1, 0, 1, 0, 0, 0, 0,
            0, 4, 0, 1, 0, 0, 2, 0, 1,
            1, 0, 1, 0, 1, 1, 0, 1, 1
        ]
        let f = SparseLUFactorizer(A)
        f.partialLU()
        
        let (P, Q, L, U) = f.result
        let S = f.target
        
        XCTAssertFalse(f.isComplete)
        XCTAssertTrue(L.isLowerTriangular)
        XCTAssertTrue(U.isUpperTriangular)
        
        let r = L.size.cols
        XCTAssertEqual(
            A.impl.permute(rowsBy: P, colsBy: Q),
            (L * U) + (.zero(size: (r, r)) ⊕ S)
        )
    }
    
    func testFactorizerPartialLU2() {
        let A: M<_5, _5> = [
            1, 0, 1, 0, 2,
            0, 0, 0, 3, 0,
            0, 0, 2, 0, 3,
            1, 0, 0, 1, 0,
            0, 0, 1, 0, 1,
        ]
        let f = SparseLUFactorizer(A)
        f.partialLU()
        
        let (P, Q, L, U) = f.result
        let S = f.target
        
        XCTAssertFalse(f.isComplete)
        XCTAssertTrue(L.isLowerTriangular)
        XCTAssertTrue(U.isUpperTriangular)
        
        let r = L.size.cols
        XCTAssertEqual(
            A.impl.permute(rowsBy: P, colsBy: Q),
            (L * U) + (.zero(size: (r, r)) ⊕ S)
        )
    }
    
    func testFactorizerFullLU() {
        let A: M<_5, _5> = [
            1, 0, 1, 0, 2,
            0, 0, 0, 3, 0,
            0, 0, 2, 0, 3,
            1, 0, 0, 1, 0,
            0, 0, 1, 0, 1,
        ]
        let f = SparseLUFactorizer(A)
        f.fullLU()
        
        let (P, Q, L, U) = f.result
        
        XCTAssertTrue(f.isComplete)
        XCTAssertTrue(L.isLowerTriangular)
        XCTAssertTrue(U.isUpperTriangular)
        XCTAssertEqual(A.impl.permute(rowsBy: P, colsBy: Q), L * U)
    }
    
    func testFactorizerFactorize() {
        let A: M<_5, _5> = [
            1, 0, 1, 0, 2,
            0, 0, 0, 3, 0,
            0, 0, 2, 0, 3,
            1, 0, 0, 1, 0,
            0, 0, 1, 0, 1,
        ]
        let f = SparseLUFactorizer(A)
        f.run()
        
        let (P, Q, L, U) = f.result
        
        XCTAssertTrue(f.isComplete)
        XCTAssertTrue(L.isLowerTriangular)
        XCTAssertTrue(U.isUpperTriangular)
        XCTAssertEqual(A.impl.permute(rowsBy: P, colsBy: Q), L * U)
    }
}
