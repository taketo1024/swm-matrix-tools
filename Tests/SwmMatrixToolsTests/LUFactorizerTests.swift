//
//  File.swift
//  
//
//  Created by Taketo Sano on 2021/06/01.
//

import XCTest
import SwmCore
@testable import SwmMatrixTools

class LUFactorizationTests: XCTestCase {
    
    typealias M<n: SizeType, m: SizeType> = Matrix<Int, n, m>

    func testSolveLowerTriangular() {
        let L: M<_6, _4> = [
            1, 0, 0, 0,
            1, 1, 0, 0,
            -1,0, 1, 0,
            0, 2, 0, -1,
            1, -1,2, 3,
            2, 0, 0, 1
        ]
        let b: M<_6, _1> = [1, 2, 1, 1, 7, 3]
        let x = DefaultMatrixImpl<Int>.solveLowerTriangular(L.impl, b.impl)
        
        if let x = x {
            XCTAssertEqual(L.impl * x, b.impl)
        } else {
            XCTFail()
        }
    }
    
    func testSolveLowerTriangular_noSolution() {
        let L: M<_6, _4> = [
            1, 0, 0, 0,
            1, 1, 0, 0,
            -1,0, 1, 0,
            0, 2, 0, -1,
            1, -1,2, 3,
            2, 0, 0, 1
        ]
        let b: M<_6, _1> = [1, 2, 1, 1, 7, 1]
        let x = DefaultMatrixImpl<Int>.solveLowerTriangular(L.impl, b.impl)
        XCTAssertNil(x)
    }

    func testSolveUpperTriangular() {
        let U: M<_4, _6> = [
            1, 2, 3, 4, 5, 6,
            0,-1, 0, 2, 0, 1,
            0 ,0, 1, 0, 2, 1,
            0, 0, 0, -1,3, -1
        ]
        let b: M<_4, _1> = [65, 3, 18, 17]
        let x = DefaultMatrixImpl<Int>.solveUpperTriangular(U.impl, b.impl)

        if let x = x {
            XCTAssertEqual(U.impl * x, b.impl)
        } else {
            XCTFail()
        }
    }
    
    func testPartialLU() {
        let A: M<_6, _9> = [
            1, 0, 1, 0, 0, 1, 1, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 3, 1, 0, 0, 0, 1, 1,
            0, 1, 1, 0, 1, 0, 0, 0, 0,
            0, 4, 0, 1, 0, 0, 2, 0, 1,
            1, 0, 1, 0, 1, 1, 0, 1, 1
        ]
        let (P, Q, L, U, S) = LUFactorizer.partialLU(A.impl)
        assert(L.isLowerTriangular)
        assert(U.isUpperTriangular)
        
        print(L.detailDescription)
        print(U.detailDescription)
        print(S.detailDescription)
        print(A.impl.permute(rowsBy: P, colsBy: Q).detailDescription)
        print((L * U).detailDescription)

        let r = L.size.cols
        XCTAssertEqual(
            A.impl.permute(rowsBy: P, colsBy: Q),
            (L * U) + (.zero(size: (r, r)) ⊕ S)
        )
    }
}
