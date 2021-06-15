//
//  File.swift
//
//
//  Created by Taketo Sano on 2021/06/01.
//

import XCTest
import SwmCore
import SwmEigen
@testable import SwmMatrixTools

class EigenSparseLUFactorizationTests: XCTestCase {
    
    typealias R = RationalNumber
    typealias M<n: SizeType, m: SizeType> = R.ComputationalSparseMatrix<n, m>

    func testSolveLowerTriangular() {
        let L: M<_4, _4> = [
            1, 0, 0, 0,
            1, 1, 0, 0,
            -1,0, 1, 0,
            0, 2, 0, -1,
        ]
        let b: M<_4, _1> = [1, 2, 1, 1]
        let x = M.solveLowerTriangular(L, b)
        
        XCTAssertEqual(L * x, b)
    }
    
    func testSolveLowerTrapezoidal() {
        let L: M<_6, _4> = [
            1, 0, 0, 0,
            1, 1, 0, 0,
            -1,0, 1, 0,
            0, 2, 0, -1,
            1, -1,2, 3,
            2, 0, 0, 1
        ]
        let b: M<_6, _1> = [1, 2, 1, 1, 7, 3]
        let x = M.solveLowerTrapezoidal(L, b)
        
        if let x = x {
            XCTAssertEqual(L * x, b)
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
        let x = M.solveLowerTrapezoidal(L, b)
        XCTAssertNil(x)
    }

    func testSolveUpperTriangular() {
        let U: M<_4, _4> = [
            1, 2, 3, 4,
            0,-1, 0, 2,
            0 ,0, 1, 0,
            0, 0, 0, -1,
        ]
        let b: M<_4, _1> = [65, 3, 18, 17]
        let x = M.solveUpperTriangular(U, b)

        XCTAssertEqual(U * x, b)
    }
    
    func testSolveUpperTrapezoidal() {
        let U: M<_4, _6> = [
            1, 2, 3, 4, 5, 6,
            0,-1, 0, 2, 0, 1,
            0 ,0, 1, 0, 2, 1,
            0, 0, 0, -1,3, -1
        ]
        let b: M<_4, _1> = [65, 3, 18, 17]
        let x = M.solveUpperTrapezoidal(U, b)

        XCTAssertEqual(U * x, b)
    }
    
    func testFactorize() {
        let A: M<_5, _5> = [
            1, 0, 1, 0, 2,
            0, 0, 0, 3, 0,
            0, 0, 2, 0, 3,
            1, 0, 0, 1, 0,
            0, 0, 1, 0, 1,
        ]
        let e = A.LUfactorize()
        let (P, Q, L, U) = e.PQLU
        
        XCTAssertTrue(L.isLowerTriangular)
        XCTAssertTrue(U.isUpperTriangular)
        XCTAssertEqual(A.permute(rowsBy: P, colsBy: Q), L * U)
    }

    func testFactorizeRandom() {
        let n = 10
        let grid = (0 ..< n * n).map{ _ in R(Int.random(in: 0 ..< 5)) }
        let A: M<anySize, anySize> = .init(size: (n, n), grid: grid)
        let e = A.LUfactorize()
        let (P, Q, L, U) = e.PQLU
        
        XCTAssertTrue(L.isLowerTriangular)
        XCTAssertTrue(U.isUpperTriangular)
        XCTAssertEqual(A.permute(rowsBy: P, colsBy: Q), L * U)
    }
    
    func testSolve() {
        let A: M<_5, _5> = [
            1, 0, 1, 0, 2,
            0, 0, 0, 3, 0,
            0, 0, 2, 0, 3,
            1, 0, 0, 1, 0,
            0, 0, 1, 0, 1,
        ]
        let b: M<_5, _1> = [14,12,21,5,8]
        let e = A.LUfactorize()
        
        if let x = e.solve(b) {
            XCTAssertEqual(A * x, b)
        } else {
            XCTFail()
        }
    }

    func testSolve_noSolution() {
        let A: M<_5, _5> = [
            1, 0, 1, 0, 2,
            0, 0, 0, 3, 0,
            0, 0, 2, 0, 3,
            1, 0, 0, 1, 0,
            0, 0, 1, 0, 1,
        ]
        let b: M<_5, _1> = [14, 12, 21, 5, 7]
        let e = A.LUfactorize()
        
        if let _ = e.solve(b) {
            XCTFail()
        } else {
            // OK
        }
    }

    func testSolve2() {
        let A: M<_4, _6> = [
            1, 0, 1, 2, 2, 3,
            0, 2,-1, 3, 0,-1,
            3,-1, 2, 0, 3, 5,
            1, 3, 0, 1, 0, 3
        ]
        let b: M<_4, _1> = [40, 7, 52, 29]
        let e = A.LUfactorize()
        
        if let x = e.solve(b) {
            XCTAssertEqual(A * x, b)
        } else {
            XCTFail()
        }
    }

    func testKernel() {
        let A: M<_4, _6> = [
            1, 0, 1, 2, 2, 3,
            0, 2,-1, 3, 0,-1,
            3,-1, 2, 0, 3, 5,
            1, 3, 0, 1, 0, 3
        ]
        let e = A.LUfactorize()
        XCTAssertEqual(e.nullity, 2)
        
        let K = e.kernel
        XCTAssertEqual(K.size.cols, 2)
        XCTAssertTrue((A * K).isZero)
    }
}
