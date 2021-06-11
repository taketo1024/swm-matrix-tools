//
//  File.swift
//  
//
//  Created by Taketo Sano on 2021/05/28.
//

import XCTest
import SwmCore
@testable import SwmMatrixTools

class MatrixPivotFinderTests: XCTestCase {
    
    typealias M<n: SizeType, m: SizeType> = Matrix<Int, n, m>
    
    func testFindPivots() {
        let A: M<_6, _9> = [
            1, 0, 1, 0, 0, 1, 1, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 1, 0, 1, 0, 0, 0, 0,
            0, 1, 0, 1, 0, 0, 1, 0, 1,
            1, 0, 1, 0, 1, 1, 0, 1, 1
        ]
        let (pivs, p, q) = A.findPivots()
        let B = A.permute(rowsBy: p, colsBy: q)
        
        XCTAssertTrue(pivs.count >= 3)
        XCTAssertTrue(
            B.submatrix(rowRange: 0 ..< pivs.count, colRange: 0 ..< pivs.count)
                .nonZeroEntries
                .allSatisfy{ (i, j, a) in i < j || (i == j && a.isIdentity) }
        )
    }
    
    func testPermuteByPivotsColBased() {
        let A: M<_6, _9> = [
            1, 0, 1, 0, 0, 1, 1, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 1, 0, 1, 0, 0, 0, 0,
            0, 1, 0, 1, 0, 0, 1, 0, 1,
            1, 0, 1, 0, 1, 1, 0, 1, 1
        ]
        let (pivs, p, q) = A.findPivots(mode: .colBased)
        let B = A.permute(rowsBy: p, colsBy: q)

        print(B.detailDescription)
        
        XCTAssertTrue(pivs.count >= 3)
        XCTAssertTrue(
            B.submatrix(rowRange: 0 ..< pivs.count, colRange: 0 ..< pivs.count)
                .nonZeroEntries
                .allSatisfy{ (i, j, a) in i > j || (i == j && a.isIdentity) }
        )
    }
}
