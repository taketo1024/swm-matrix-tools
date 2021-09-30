//
//  File.swift
//  
//
//  Created by Taketo Sano on 2021/09/30.
//

import SwmCore

public typealias SparseMatrix<R, n, m> = MatrixIF<DefaultSparseMatrixImpl<R>, n, m> where R: Ring, n: SizeType, m: SizeType
public typealias SparseVector<R, n> = SparseMatrix<R, n, _1> where R: Ring, n: SizeType
