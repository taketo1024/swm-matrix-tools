// swift-tools-version:5.3
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
    name: "swm-matrix-tools",
    products: [
        .library(
            name: "SwmMatrixTools",
            targets: ["SwmMatrixTools"]),
    ],
    dependencies: [
        .package(
            url: "https://github.com/taketo1024/swm-core.git",
            from: "1.2.1"
//            path: "../swm-core/"
        ),
    ],
    targets: [
        .target(
            name: "SwmMatrixTools",
            dependencies: [
                .product(name: "SwmCore", package: "swm-core")
            ]
		),
        .testTarget(
            name: "SwmMatrixToolsTests",
            dependencies: ["SwmMatrixTools"]
		),
    ]
)
