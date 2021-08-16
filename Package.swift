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
            from: "1.2.9"
//            path: "../swm-core/"
        ),
        .package(
            url: "https://github.com/taketo1024/swm-eigen.git",
            from: "1.0.0"
//            path: "../swm-eigen/"
        ),
    ],
    targets: [
        .target(
            name: "SwmMatrixTools",
            dependencies: [
                .product(name: "SwmCore", package: "swm-core"),
                .product(name: "SwmEigen", package: "swm-eigen"),
            ],
            swiftSettings: [
                .define("USE_EIGEN")
            ]
        ),
        .testTarget(
            name: "SwmMatrixToolsTests",
            dependencies: ["SwmMatrixTools"]
		),
    ]
)
