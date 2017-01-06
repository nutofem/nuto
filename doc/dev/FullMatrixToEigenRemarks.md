@page FullMatrixToEigenRemarks

## Some common syntax changes

| `NuTo::FullMatrix` | `Eigen::Matrix` |
|:------:|:-------:|
| GetValue(i,j) / GetValue(i) | (i,j) / (i) |
| SetValue(i,j,x) / SetValue(i,x) | (i,j) = x / (i) = x |
| GetNumRows() /  GetNumColumns() | rows() / cols() |
| GetRow(i) / GetColumn(i) | row(i) / col(i)|
| SetRow(i, x) / SetColumn(i, x) | row(i) = x / col(i) = x | 
| Resize() | resize() |
| Sum() | sum()|
| Norm() | norm() |
| Abs() | cwiseAbs()|
| Max() | maxCoeff() |
| Dot() | dot() |
| matrix.ConservativeResizeRows(newNumRows) | matrix.conservativeResize(newNumRows, matrix.cols()) |

## *PITFALLS*

- `Eigen::Matrix::resize(...)` **does not** perform a `SetZero` like `NuTo::FullMatrix::Resize(...)` did. Think about it. 

Example 1)
~~~cpp
eigenVector.resize(3);
eigenVector << 1,2,3;
~~~
`setZero` not needed, since you set each value. Consider [RAII](https://en.wikipedia.org/wiki/Resource_acquisition_is_initialization)...

Example 2)
~~~cpp
eigenVector.resize(3);
for (...)
    eigenVector += someOtherVectors
~~~
`setZero` required.