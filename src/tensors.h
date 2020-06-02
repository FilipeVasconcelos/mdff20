typedef struct TENSOR_RK0 {
    double sca;
}TENSOR_RK0;
typedef struct TENSOR_RK1 {
    double a[3];
}TENSOR_RK1;
typedef struct TENSOR_RK2 {
    double ab[3][3];
}TENSOR_RK2;
typedef struct TENSOR_RK3 {
    double abc[3][3][3];
}TENSOR_RK3;
typedef struct TENSOR_RK4 {
    double abcd[3][3][3][3];
} TENSOR_RK4;
